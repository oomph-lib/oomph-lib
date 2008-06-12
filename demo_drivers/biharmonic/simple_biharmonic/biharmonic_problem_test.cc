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

// This will move to meshes
#include "hermite_element_quad_mesh.h"
#include "one_d_hermite_element_mesh.h"  
#include "topologically_rectangular_domain.h"


// use oomph namespace
using namespace oomph;






//=============================================================================
// Non Uniform Mesh Spacing Namespace - for a topologically rectangular 
// domain. In macro element coordinates the right most (or upper most) element
// will be M_i times as big as the left most (or lower) element. The elements
// in between scale linearly.
// Note - the number of elements in the each direction is required - the 
// n_element must be set prior to use
//=============================================================================
namespace NonUniformSpacing
{
 // Number of elements in the x0 and x1 direction
 Vector<unsigned> n_element(2);
 
 // The multiplication factor in the x0 and x1 direction
 Vector<double> M(2);
 
 // The mesh spacing functions - takes the position of a uniformly spaced node
 // and returns the position of the non-uniformly spaced node (both in macro
 // element coordinates)
 void spacing(const Vector<double>& s_uniform, Vector<double>& s_non_uniform)
 {  
  for (unsigned i = 0; i < 2; i++)
   {
    double x0 = 2 / (double(n_element[i])*(1+0.5*(M[i]-1)));
    double s = (M[i]-1)*x0/double(n_element[i]-1);
    double n = 0.5*(s_uniform[i]+1)*double(n_element[i]);
    s_non_uniform[i] = n*x0+0.5*s*n*(n-1)-1;
   }
 }
}




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////





//=============================================================================
// TWO DIMENSIONAL BIHARMONIC TEST 1 (PLATE PROBLEM)
// u=x^4*y^4
//=============================================================================





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
// Two Dimensional Biharmonic Test Problem 1a (plate)
// All Edges Clamped
//=============================================================================
class BiharmonicTestProblem1a : public BiharmonicPlateProblem<2>
{

private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

public:

 // constructor
 BiharmonicTestProblem1a(const unsigned n_x, const unsigned n_y)
  {
   // use test functions 1 namespace
   using namespace BiharmonicTestFunctions1;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(1,1);

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<BiharmonicElement<2> >(n_x, n_y, Domain_pt);

   // set the bulk element source function
   impose_surface_load(surface_load);

   // clamp edge on all boundaries
   impose_clamped_edge(0);
   impose_clamped_edge(1,u_NE,dudn_NE);
   impose_clamped_edge(2,u_NE,dudn_NE);
   impose_clamped_edge(3);

   // assign equation numbers
   assign_eqn_numbers();
  }

 /// Destructor - just deletes domain pt
 virtual ~BiharmonicTestProblem1a()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};



//=============================================================================
// Two Dimensional Biharmonic Test Problem 1a
// All Edges Simply Supported
//=============================================================================
class BiharmonicTestProblem1b : public BiharmonicPlateProblem<2>
{

 private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

 public:

 // constructor
 BiharmonicTestProblem1b(const unsigned n_x, const unsigned n_y)
  {

   // use test functions 1 namespace
   using namespace BiharmonicTestFunctions1;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(1,1);

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<BiharmonicElement<2> >(n_x, n_y, Domain_pt);

   // set the bulk element source function
   impose_surface_load(surface_load);

   // all edges simply supported
   impose_simply_supported_edge(0);
   impose_simply_supported_edge(1,u_NE,flux0_NE);
   impose_simply_supported_edge(2,u_NE,flux0_NE);
   impose_simply_supported_edge(3);

   // assign equation numbers
   assign_eqn_numbers();
  }

 /// Destructor just deletes the domain
 ~BiharmonicTestProblem1b()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};



//=============================================================================
// Two Dimensional Biharmonic Test Problem 1c (plate)
// S & W Edges Clamped
// N & E Edges Free
//=============================================================================
class BiharmonicTestProblem1c : public BiharmonicPlateProblem<2>
{

 private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

 public:

 // constructor
 BiharmonicTestProblem1c(const unsigned n_x, const unsigned n_y)
  {   
   
   // use test functions 1 namespace
   using namespace BiharmonicTestFunctions1;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(1,1);

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<BiharmonicElement<2> >(n_x, n_y, Domain_pt);

   // set the bulk element source function
   impose_surface_load(surface_load);

   // impose boundary conditions
   impose_clamped_edge(0);
   impose_free_edge(1,flux0_NE,flux1_NE);
   impose_free_edge(2,flux0_NE,flux1_NE); 
   impose_clamped_edge(3);

   // assign equation numbers
   assign_eqn_numbers();
  }


 /// Destructor just deletes the domain
 ~BiharmonicTestProblem1c()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};



//=============================================================================
// Two Dimensional Biharmonic Test Problem 1d (plate)
// S & W Edges Clamped
// E Edge Free
// N Edge Simply Supported
// Non Uniform Mesh Spacing
//=============================================================================
class BiharmonicTestProblem1d : public BiharmonicPlateProblem<2>
{
 
 private:
 
  // Domain pointer
  TopologicallyRectangularDomain* Domain_pt;

 public:

 // constructor
 BiharmonicTestProblem1d(const unsigned n_x, const unsigned n_y)
  {   
   
   // use test functions 1 namespace
   using namespace BiharmonicTestFunctions1;

   // set problem parameters in non-uniform mesh spacing namespace
   NonUniformSpacing::n_element[0] = n_x;
   NonUniformSpacing::n_element[1] = n_y;
   NonUniformSpacing::M[0] = 0.25;
   NonUniformSpacing::M[1] = 0.25;  

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(1,1);

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<BiharmonicElement<2> >(n_x, n_y, Domain_pt,
                                           NonUniformSpacing::spacing); 
   
   // set the bulk element source function
   impose_surface_load(surface_load);

   // impose boundary conditions
   impose_clamped_edge(0);
   impose_free_edge(1,flux0_NE,flux1_NE);
   impose_simply_supported_edge(2,u_NE,flux0_NE);
   impose_clamped_edge(3);
   
   // assign equation numbers
   assign_eqn_numbers();
  }


 /// Destructor - just deletes the domain
 ~BiharmonicTestProblem1d()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};





///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////





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


 // DIRICHLET BCs - u

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


 // DIRICHLET BCs - dudn

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


 // NEUMANN BOUNDARY CONDITIONS - laplacian(u)

 double laplacian_u(const Vector<double> x)
 {
  return(-4*sin(x[0]*x[0]+x[1]*x[1])*x[0]*x[1]
         +4*cos(x[0]*x[0]+x[1]*x[1])*x[1]/x[0]
         +2*sin(x[0]*x[0]+x[1]*x[1])*x[1]/(x[0]*x[0]*x[0])
         -4*sin(x[0]*x[0]+x[1]*x[1])*x[1]*x[1]*x[1]/x[0]);
 }
 void flux0_N(const double& s, double& laplacian)
 {
  Vector<double> x(2);
  boundary_N(s,x);
  laplacian = laplacian_u(x);
 }
 void flux0_E(const double& s, double& laplacian)
 {      
  Vector<double> x(2);
  boundary_E(s,x);
  laplacian = laplacian_u(x);
 }
 void flux0_S(const double& s, double& laplacian)
 {  
  laplacian = 0;
 }
 void flux0_W(const double& s, double& laplacian)
 {      
  Vector<double> x(2);
  boundary_W(s,x);
  laplacian = laplacian_u(x);
 }


 // NEUMANN BOUNDARY CONDITIONS - dlaplacian(u)/dn
 
 double dlaplacianudx_0(const Vector<double> x)
 {
  return(-8*cos(x[0]*x[0]+x[1]*x[1])*x[0]*x[0]*x[1]
         -12*sin(x[0]*x[0]+x[1]*x[1])*x[1]
         -6*sin(x[0]*x[0]+x[1]*x[1])*x[1]/(x[0]*x[0]*x[0]*x[0])
         -8*cos(x[0]*x[0]+x[1]*x[1])*(x[1]*x[1]*x[1])
         +4*sin(x[0]*x[0]+x[1]*x[1])*x[1]*x[1]*x[1]/(x[0]*x[0]));
 }
 double dlaplacianudx_1(const Vector<double> x)
 {
  return(-8*cos(x[0]*x[0]+x[1]*x[1])*x[1]*x[1]*x[0]
         -4*sin(x[0]*x[0]+x[1]*x[1])*x[0]
         -20*sin(x[0]*x[0]+x[1]*x[1])*x[1]*x[1]/x[0]
         +4*cos(x[0]*x[0]+x[1]*x[1])/x[0]
         +4*cos(x[0]*x[0]+x[1]*x[1])*x[1]*x[1]/(x[0]*x[0]*x[0])
         +2*sin(x[0]*x[0]+x[1]*x[1])/(x[0]*x[0]*x[0])
         -8*cos(x[0]*x[0]+x[1]*x[1])*x[1]*x[1]*x[1]*x[1]/x[0]);
 }
 void flux1_N(const double& s, double& du)
 {   
  Vector<double> x(2);
  boundary_N(s,x);
  Vector<double> n(2);
  normal_N(s,n);
  du = dlaplacianudx_0(x)*n[0] + dlaplacianudx_1(x)*n[1];
 }
 void flux1_E(const double& s, double& du)
 {   
  Vector<double> x(2);
  boundary_E(s,x);
  Vector<double> n(2);
  normal_E(s,n);
  du = dlaplacianudx_0(x)*n[0] + dlaplacianudx_1(x)*n[1];
 }
 void flux1_S(const double& s, double& du)
 {   
  Vector<double> x(2);
  boundary_S(s,x);
  Vector<double> n(2);
  normal_S(s,n);
  du = dlaplacianudx_0(x)*n[0] + dlaplacianudx_1(x)*n[1];
 }
 void flux1_W(const double& s, double& du)
 {   
  Vector<double> x(2);
  boundary_W(s,x);
  Vector<double> n(2);
  normal_W(s,n);
  du = dlaplacianudx_0(x)*n[0] + dlaplacianudx_1(x)*n[1];
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
// Two Dimensional Biharmonic Test Problem 2a (plate)
// All Edges Clamped
//=============================================================================
class BiharmonicTestProblem2a : public BiharmonicPlateProblem<2>
{

 private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

 public:

 // constructor
 BiharmonicTestProblem2a(const unsigned n_x, const unsigned n_y)
  {
   
   // namespace from problem functions
   using namespace BiharmonicTestFunctions2;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(boundary_N,boundary_E,
                                                  boundary_S,boundary_W);

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<BiharmonicElement<2> >(n_x, n_y, Domain_pt);
   
   // set the bulk element source function
   impose_surface_load(surface_load);

   // clamp edge on all boundaries
   impose_clamped_edge(0,u_S,dudn_S);
   impose_clamped_edge(1,u_E,dudn_E);
   impose_clamped_edge(2,u_N,dudn_N);
   impose_clamped_edge(3,u_W,dudn_W);

   // assign equation numbers
   assign_eqn_numbers();
  }

 /// Destructor - just deletes the domain, all other cleanup done in base class
 ~BiharmonicTestProblem2a()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};




//=============================================================================
// Two Dimensional Biharmonic Test Problem 2b (plate)
// All Edges Simply Supported
//=============================================================================
class BiharmonicTestProblem2b : public BiharmonicPlateProblem<2>
{
 
private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

 public:

 // constructor
 BiharmonicTestProblem2b(const unsigned n_x, const unsigned n_y)
  {
   
   // namespace for problem functions
   using namespace BiharmonicTestFunctions2;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(boundary_N,boundary_E,
                                                  boundary_S,boundary_W);

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<BiharmonicElement<2> >(n_x, n_y, Domain_pt);
   
   // set the bulk element source function
   impose_surface_load(surface_load);

   // simply support all edges
   impose_simply_supported_edge(0,u_S,flux0_S);
   impose_simply_supported_edge(1,u_E,flux0_E);
   impose_simply_supported_edge(2,u_N,flux0_N);
   impose_simply_supported_edge(3,u_W,flux0_W);   

   // assign equation numbers
   assign_eqn_numbers();
  }

 /// Destructor - just deletes the domain, all other cleanup done in bass class
 ~BiharmonicTestProblem2b()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};




//=============================================================================
// Two Dimensional Biharmonic Test Problem 2c (plate)
// S & E Edges Clamped
// N & W Edges Free
//=============================================================================
class BiharmonicTestProblem2c : public BiharmonicPlateProblem<2>
{
 
private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

public:

 // constructor
 BiharmonicTestProblem2c(const unsigned n_x, const unsigned n_y)
  {
   
   // namespace for problem functions
   using namespace BiharmonicTestFunctions2;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(boundary_N,boundary_E,
                                                  boundary_S,boundary_W);

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<BiharmonicElement<2> >(n_x, n_y, Domain_pt);

   // set the bulk element source function
   impose_surface_load(surface_load);

   // impose boundary conditions
   impose_free_edge(0,flux0_S,flux1_S);
   impose_free_edge(1,flux0_E,flux1_E);
   impose_clamped_edge(2,u_N,dudn_N);
   impose_clamped_edge(3,u_W,dudn_W);

   // assign equation numbers
   assign_eqn_numbers();
  }

 /// Destructor - just deletes the domain, all other cleanup done in base class
 ~BiharmonicTestProblem2c()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};




//=============================================================================
// Two Dimensional Biharmonic Test Problem 2d (plate)
// E & W edges simply supported
// S edge free
// N edge clamped
// Non-Uniform Spacing
//=============================================================================
class BiharmonicTestProblem2d : public BiharmonicPlateProblem<2>
{
 
private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

public:

 // constructor
 BiharmonicTestProblem2d(const unsigned n_x, const unsigned n_y)
  {
   
   // namespace for problem functions
   using namespace BiharmonicTestFunctions2;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(boundary_N,boundary_E,
                                                  boundary_S,boundary_W);

   // Non uniform mesh spacing
   NonUniformSpacing::n_element[0] = n_x;
   NonUniformSpacing::n_element[1] = n_y;
   NonUniformSpacing::M[0] = 1;
   NonUniformSpacing::M[1] = 0.5;  

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<BiharmonicElement<2> >(n_x, n_y, Domain_pt,
                                           NonUniformSpacing::spacing); 
   
   // set the bulk element source function
   impose_surface_load(surface_load);
   
   // impose boundary conditions
   impose_simply_supported_edge(1,u_E,flux0_E);
   impose_simply_supported_edge(3,u_W,flux0_W);
   impose_free_edge(0,flux0_S,flux1_S);
   impose_clamped_edge(2,u_N,dudn_N);

   // assign equation numbers
   assign_eqn_numbers();
  }


 /// Destructor - just deletes the domain, all other cleanup done in base class
 ~BiharmonicTestProblem2d()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};






///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////





//=============================================================================
// TWO DIMENSIONAL BIHARMONIC TEST 3 (FLUID PROBLEM)
// Poiseuille Flow in Rectangular Domain
//=============================================================================




//=============================================================================
// two dimensional biharmonic problem 3 namespace - contains all the problem
// functions
//=============================================================================
namespace BiharmonicTestFunctions3
{

 // FLUID VELOCITY FUNCTIONS

 void v_W(const double& s, Vector<double>& u)
 {
  double x = 0.5*(s+1);
  u[0] = -(6*x - 6*x*x);
  u[1] = 0.0;
 } 
 void v_E(const double& s, Vector<double>& u)
 {
  double x = 0.5*(s+1);
  u[0] = (6*x - 6*x*x);
  u[1] = 0.0;
 } 

 // SOLUTIONS

 void solution(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = (3*x[1]*x[1] - 2*x[1]*x[1]*x[1]);
 }
}



//=============================================================================
// Two Dimensional Biharmonic Test Problem 3a
// Rectangular Domain
// N & S edges - solid boundary
// E & W edges - imposed flow
//=============================================================================
class BiharmonicTestProblem3a : public BiharmonicFluidProblem<2>
{
  
private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

public:

 // constructor
 BiharmonicTestProblem3a(const unsigned n_x, const unsigned n_y)
  {
   
   // use test functions 3 namespace
   using namespace BiharmonicTestFunctions3;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(3,1);

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<BiharmonicElement<2> >(n_x, n_y, Domain_pt); 
   
   // impose boundary conditions
   impose_solid_boundary_on_edge(0);
   impose_fluid_flow_on_edge(1,v_E);
   impose_solid_boundary_on_edge(2,1.0);
   impose_fluid_flow_on_edge(3,v_W);

   // assign equation numbers
   assign_eqn_numbers();
  }

 /// Destructor - just deletes the domain, all other cleanup done in base class
 ~BiharmonicTestProblem3a()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};




//=============================================================================
// Two Dimensional Biharmonic Test Problem 3b
// Rectangular Domain
// N & S edges - solid boundary
// W edge - imposed flow
// E edge - traction free
//=============================================================================
class BiharmonicTestProblem3b : public BiharmonicFluidProblem<2>
{
  
private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

public:

 // constructor
 BiharmonicTestProblem3b(const unsigned n_x, const unsigned n_y)
  {
   
   // use test functions 3 namespace
   using namespace BiharmonicTestFunctions3;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(3,1);

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<BiharmonicElement<2> >(n_x, n_y, Domain_pt); 
   
   // impose boundary conditions
   impose_solid_boundary_on_edge(0);
   impose_traction_free_edge(1);
   impose_solid_boundary_on_edge(2,1.0);
   impose_fluid_flow_on_edge(3,v_W);

   // assign equation numbers
   assign_eqn_numbers();
  }

 /// Destructor - just deletes the domain, all other cleanup done in base class
 ~BiharmonicTestProblem3b()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};





///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////





//=============================================================================
// TWO DIMENSIONAL BIHARMONIC TEST 4 (FLUID PROBLEM)
// Couette Flow in Rectangular Domain
//=============================================================================




//=============================================================================
// two dimensional biharmonic problem 4 namespace - contains all the problem
// functions
//=============================================================================
namespace BiharmonicTestFunctions4
{

 // FLUID VELOCITY FUNCTIONS

 void v_W(const double& s, Vector<double>& u)
 {
  double x = 0.5*(s+1);
  u[0] = -2*x;
  u[1] = 0.0;
 } 
 void v_E(const double& s, Vector<double>& u)
 {
  double x = 0.5*(s+1);
  u[0] = 2*x;
  u[1] = 0.0;
 } 
 void v_N(const double& s, Vector<double>& u)
 {
  u[0] = 0.0;
  u[1] = 2.0;
 } 

 // SOLUTIONS

 void solution(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = x[1]*x[1];
 }
}




//=============================================================================
// Two Dimensional Biharmonic Test Problem 4a
// Rectangular Domain
// S edge - solid boundary
// N, E & W edges - imposed flow
//=============================================================================
class BiharmonicTestProblem4a : public BiharmonicFluidProblem<2>
{
  
private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

public:

 // constructor
 BiharmonicTestProblem4a(const unsigned n_x, const unsigned n_y)
  {
   
   // use test functions 4 namespace
   using namespace BiharmonicTestFunctions4;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(3,1);

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<BiharmonicElement<2> >(n_x, n_y, Domain_pt); 
   
   // impose boundary conditions
   impose_solid_boundary_on_edge(0);
   impose_fluid_flow_on_edge(1,v_E);   
   impose_fluid_flow_on_edge(2,v_N);
   impose_fluid_flow_on_edge(3,v_W);
   // assign equation numbers
   assign_eqn_numbers();
  }

 /// Destructor - just deletes the domain, all other cleanup done in base class
 ~BiharmonicTestProblem4a()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};




//=============================================================================
// Two Dimensional Biharmonic Test Problem 4b
// Rectangular Domain
// S edge - solid boundary
// N & W edges - imposed flow
// E edge - traction free
// Non Uniform Spacing
//=============================================================================
class BiharmonicTestProblem4b : public BiharmonicFluidProblem<2>
{
  
private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

public:

 // constructor
 BiharmonicTestProblem4b(const unsigned n_x, const unsigned n_y)
  {
   
   // use test functions 4 namespace
   using namespace BiharmonicTestFunctions4;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(3,1);

   // Non uniform mesh spacing
   NonUniformSpacing::n_element[0] = n_x;
   NonUniformSpacing::n_element[1] = n_y;
   NonUniformSpacing::M[0] = 1;
   NonUniformSpacing::M[1] = 0.25;  

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<BiharmonicElement<2> >(n_x, n_y, Domain_pt,
                                           NonUniformSpacing::spacing); 

   
   // impose boundary conditions
   impose_solid_boundary_on_edge(0);
   impose_traction_free_edge(1);
   impose_fluid_flow_on_edge(2,v_N);
   impose_fluid_flow_on_edge(3,v_W);
   // assign equation numbers
   assign_eqn_numbers();
  }

 /// Destructor - just deletes the domain, all other cleanup done in base class
 ~BiharmonicTestProblem4b()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};





///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////





//=============================================================================
// TWO DIMENSIONAL BIHARMONIC TEST 5 (FLUID PROBLEM)
// Couette Flow in a Curved Domain
//=============================================================================



// two dimensional biharmonic fluid problem 5 namespace - contains all the 
// problem functions
namespace BiharmonicTestFunctions5
{


 // PARAMETERS

 double Pi = MathematicalConstants::Pi;
 double theta = Pi/4;
 double r_min = 1;
 double r_max = 2;


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


 // FLUID VELOCITY FUNCTIONS

 void v_N(const double& s, Vector<double>& u)
 {
  double r = r_min + 0.5*(s+1)*(r_max-r_min);
  u[0] = 1/r-r;
  u[1] = 0.0;
 } 
 void v_S(const double& s, Vector<double>& u)
 {
  double r = r_min + 0.5*(s+1)*(r_max-r_min);
  u[0] = r-1/r;
  u[1] = 0.0;
 } 
 void v_E(const double& s, Vector<double>& u)
 {
  u[0] = 0.0;
  u[1] = -1.5;
 } 


 // SOLUTIONS

 void solution(const Vector<double>& x, Vector<double>& u)
 {
  double r = sqrt(x[0]*x[0] + x[1]*x[1]);
  u[0] = 0.5*r*r - log(r);
 }
}




//=============================================================================
// Two Dimensional Biharmonic Test Problem 5a
// Rectangular Domain
// W edge - solid boundary
// N, E & S edges - imposed flow
//=============================================================================
class BiharmonicTestProblem5a : public BiharmonicFluidProblem<2>
{
  
private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

public:

 // constructor
 BiharmonicTestProblem5a(const unsigned n_x, const unsigned n_y)
  {
   
   // use test functions 5 namespace
   using namespace BiharmonicTestFunctions5;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(boundary_N,boundary_E,
                                                  boundary_S,boundary_W);

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<BiharmonicElement<2> >(n_x, n_y, Domain_pt); 
   
   // impose boundary conditions
   impose_fluid_flow_on_edge(0,v_S);
   impose_fluid_flow_on_edge(1,v_E);
   impose_fluid_flow_on_edge(2,v_N);
   impose_solid_boundary_on_edge(3,0.5);

   // assign equation numbers
   assign_eqn_numbers();
  }

 /// Destructor - just deletes the domain, all other cleanup done in base class
 ~BiharmonicTestProblem5a()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};




//=============================================================================
// Two Dimensional Biharmonic Test Problem 5b
// Rectangular Domain
// N & S edges - solid boundary
// W edge - imposed flow
// E edge - traction free
// Non Uniform Mesh
//=============================================================================
class BiharmonicTestProblem5b : public BiharmonicFluidProblem<2>
{
  
private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

public:

 // constructor
 BiharmonicTestProblem5b(const unsigned n_x, const unsigned n_y)
  {
      
   // use test functions 5 namespace
   using namespace BiharmonicTestFunctions5;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(boundary_N,boundary_E,
                                                  boundary_S,boundary_W);

   // Non uniform mesh spacing
   NonUniformSpacing::n_element[0] = n_x;
   NonUniformSpacing::n_element[1] = n_y;
   NonUniformSpacing::M[0] = 1;
   NonUniformSpacing::M[1] = 0.5;  

   // assemble mesh	
   mesh_pt() = new 
    HermiteQuadMesh<Hijacked<BiharmonicElement<2> > >(n_x, n_y, Domain_pt,
                                           NonUniformSpacing::spacing); 

   
   // impose boundary conditions
   impose_traction_free_edge(0);
   impose_fluid_flow_on_edge(1,v_E);
   impose_fluid_flow_on_edge(2,v_N);
   impose_solid_boundary_on_edge(3,0.5);

   // assign equation numbers
   assign_eqn_numbers();

  }

 /// Destructor - just deletes the domain, all other cleanup done in base class
 ~BiharmonicTestProblem5b()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};




//=============================================================================
/// main
//=============================================================================
int main()
{

 // number of elements in each direction
 unsigned n_element;

 // Set doc_info
 DocInfo doc_info;
 doc_info.set_directory("RESULTS");


 // BIHARMONIC TEST PROBLEM 1a (PLATE)
 // ==================================
 {
  n_element = 30;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 1a (Plate)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem1a problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "1a";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions1::solution);
 }


 // BIHARMONIC TEST PROBLEM 1b (PLATE)
 // ==================================
 {  
  n_element = 20;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 1b (Plate)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem1b problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "1b";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions1::solution);
 }


 // BIHARMONIC TEST PROBLEM 1c (PLATE)
 // ==================================
 {
  n_element = 30;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 1c (Plate)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem1c problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "1c";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions1::solution);
 }


 // BIHARMONIC TEST PROBLEM 1d (PLATE)
 // ==================================
 { 
  n_element = 80;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 1d (Plate)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem1d problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "1d";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions1::solution);
 }


 // BIHARMONIC TEST PROBLEM 2a (PLATE)
 // ==================================
 {
  n_element = 10;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 2a (Plate)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem2a problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "2a";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions2::solution);
 }


 // BIHARMONIC TEST PROBLEM 2b (PLATE)
 // ==================================
 {  
  n_element = 20;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 2b (Plate)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem2b problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "2b";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions2::solution);
 }


 // BIHARMONIC TEST PROBLEM 2c (PLATE)
 // ==================================
 {
  n_element = 30;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 2c (Plate)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem2c problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "2c";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions2::solution);
 }


 // BIHARMONIC TEST PROBLEM 2d (PLATE)
 // ==================================
 { 
  n_element = 60;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 2d (Plate)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem2d problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "2d";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions2::solution);
 }


 // BIHARMONIC TEST PROBLEM 3a (FLUID)
 // ==================================
 { 
  n_element = 5;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 3a (Fluid)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem3a problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "3a";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions3::solution);
 }


 // BIHARMONIC TEST PROBLEM 3b (FLUID)
 // ==================================
 { 
  n_element = 10;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 3b (Fluid)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem3b problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "3b";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions3::solution);
 }


 // BIHARMONIC TEST PROBLEM 4a (FLUID)
 // ==================================
 { 
  n_element = 5;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 4a (Fluid)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem4a problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "4a";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions4::solution);
 }


 // BIHARMONIC TEST PROBLEM 4b (FLUID)
 // ==================================
 { 
  n_element = 40;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 4b (Fluid)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem4b problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "4b";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions4::solution);
 }


 // BIHARMONIC TEST PROBLEM 5a (FLUID)
 // ==================================
 { 
  n_element = 25;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 5a (Fluid)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem5a problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "5a";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions5::solution);
 }


 // BIHARMONIC TEST PROBLEM 5b (FLUID)
 // ==================================
 { 
  n_element = 25;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 5b (Fluid)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem5b problem(n_element,n_element); 
  problem.newton_solve();
  doc_info.label() = "5b";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions5::solution);
 }

/*
 for (unsigned i = 0; i < 4; i++)
  {

 // BIHARMONIC TEST PROBLEM 5b (FLUID)
 // ==================================
 { 
  n_element = 30;
  std::cout << "\n\n\nBIHARMONIC TEST PROBLEM 5b (Fluid)\n"
            << "(" << n_element << " x " << n_element << " elements)\n"
            << "==================================\n\n";
  BiharmonicTestProblem5b problem(n_element,n_element); 
//  CCDoubleMatrix jacobian;
//  Vector<double> residual;
//  problem.get_jacobian(residual,jacobian);
//  std::ofstream some_file1;
//  char filename1[100];
//  sprintf(filename1, "J_%i.dat",i);
//  some_file1.open(filename1);
//  jacobian.sparse_indexed_output(some_file1);
//  some_file1.close();

  // create a point to a hijacked biharmonic element
  Hijacked<BiharmonicElement<2> >* hijacked_element_pt;
  
  
  for (unsigned e = 0; e < 3; e++)
   {
    std::cout << "e = " << e << std::endl;
//", element_pt(" << e << ") = " 
//              << problem.mesh_pt()->element_pt(e) << std::endl;
    hijacked_element_pt =
     dynamic_cast<Hijacked<BiharmonicElement<2> > *>
     (problem.mesh_pt()->element_pt(e));

    for (unsigned n = 0; n < 2; n++)
     {
      std::cout << " node_pt(" << n << ") = "
                << hijacked_element_pt->node_pt(n) << std::endl;
      for (unsigned k = 0; k < 4; k++)
       {
        std::cout << "   k = " << k << " -      " 
                  << hijacked_element_pt->is_nodal_value_hijacked(n,k)
                  << "   " << hijacked_element_pt->node_pt(n)->eqn_number_pt(k)
                  << "   " << hijacked_element_pt->node_pt(n)->eqn_number(k)
                  << std::endl;
       }
     }
   }
  for (unsigned n = 0; n < 9; n++)
   {
    std::cout << "node_pt(" << n << ") = " << problem.mesh_pt()->node_pt(n) 
              << std::endl;
   }
  for (unsigned e = 0; e < 9; e++)
   {
    std::cout << "element_pt(" << e << ") = " 
              << problem.mesh_pt()->finite_element_pt(e) 
              << std::endl;
   }
  for (unsigned b = 0; b < 4; b++)
   {
    for (unsigned e = 0; e < 3; e++)
     {
      std::cout << "boundary_element_pt(" << b << "," << e << ") = " 
              << problem.mesh_pt()->boundary_element_pt(b,e) 
              << std::endl;
     }
   }  

  problem.newton_solve();
  doc_info.label() = "5b";  
  problem.doc_solution(doc_info,BiharmonicTestFunctions5::solution);
 }
  }
*/
}
