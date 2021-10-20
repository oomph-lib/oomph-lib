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
//Driver for a simple 2D poisson problem

//Generic routines
#include "generic.h"

// The Poisson equation
#include "poisson.h"

// The mesh
#include "meshes/simple_rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

// The shiny new Poisson with singularity elements
#include "poisson_elements_with_singularity.h"


//===== start_of_namespace=============================================
/// Namespace 
//=====================================================================
namespace Global_parameters
{

 /// element multiplier for convergence tess
 unsigned Element_multiplier=1;
 
 /// Boolean that imposes the blending or not
 bool Blend = false; 
 
 /// Limit of the blending region
 double R_blend = 0.5;

 // Pointer to the direction at which the derivative used in the residuals
 // of C1 and C2 will be computed
 unsigned Direction = 1;

 /// Blending function based on the rectangular function but made smoother
 double b(const double& r)
 {
  using namespace MathematicalConstants;
  if (Blend)
   {
    if (r<R_blend)
     {
      return 0.5+0.5*cos(Pi*r/R_blend);
     }
    else
     {
      return 0.0;
     }
   }
  else
   {
    return 1.0;
   }
 }
 
 /// Cartesian coordinates centered at the point (0.5,1)
 Vector<double> x1(const Vector<double>& coord)
 {
  Vector<double> new_x(2);
  new_x[0] = 0.5-coord[0];
  if (coord[1]>1.0)
   {
    new_x[1] = 0.0;
   }
  else
   {
    new_x[1] = 1.0-coord[1];
   }
  return new_x;
 }

 /// Cartesian coordinates centered at the point (1.5,1)
 Vector<double> x2(const Vector<double>& coord)
 {
  Vector<double> new_x(2);
  new_x[0] = 1.5-coord[0];
  new_x[1] = 1.0-coord[1];
  return new_x;
 }

 /// Polar coordinates (r,phi) centered at the point x
 Vector<double> polar(const Vector<double>& coord)
 {
  Vector<double> polar_coord(2);
  
  // r
  polar_coord[0] = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);

  // phi
  polar_coord[1] = atan2(coord[1],coord[0]);

  return polar_coord;
 }

 /// Polar coordinates (r,phi) centered at the point (0.5,1)
 Vector<double> polar1(const Vector<double>& coord)
 {
  return polar(x1(coord));
 }

 /// Polar coordinates (r,phi) centered at the point (1.5,1)
 Vector<double> polar2(const Vector<double>& coord)
 {
  return polar(x2(coord));
 }
 
 /// function that contributs to u_exact
 double f1_exact(const Vector<double>& coord)
 {
  // Polar coordinates centered at the point (0.5,1)
  double r1 = polar1(coord)[0];
  double phi1 = polar1(coord)[1];
  
  return 1.0-sqrt(r1)*abs(sin(phi1/2.0));
 } // End of function
 
 /// f1 function, in front of the C1 unknown
 double f1(const Vector<double>& coord)
 {
  // Compute r1
  double r1 = polar1(coord)[0];
  return f1_exact(coord)*b(r1);
 }
 
 // gradient of f1 function, computed by sympy
 Vector<double> grad_f1(const Vector<double>& coord)
 {
  using namespace MathematicalConstants;
  Vector<double> df1(2);
  
  // Cartesian coordinates
  double x = coord[0];
  double y = coord[1];
  if (y>1.0)
   {
    y=1.0;
   }
  if (not(Blend))
   {
    /// Without blending
    df1[0] =  -(x/2.0 - 0.25)*sin(0.5*atan2(-y + 1.0, -x + 0.5))
     /pow((-x + 0.5)*(-x + 0.5) +(-y + 1.0)*(-y + 1.0),3.0/4.0) + 0.5*(y - 1.0)
     *cos(0.5*atan2(-y + 1.0, -x + 0.5))
     /pow((-x + 0.5)*(-x + 0.5)+ (-y + 1.0)*(-y + 1.0),3.0/4.0);
    
    df1[1] = 0.5*(-x + 0.5)*cos(0.5*atan2(-y + 1.0, -x + 0.5))
     /pow((-x + 0.5)*(-x + 0.5)
          + (-y + 1.0)*(-y + 1.0),3.0/4.0) - (y/2.0 - 1.0/2.0)
     *sin(0.5*atan2(-y + 1.0, -x + 0.5))
     /pow((-x + 0.5)*(-x + 0.5) + (-y + 1.0)*(-y + 1.0),3.0/4.0);
   }
  
  else
   {
    /// With blending
    // Compute r1
    double r1 = polar1(coord)[0];
    
    if (r1>R_blend)
     {
      df1[0] = 0.0;
      df1[1] = 0.0;
     }
    
    else
     {
      df1[0] = (-(x/2.0 - 0.25)*sin(0.5*atan2(-y + 1.0, -x + 0.5))
                /pow((-x + 0.5)*(-x + 0.5) + (-y + 1.0)*(-y + 1.0),3.0/4.0) + 0.5
                *(y - 1.0)*cos(0.5*atan2(-y + 1.0, -x + 0.5))
                /pow((-x + 0.5)*(-x + 0.5) + (-y + 1.0)*(-y + 1.0),3.0/4.0))
       *(0.5*cos(Pi*sqrt((-x + 0.5)*(-x + 0.5) + (-y + 1.0)*(-y + 1.0))/R_blend) + 0.5)
       - 0.5*Pi*(x - 0.5)*(-pow((-x + 0.5)*(-x + 0.5) + (-y + 1.0)*(-y + 1.0),1.0/4.0)
                      *sin(0.5*atan2(-y + 1.0, -x + 0.5)) + 1.0)
       *sin(Pi*sqrt((-x + 0.5)*(-x + 0.5) + (-y + 1.0)*(-y + 1.0))/R_blend)
       /(R_blend*sqrt((-x + 0.5)*(-x + 0.5) + (-y + 1.0)*(-y + 1.0)));
      
      df1[1] = (0.5*(-x + 0.5)*cos(0.5*atan2(-y + 1.0, -x + 0.5))
                /pow(pow(-x + 0.5,2) + pow(-y + 1.0,2),3.0/4.0)
                - (y/2.0 - 1.0/2.0)*sin(0.5*atan2(-y + 1.0, -x + 0.5))
                /pow(pow(-x + 0.5,2) + pow(-y + 1.0,2),3.0/4.0))
       *(0.5*cos(Pi*sqrt(pow(-x + 0.5,2) + pow(-y + 1.0,2))/R_blend) + 0.5)
       - 0.5*Pi*(y - 1.0)*(-pow(pow(-x + 0.5,2) + pow(-y + 1.0,2),1.0/4.0)
                      *sin(0.5*atan2(-y + 1.0, -x + 0.5)) + 1.0)
       *sin(Pi*sqrt(pow(-x + 0.5,2) + pow(-y + 1.0,2))/R_blend)
       /(R_blend*sqrt(pow(-x + 0.5,2) + pow(-y + 1.0,2)));
     }
     
   }
  return df1;
 }

 /// function that contributes to u_exact
 double f2_exact(const Vector<double>& coord)
 {

  // Polar coordinates centered at the point (1.5,1)
  double r2 = polar2(coord)[0];
  double phi2 = polar2(coord)[1];

  return 1.0-sqrt(r2)*cos(phi2/2.0);
 } // End of function

 /// f2 function, in front of the C2 unknown
 double f2(const Vector<double>& coord)
 {
  // Compute r2
  double r2 = polar2(coord)[0];

  return f2_exact(coord)*b(r2);
 }

 /// gradient of f2 function
 Vector<double> grad_f2(const Vector<double>& coord)
 {
  using namespace MathematicalConstants;
  Vector<double> df2(2);

  // Cartesian coordinates
  double x = coord[0];
  double y = coord[1];

  if (not(Blend))
   {
    /// Without blending
    df2[0] = -(x/2.0 - 0.75)*cos(0.5*atan2(-y + 1.0, -x + 1.5))
     /pow(pow(-x + 1.5,2) + pow(-y + 1.0,2),3.0/4.0) - 0.5*(y - 1.0)
     *sin(0.5*atan2(-y + 1.0, -x + 1.5))
     /pow(pow(-x + 1.5,2) + pow(-y + 1.0,2),3.0/4.0);
    
    df2[1] = -0.5*(-x + 1.5)*sin(0.5*atan2(-y + 1.0, -x + 1.5))
     /pow(pow(-x + 1.5,2) + pow(-y + 1.0,2),3.0/4.0) - (y/2.0 - 1.0/2.0)
     *cos(0.5*atan2(-y + 1.0, -x + 1.5))
     /pow(pow(-x + 1.5,2) + pow(-y + 1.0,2),3.0/4.0);
   }

  else
   {
    /// With blending
    
    // Compute r2
    double r2 = polar2(coord)[0];
    
    if (r2>R_blend)
     {
      df2[0] = 0.0;
      df2[1] = 0.0;
     }
    else
     {
      df2[0] = (-(x/2.0 - 0.75)*cos(0.5*atan2(-y + 1.0, -x + 1.5))
                /pow(pow(-x + 1.5,2) + pow(-y + 1.0,2),3.0/4.0) - 0.5*(y - 1.0)
                *sin(0.5*atan2(-y + 1.0, -x + 1.5))
                /pow(pow(-x + 1.5,2)+ pow(-y + 1.0,2),3.0/4.0))*(0.5*cos(Pi*sqrt(pow(-x + 1.5,2) + pow(-y + 1.0,2))/R_blend) + 0.5) - 0.5*Pi*(x - 1.5)*(-pow(pow(-x + 1.5,2) + pow(-y + 1.0,2),1.0/4.0)*cos(0.5*atan2(-y + 1.0, -x + 1.5)) + 1.0)*sin(Pi*sqrt(pow(-x + 1.5,2) + pow(-y + 1.0,2))/R_blend)/(R_blend*sqrt(pow(-x + 1.5,2) + pow(-y + 1.0,2)));
      
      df2[1] = (-0.5*(-x + 1.5)*sin(0.5*atan2(-y + 1.0, -x + 1.5))/pow(pow(-x + 1.5,2) + pow(-y + 1.0,2),3.0/4.0) - (y/2.0 - 1.0/2.0)*cos(0.5*atan2(-y + 1.0, -x + 1.5))/pow(pow(-x + 1.5,2) + pow(-y + 1.0,2),3.0/4.0))*(0.5*cos(Pi*sqrt(pow(-x + 1.5,2) + pow(-y + 1.0,2))/R_blend) + 0.5) - 0.5*Pi*(y - 1.0)*(-pow(pow(-x + 1.5,2) + pow(-y + 1.0,2),1.0/4.0)*cos(0.5*atan2(-y + 1.0, -x + 1.5)) + 1.0)*sin(Pi*sqrt(pow(-x + 1.5,2) + pow(-y + 1.0,2))/R_blend)/(R_blend*sqrt(pow(-x + 1.5,2) + pow(-y + 1.0,2)));
     }
   }
  return df2;
 }

 // Exact solution of the problem
 void get_exact_u(const Vector<double>& coord, Vector<double>& u)
 {
  using namespace MathematicalConstants;
  u[0] = sin((Pi/2.0)*coord[1])*cos(coord[0]) + f1_exact(coord) + f2_exact(coord);
 }

 /// Source function required to make the solution above an exact solution 
 void source_function(const Vector<double>& coord, double& source)
 {
  using namespace MathematicalConstants;

  /// u_FE = sin((Pi/2)y)cos(x)
  double x = coord[0];
  double y = coord[1];

  source = (x - 1.5)*(3.0*x - 4.5)*cos(0.5*atan2(-y + 1.0, -x + 1.5))
   /(4.0*pow((x - 1.5)*(x - 1.5) + (y - 1.0)*(y - 1.0),7.0/4.0)) - 0.25
   *(x - 1.5)*(y - 1.0)*sin(0.5*atan2(-y + 1.0, -x + 1.5))
   /pow((x - 1.5)*(x - 1.5) + (y - 1.0)*(y - 1.0),7.0/4.0) + (x - 0.5)
   *(3.0*x - 1.5)*sin(0.5*atan2(-y + 1.0, -x + 0.5))
   /(4.0*pow((x - 0.5)*(x - 0.5) + (y - 1.0)*(y - 1.0),7.0/4.0)) + (x - 0.5)
   *(y - 1.0)*cos(0.5*atan2(-y + 1.0, -x + 0.5))
   /(4.0*pow((x - 0.5)*(x - 0.5) + (y - 1.0)*(y - 1.0),7.0/4.0)) + (3.0*x - 4.5)
   *(y - 1.0)*sin(0.5*atan2(-y + 1.0, -x + 1.5))
   /(4.0*pow((x - 1.5)*(x - 1.5) + (y - 1.0)*(y - 1.0),7.0/4.0)) - 0.25
   *(3.0*x - 1.5)*(y - 1.0)*cos(0.5*atan2(-y + 1.0, -x + 0.5))
   /pow((x - 0.5)*(x - 0.5) + (y - 1.0)*(y - 1.0),7.0/4.0) + (y - 1.0)
   *(y - 1.0)*sin(0.5*atan2(-y + 1.0, -x + 0.5))
   /(4.0*pow((x - 0.5)*(x - 0.5) + (y - 1.0)*(y - 1.0),7.0/4.0)) + (y - 1.0)
   *(y - 1.0)*cos(0.5*atan2(-y + 1.0, -x + 1.5))
   /(4.0*pow((x - 1.5)*(x - 1.5) + (y - 1.0)*(y - 1.0),7.0/4.0))
   + ((x - 1.5)*(x - 1.5)*cos(0.5*atan2(-y + 1.0, -x + 1.5))
      /pow((x - 1.5)*(x - 1.5) + (y - 1.0)*(y - 1.0),7.0/4.0) - 2.0*(x - 1.5)
      *(y - 1.0)*sin(0.5*atan2(-y + 1.0, -x + 1.5))
      /pow((x - 1.5)*(x - 1.5) + (y - 1.0)*(y - 1.0),7.0/4.0) + (x - 0.5)
      *(x - 0.5)*sin(0.5*atan2(-y + 1.0, -x + 0.5))
      /pow((x - 0.5)*(x - 0.5)+ (y - 1.0)*(y - 1.0),7.0/4.0) + 2.0
      *(x - 0.5)*(y - 1.0)*cos(0.5*atan2(-y + 1.0, -x + 0.5))
      /pow((x - 0.5)*(x - 0.5) + (y - 1.0)*(y - 1.0),7.0/4.0) + 3.0
      *(y - 1.0)*(y - 1.0)*sin(0.5*atan2(-y + 1.0, -x + 0.5))
      /pow((x - 0.5)*(x - 0.5) + (y - 1.0)*(y - 1.0),7.0/4.0) + 3.0
      *(y - 1.0)*(y - 1.0)*cos(0.5*atan2(-y + 1.0, -x + 1.5))
      /pow((x - 1.5)*(x - 1.5) + (y - 1.0)*(y - 1.0),7.0/4.0) - 1.0
      *Pi*Pi*sin(0.5*Pi*y)*cos(x) - 2.0*sin(0.5*atan2(-y + 1.0, -x + 0.5))
      /pow((x - 0.5)*(x - 0.5) + (y - 1.0)*(y - 1.0),3.0/4.0) - 2.0
      *cos(0.5*atan2(-y + 1.0, -x + 1.5))
      /pow((x - 1.5)*(x - 1.5) + (y - 1.0)*(y - 1.0),3.0/4.0))
   /4.0 - sin(0.5*Pi*y)*cos(x) - sin(0.5*atan2(-y + 1.0, -x + 0.5))
   /(2.0*pow((x - 0.5)*(x - 0.5) + (y - 1.0)*(y - 1.0),3.0/4.0))
   - cos(0.5*atan2(-y + 1.0, -x + 1.5))
   /(2.0*pow((x - 1.5)*(x - 1.5) + (y - 1.0)*(y - 1.0),3.0/4.0));
 }
} 

//====== start_of_problem_class=======================================
/// 2D Poisson problem on rectangular domain, discretised with
/// 2D QPoisson elements. The specific type of element is
/// specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class PoissonProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source function
 PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt);

 /// Destructor (empty)
 ~PoissonProblem(){}

 ///  Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve(){} 

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 ///  Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);

private:

 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

 ///  Pointer to the element defining the unknown C1 (Note: eqn element is
 /// templated by the wrapped element!)
 SingularPoissonSolutionElement<ELEMENT>* Singular_poisson_solution1_element_pt;
 
 ///  Pointer to the element defining the unknown C2 (Note: eqn element is
 /// templated by the wrapped element!)
 SingularPoissonSolutionElement<ELEMENT>* Singular_poisson_solution2_element_pt;

 ///  Pointer to the C mesh associated with the elements defining the
 /// unknowns C1 and C2
 Mesh* Singular_poisson_solution_mesh_pt;

 /// Pointer to mesh with Poisson elements
 Mesh* Poisson_mesh_pt;

}; // end of problem class




//=====start_of_constructor===============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT>
PoissonProblem<ELEMENT>::
      PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt)
       :  Source_fct_pt(source_fct_pt)
{
 
 // Setup mesh
 
 // # of elements in x-direction
 unsigned n_x= 8*Global_parameters::Element_multiplier; 

 // # of elements in y-direction
 unsigned n_y= 4*Global_parameters::Element_multiplier;

 // Domain length in x-direction
 double l_x=2.0; 
 
 // Domain length in y-direction
 double l_y=1.0;

 // Create a SimpleRectangularQuadMesh object
 Poisson_mesh_pt = new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Element width
 double dx=l_x/double(n_x);

 // Find the x-coordinate of the points subject to a singularity
 double unpinned_fraction=0.5;
 double x_left=unsigned(0.5*l_x*(1.0-unpinned_fraction)/dx)*dx;
 double x_right=unsigned((l_x-0.5*l_x*(1.0-unpinned_fraction))/dx)*dx;
 oomph_info << "x_left, x_right " << x_left << " " << x_right << std::endl;


#ifdef DONT_USE_SINGULARITY

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet conditions
 // here. 
 unsigned n_bound = Poisson_mesh_pt->nboundary();
 for(unsigned i=0;i<n_bound;i++)
  {
   unsigned n_node = Poisson_mesh_pt->nboundary_node(i);
   for (unsigned n=0;n<n_node;n++)
    {
     Node* nod_pt=Poisson_mesh_pt->boundary_node_pt(i,n);

     // Pin
     nod_pt->pin(0); 

     // Get exact solution
     Vector<double> x(2);
     x[0] = nod_pt->x(0);
     x[1] = nod_pt->x(1);
     Vector<double> value_imposed(1);
     Global_parameters::get_exact_u(x,value_imposed);

     // Impose that value for the node
     nod_pt->set_value(0,value_imposed[0]);
    }
  }

 // Very hacky... Unpin central nodes on top boundary (boundary 2)
 //---------------------------------------------------------------

 // Unpin nodes that are genuinely inside elements
 unsigned i=2;
 unsigned n_node = Poisson_mesh_pt->nboundary_node(i);
 for (unsigned n=0;n<n_node;n++)
  {
   double x=Poisson_mesh_pt->boundary_node_pt(i,n)->x(0); 
   if ((x>x_left)&&(x<x_right))
    {
     oomph_info << "Unpinning at x = " << x << std::endl;
     Poisson_mesh_pt->boundary_node_pt(i,n)->unpin(0); 
    }
  }

#else

 //--------------------------------------------------------------------
 /// Boundary conditions
 //--------------------------------------------------------------------
 
 // Number of nodes along element edge
 unsigned nnode_1d=dynamic_cast<ELEMENT*>(Poisson_mesh_pt->element_pt(0))
  ->nnode_1d();
 
 // Boundary 0
 
 // Loop over the elements
 for (unsigned e=0;e<n_x;e++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Poisson_mesh_pt->element_pt(e));
   // Loop over the boundary nodes
   for (unsigned j=0;j<nnode_1d;j++)
    {
     // The node j is subject to Dirichlet BC
     el_pt->impose_dirichlet_bc_on_node(j);
     // Find the global coordinate at node j
     Vector<double> x(2);
     x[0] = el_pt->nodal_position(j,0);
     x[1] = el_pt->nodal_position(j,1);
     // Find the Dirichlet value at this local coordinate
     Vector<double> value_imposed(1);
     Global_parameters::get_exact_u(x,value_imposed);
     // Impose that value for the node
     el_pt->set_dirichlet_value_on_node(j,value_imposed[0]);
    }
  }

 // Boundary one

 // Loop over the elements
 for (unsigned i=0;i<n_y;i++)
  {
   unsigned e = i*n_x + n_x - 1;
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Poisson_mesh_pt->element_pt(e));
   // Loop over the nodes
   for (unsigned k=0;k<nnode_1d;k++)
    {
     unsigned j =
      k*nnode_1d + nnode_1d - 1;
     // The node j is subject to Dirichlet BC
     el_pt->impose_dirichlet_bc_on_node(j);
     // Find the global coordinate at node j
     Vector<double> x(2);
     x[0] = el_pt->nodal_position(j,0);
     x[1] = el_pt->nodal_position(j,1);
     // Find the Dirichlet value at this local coordinate
     Vector<double> value_imposed(1);
     Global_parameters::get_exact_u(x,value_imposed);
     // Impose that value for the node
     el_pt->set_dirichlet_value_on_node(j,value_imposed[0]);
    }
  }
 
 // Boundary three

 // Loop over the elements
 for (unsigned i=0;i<n_y;i++)
  {
   unsigned e = i*n_x;           
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Poisson_mesh_pt->element_pt(e));
   // Loop over the nodes
   for (unsigned k=0;k<nnode_1d;k++)
    {
     unsigned j = k*nnode_1d;
     // The node j is subject to Dirichlet BC
     el_pt->impose_dirichlet_bc_on_node(j);
     // Find the global coordinate at node j
     Vector<double> x(2);
     x[0] = el_pt->nodal_position(j,0);
     x[1] = el_pt->nodal_position(j,1);
     // Find the Dirichlet value at this local coordinate
     Vector<double> value_imposed(1);
     Global_parameters::get_exact_u(x,value_imposed);
     // Impose that value for the node
     el_pt->set_dirichlet_value_on_node(j,value_imposed[0]);
    }
  }

 
 // Boundary two
 // For the moment we set the Dirichlet BC for all the nodes of this boundary
 
 // Loop over the elements
 for (unsigned i=0;i<n_x;i++)
  {
   unsigned e = (n_y-1)*n_x + i;
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Poisson_mesh_pt->element_pt(e));
   // Loop over the boundary nodes
   for (unsigned k=0;k<nnode_1d;k++)
    {
     unsigned j =
      (nnode_1d-1)*nnode_1d + k;
     // The node j is subject to Dirichlet BC
     el_pt->impose_dirichlet_bc_on_node(j);
     // Find the global coordinate at node j
     Vector<double> x(2);
     x[0] = el_pt->nodal_position(j,0);
     x[1] = el_pt->nodal_position(j,1);
     // Find the Dirichlet value at this local coordinate
     Vector<double> value_imposed(1);
     Global_parameters::get_exact_u(x,value_imposed);
     // Impose that value for the node
     el_pt->set_dirichlet_value_on_node(j,value_imposed[0]);
    }
  }

 // Find the numbers of the elements with the boundary nodes subject to the
 // Neumann BC

 // We begin by the first element on the left of the upper line of element
 unsigned n_el = n_x*(n_y-1);
 // Upcast from GeneralisedElement to the present element 
 ELEMENT *first_el_pt = dynamic_cast<ELEMENT*>(Poisson_mesh_pt->element_pt(n_el));
 // Initialise the x coordinate of the top right node of the first neumann element
 double xstart = first_el_pt->nodal_position((nnode_1d-1)
                                             *nnode_1d,0);
 // Initialise the number of the first neumann element on the upper row
 unsigned istart = 0;
 // Loop that increases istart until n_el+istart is the number of the first
 // neumann element
 while (xstart < x_left)
  {
   // Increases the value of i
   istart++;
   // Upcast from GeneralisedElement to the present element 
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>
    (Poisson_mesh_pt->element_pt(n_el+istart));
   // Update the value of the x coordinate of the top right node of the right
   // element
   xstart = el_pt->nodal_position((nnode_1d-1)
                                  *nnode_1d,0);
  }
 // The number of the first neumann element is n_el+istart
 // We now want to find the number of the last neumann element

 // We begin by the istart-th element on the left of the upper line of elements
 // Upcast from GeneralisedElement to the present element 
 ELEMENT* istartth_el_pt = dynamic_cast<ELEMENT*>
  (Poisson_mesh_pt->element_pt(n_el+istart));
 // Initialise the x coordinate of the top right node of the last neumann element
 double xlast = istartth_el_pt->nodal_position((nnode_1d-1)
                                               *nnode_1d,0);
 // Initialise the number of the last neumann element on the upper row
 unsigned ilast = istart;
 // Loop that increases ilast until n_el+istart is the number of the last element
 while (xlast < x_right)
  {
   // Increases the value of i
   ilast++;
   // Upcast from GeneralisedElement to the present element 
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>
    (Poisson_mesh_pt->element_pt(n_el+ilast));
   // Update the value of the x coordinate of the top right node of the right
   // element
   xlast = el_pt->nodal_position((nnode_1d-1)
                                 *nnode_1d,0);
  }
 ilast--;
 // The number of the last neumann element is n_el+ilast
 // We will now remove the Dirichlet BC on the nodes of the neumann elements
 
 // Boundary two between x=0.5 and x=1.5

 // Loop over the elements
 for (unsigned i=istart;i<=ilast;i++)
  {
   unsigned e = n_x*(n_y-1) + i;
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Poisson_mesh_pt->element_pt(e));
   if (i==istart)
    {
     // Loop over the boundary nodes except the first one
     for (unsigned k=1;k<=nnode_1d-1;k++)
      {
       unsigned j =
        (nnode_1d-1)*nnode_1d + k;
       // Undo BC
       el_pt->undo_dirichlet_bc_on_node(j);
      }
    }
   else if (i==ilast)
    {
     // Loop over the boundary nodes except the last one
     for (unsigned k=0;k<=nnode_1d-2;k++)
      {
       unsigned j =
        (nnode_1d-1)*nnode_1d + k;
       // Undo BC
       el_pt->undo_dirichlet_bc_on_node(j);
      }
    }
   else
    {
     el_pt->undo_dirichlet_bc_on_all_nodes();
    }
  }
 
#endif

  //--------------------------------------------------------------------
 /// End of boundary conditions
 //--------------------------------------------------------------------

 // Add Poisson mesh to global mesh
 add_sub_mesh(Poisson_mesh_pt);

#ifndef DONT_USE_SINGULARITY

 // Create a C mesh 
 Singular_poisson_solution_mesh_pt=new Mesh;

 // Create a pointer to a C1EquationElement object
 Singular_poisson_solution1_element_pt = new SingularPoissonSolutionElement<ELEMENT>;
 
 // Set the pointer to the f1 singular function for this element, defined in
 // Global_parameters namespace
 Singular_poisson_solution1_element_pt->singular_fct_pt() = &Global_parameters::f1;

 // Set the pointer to the gradient of the f1 singular function for
 // this element, defined in Global_parameters namespace
 Singular_poisson_solution1_element_pt->grad_singular_fct_pt() =
  &Global_parameters::grad_f1;
 
 // Singular function satisfies Laplace's eqn
 if (!Global_parameters::Blend)
  {
   Singular_poisson_solution1_element_pt->singular_function_satisfies_laplace_equation()=true;
  }

 // Tell the C1EquationElement object about its associated Poisson element
 // and the local coordinate in the poisson element at which the residual
 // will be computed
 Vector<double> s1_pt(2);
 s1_pt[0] = -1.0;
 s1_pt[1] = 1.0;
 ELEMENT *el1_pt = dynamic_cast<ELEMENT*>
  (Poisson_mesh_pt->element_pt(n_el+istart));
 Singular_poisson_solution1_element_pt->set_wrapped_poisson_element_pt
  (el1_pt,s1_pt,&Global_parameters::Direction);

 // Add element to the C mesh
 Singular_poisson_solution_mesh_pt->add_element_pt(Singular_poisson_solution1_element_pt);

 // Create a pointer to a C2EquationElement object
 Singular_poisson_solution2_element_pt = new SingularPoissonSolutionElement<ELEMENT>;
 
 // Set the pointer to the f2 singular function for this element, defined in
 // Global_parameters namespace
 Singular_poisson_solution2_element_pt->singular_fct_pt() = &Global_parameters::f2;

 // Set the pointer to the gradient of the f2 singular function for
 // this element, defined in Global_parameters namespace
 Singular_poisson_solution2_element_pt->grad_singular_fct_pt() =
  &Global_parameters::grad_f2;
 
 // Singular function satisfies Laplace's eqn
 if (!Global_parameters::Blend)
  {
   Singular_poisson_solution2_element_pt->singular_function_satisfies_laplace_equation()=true;
  }

 // Tell the C2EquationElement object about its associated Poisson element
 // and the local coordinate in the poisson element at which the residual
 // will be computed
 Vector<double> s2_pt(2);
 s2_pt[0] = -1.0;
 s2_pt[1] = 1.0;
 ELEMENT *el2_pt = dynamic_cast<ELEMENT*>
  (Poisson_mesh_pt->element_pt(n_el+ilast+1));
 Singular_poisson_solution2_element_pt->set_wrapped_poisson_element_pt
  (el2_pt,s2_pt,&Global_parameters::Direction);
 
 // Add element to its own mesh
 Singular_poisson_solution_mesh_pt->add_element_pt(Singular_poisson_solution2_element_pt);

 // Add C mesh to global mesh
 add_sub_mesh(Singular_poisson_solution_mesh_pt);

#endif

 // Build global mesh
 build_global_mesh();

 // Complete the build of all elements so they are fully functional

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor: Pass pointer to source function
 unsigned n_element = Poisson_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Poisson_mesh_pt->element_pt(i));

   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;

#ifndef DONT_USE_SINGULARITY

   // Set the pointer to the element that determines the amplitude
   // of the singular fct
   el_pt->add_c_equation_element_pt(Singular_poisson_solution1_element_pt);

   // Set the pointer to the element that determines the amplitude
   // of the singular fct
   el_pt->add_c_equation_element_pt(Singular_poisson_solution2_element_pt);   

#endif

  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;

} // end of constructor



//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts;

 // Output solution with few plot points
 //-------------------------------------
 npts=5;
 sprintf(filename,"%s/soln_coarse%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Poisson_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output solution with lots of plot points
 //-----------------------------------------
 npts=5;
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Poisson_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Poisson_mesh_pt->output_fct(some_file,npts,Global_parameters::get_exact_u); 
 some_file.close();

 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Poisson_mesh_pt->compute_error(some_file,Global_parameters::get_exact_u,
                          error,norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

#ifndef DONT_USE_SINGULARITY

 // Print the value of C1
 oomph_info << "Value of C1 = "
               << Singular_poisson_solution1_element_pt->c() << std::endl;

 // Print the value of C2
 oomph_info << "Value of C2 = "
               << Singular_poisson_solution2_element_pt->c() << std::endl;
    
#endif
 
} // end of doc

 




//===== start_of_main=====================================================
/// Driver code for 2D Poisson problem
//========================================================================
int main(int argc, char **argv)
{
 
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Multiplier for elements
 CommandLineArgs::specify_command_line_flag(
  "--element_multiplier",
  &Global_parameters::Element_multiplier);

 // Impose blending or not
 CommandLineArgs::specify_command_line_flag("--blend");

 // radius of the blending
 CommandLineArgs::specify_command_line_flag(
  "--r_blend",&Global_parameters::R_blend);
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();


 if (CommandLineArgs::command_line_flag_has_been_set
     ("--blend"))
  {
   Global_parameters::Blend=true;
  }
 
 //Set up the problem
 //------------------

 // Create the problem with 2D nine-node elements from the
 // QPoissonElement family. Pass pointer to source function. 
 PoissonProblem<PoissonElementWithSingularity<
  QPoissonElement<2,3> > > 
  problem(&Global_parameters::source_function);

 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Solve the problem
 problem.newton_solve();
 
 //Output the solution
 problem.doc_solution(doc_info);

} //end of main









