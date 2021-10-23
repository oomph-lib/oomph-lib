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
// Driver code for a simple test poisson problem using a mesh
// generated inline by the triangle mesh generator Triangle.

//Generic routines
#include "generic.h"

// Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;

using namespace oomph;

//====================================================================
/// Namespace for exact solution for Poisson equation with sharp step 
//====================================================================
namespace TanhSolnForPoisson
{

 /// Parameter for steepness of step
 double Alpha;

 /// Parameter for angle of step
 double Beta;

 
 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=tanh(1.0-Alpha*(Beta*x[0]-x[1]));
 }


 /// Exact solution as a scalar
 void get_exact_u(const Vector<double>& x, double& u)
 {
  u=tanh(1.0-Alpha*(Beta*x[0]-x[1]));
 }


 /// Source function to make it an exact solution 
 void get_source(const Vector<double>& x, double& source)
 {
  source = 2.0*tanh(-1.0+Alpha*(Beta*x[0]-x[1]))*
             (1.0-pow(tanh(-1.0+Alpha*(Beta*x[0]-x[1])),2.0))*
             Alpha*Alpha*Beta*Beta+2.0*tanh(-1.0+Alpha*(Beta*x[0]-x[1]))*
            (1.0-pow(tanh(-1.0+Alpha*(Beta*x[0]-x[1])),2.0))*Alpha*Alpha;
 }

}


//====================================================================
/// Micky mouse Poisson problem.
//====================================================================

// Poisson problem
template<class ELEMENT> 
class PoissonProblem : public Problem
{

public:


 /// Constructor: Pass pointer to source function and names of
 /// two triangle input files
 PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt);

 /// Destructor (empty)
 ~PoissonProblem(){

 }

 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve()
  {
   apply_boundary_conditions();
  }

 /// Update the problem specs before solve (empty)
 void actions_after_newton_solve(){}

 void actions_before_adapt() {}
 void actions_after_adapt()
  {
   complete_problem_setup();
  }

 /// Access function for the specific mesh
 RefineableTriangleMesh<ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RefineableTriangleMesh<ELEMENT>*>(Problem::mesh_pt());
  }

  /// Doc the solution
 void doc_solution(DocInfo& doc_info, ofstream &trace_file);

 /// Error estimator
 Z2ErrorEstimator* error_estimator_pt;


private:

 RefineableTriangleMesh<ELEMENT> *My_mesh_pt;

 void apply_boundary_conditions();

 void complete_problem_setup();

 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

};



//========================================================================
/// Constructor for Poisson problem
//========================================================================
template<class ELEMENT>
PoissonProblem<ELEMENT>::
 PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt)
        : Source_fct_pt(source_fct_pt)
{ 

 // >> ====================================================================
 // >> Setup boundary (Begin)
 // >> ====================================================================

 unsigned boundary_id;
 
 // ************************************************************************
 // Outer boundary
 // ************************************************************************

 Vector<TriangleMeshClosedCurve*> outer_boundaries_pt(2);
 
 // ************************************************************************
 // First outer boundary
 // ************************************************************************

 Vector<TriangleMeshCurveSection*> squares_pt(8);

 // ------------------------------------------------------------------------
 Vector<Vector<double> >verticesb1(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb1[i].resize(2);
  }

 verticesb1[0][0] = 1.0;
 verticesb1[0][1] = 1.0;

 verticesb1[1][0] = 1.0;
 verticesb1[1][1] = 1.5;
 
 verticesb1[2][0] = 1.0;
 verticesb1[2][1] = 2.0;
 
 boundary_id = 1;

 TriangleMeshPolyLine *boundaryp1 =
  new TriangleMeshPolyLine(verticesb1, boundary_id);

 squares_pt[0] = boundaryp1;

 // ------------------------------------------------------------------------

 Vector<Vector<double> >verticesb2(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb2[i].resize(2);
  }

 verticesb2[0][0] = 1.0;
 verticesb2[0][1] = 2.0;

 verticesb2[1][0] = 1.5;
 verticesb2[1][1] = 2.0;

 verticesb2[2][0] = 2.0;
 verticesb2[2][1] = 2.0;

 boundary_id = 2;

 TriangleMeshPolyLine *boundaryp2 =
  new TriangleMeshPolyLine(verticesb2, boundary_id);

 squares_pt[1] = boundaryp2;

 // ------------------------------------------------------------------------

 Vector<Vector<double> >verticesb3(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb3[i].resize(2);
  }

 verticesb3[0][0] = 2.0;
 verticesb3[0][1] = 2.0;

 verticesb3[1][0] = 2.0;
 verticesb3[1][1] = 1.5;

 verticesb3[2][0] = 2.0;
 verticesb3[2][1] = 1.0;

 boundary_id = 3;

 TriangleMeshPolyLine *boundaryp3 =
  new TriangleMeshPolyLine(verticesb3, boundary_id);

 squares_pt[2] = boundaryp3;

 // ------------------------------------------------------------------------

 Vector<Vector<double> >verticesb4(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb4[i].resize(2);
  }

 verticesb4[0][0] = 2.0;
 verticesb4[0][1] = 1.0;

 verticesb4[1][0] = 1.5;
 verticesb4[1][1] = 1.0;

 verticesb4[2][0] = 1.0;
 verticesb4[2][1] = 1.0;

 boundary_id = 4;

 TriangleMeshPolyLine *boundaryp4 =
  new TriangleMeshPolyLine(verticesb4, boundary_id);

 unsigned vertex_to_connect = 0;
 boundaryp4->connect_final_vertex_to_polyline(boundaryp1, vertex_to_connect);

 squares_pt[3] = boundaryp4;

 // ------------------------------------------------------------------------

 Vector<Vector<double> >verticesb5(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb5[i].resize(2);
  }

 verticesb5[0][0] = 1.0;
 verticesb5[0][1] = 1.0;

 verticesb5[1][0] = 0.5;
 verticesb5[1][1] = 1.0;

 verticesb5[2][0] = 0.0;
 verticesb5[2][1] = 1.0;

 boundary_id = 5;

 TriangleMeshPolyLine *boundaryp5 =
  new TriangleMeshPolyLine(verticesb5, boundary_id);

 vertex_to_connect = 0;
 boundaryp5->connect_initial_vertex_to_polyline(boundaryp1, vertex_to_connect);

 squares_pt[4] = boundaryp5;

 // ------------------------------------------------------------------------

 Vector<Vector<double> >verticesb6(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb6[i].resize(2);
  }

 verticesb6[0][0] = 0.0;
 verticesb6[0][1] = 1.0;

 verticesb6[1][0] = 0.0;
 verticesb6[1][1] = 0.5;

 verticesb6[2][0] = 0.0;
 verticesb6[2][1] = 0.0;

 boundary_id = 6;

 TriangleMeshPolyLine *boundaryp6 =
  new TriangleMeshPolyLine(verticesb6, boundary_id);

 squares_pt[5] = boundaryp6;

 // ------------------------------------------------------------------------

 Vector<Vector<double> >verticesb7(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb7[i].resize(2);
  }

 verticesb7[0][0] = 0.0;
 verticesb7[0][1] = 0.0;

 verticesb7[1][0] = 0.5;
 verticesb7[1][1] = 0.0;

 verticesb7[2][0] = 1.0;
 verticesb7[2][1] = 0.0;

 boundary_id = 7;

 TriangleMeshPolyLine *boundaryp7 =
  new TriangleMeshPolyLine(verticesb7, boundary_id);

 squares_pt[6] = boundaryp7;

 // ------------------------------------------------------------------------

 Vector<Vector<double> >verticesb8(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb8[i].resize(2);
  }

 verticesb8[0][0] = 1.0;
 verticesb8[0][1] = 0.0;

 verticesb8[1][0] = 1.0;
 verticesb8[1][1] = 0.5;

 verticesb8[2][0] = 1.0;
 verticesb8[2][1] = 1.0;

 boundary_id = 8;

 TriangleMeshPolyLine *boundaryp8 =
  new TriangleMeshPolyLine(verticesb8, boundary_id);

 squares_pt[7] = boundaryp8;

 // ------------------------------------------------------------------------
 
 outer_boundaries_pt[0] = new TriangleMeshPolygon(squares_pt);
 
 // ************************************************************************
 // Second outer boundary
 // ************************************************************************

 Vector<TriangleMeshCurveSection*> triangles_pt(6);

 // ------------------------------------------------------------------------
 Vector<Vector<double> >verticesb9(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb9[i].resize(2);
  }
 
 verticesb9[0][0] = 3.0;
 verticesb9[0][1] = 1.0;

 verticesb9[1][0] = 4.0;
 verticesb9[1][1] = 1.0;

 verticesb9[2][0] = 5.0;
 verticesb9[2][1] = 1.0;

 boundary_id = 9;

 TriangleMeshPolyLine *boundaryp9 =
  new TriangleMeshPolyLine(verticesb9, boundary_id);

 triangles_pt[0] = boundaryp9;

 // ------------------------------------------------------------------------
 Vector<Vector<double> >verticesb10(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb10[i].resize(2);
  }
 
 verticesb10[0][0] = 5.0;
 verticesb10[0][1] = 1.0;

 verticesb10[1][0] = 4.5;
 verticesb10[1][1] = 1.5;

 verticesb10[2][0] = 4.0;
 verticesb10[2][1] = 2.0;

 boundary_id = 10;

 TriangleMeshPolyLine *boundaryp10 =
  new TriangleMeshPolyLine(verticesb10, boundary_id);

 triangles_pt[1] = boundaryp10;
 
// ------------------------------------------------------------------------
 Vector<Vector<double> >verticesb11(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb11[i].resize(2);
  }
 
 verticesb11[0][0] = 4.0;
 verticesb11[0][1] = 2.0;

 verticesb11[1][0] = 3.5;
 verticesb11[1][1] = 2.5;

 verticesb11[2][0] = 3.0;
 verticesb11[2][1] = 3.0;

 boundary_id = 11;

 TriangleMeshPolyLine *boundaryp11 =
  new TriangleMeshPolyLine(verticesb11, boundary_id);

 triangles_pt[2] = boundaryp11;
 
 // ------------------------------------------------------------------------
 Vector<Vector<double> >verticesb12(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb12[i].resize(2);
  }
 
 verticesb12[0][0] = 3.0;
 verticesb12[0][1] = 3.0;

 verticesb12[1][0] = 4.0;
 verticesb12[1][1] = 3.0;

 verticesb12[2][0] = 5.0;
 verticesb12[2][1] = 3.0;

 boundary_id = 12;

 TriangleMeshPolyLine *boundaryp12 =
  new TriangleMeshPolyLine(verticesb12, boundary_id);

 triangles_pt[3] = boundaryp12;

 // ------------------------------------------------------------------------
 Vector<Vector<double> >verticesb13(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb13[i].resize(2);
  }
 
 verticesb13[0][0] = 5.0;
 verticesb13[0][1] = 3.0;

 verticesb13[1][0] = 4.5;
 verticesb13[1][1] = 2.5;

 verticesb13[2][0] = 4.0;
 verticesb13[2][1] = 2.0;

 boundary_id = 13;

 TriangleMeshPolyLine *boundaryp13 =
  new TriangleMeshPolyLine(verticesb13, boundary_id);

 vertex_to_connect = 2;
 boundaryp13->connect_final_vertex_to_polyline(boundaryp10, vertex_to_connect);

 triangles_pt[4] = boundaryp13;

 // ------------------------------------------------------------------------
 Vector<Vector<double> >verticesb14(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb14[i].resize(2);
  }
 
 verticesb14[0][0] = 4.0;
 verticesb14[0][1] = 2.0;

 verticesb14[1][0] = 3.5;
 verticesb14[1][1] = 1.5;

 verticesb14[2][0] = 3.0;
 verticesb14[2][1] = 1.0;

 boundary_id = 14;

 TriangleMeshPolyLine *boundaryp14 =
  new TriangleMeshPolyLine(verticesb14, boundary_id);

 vertex_to_connect = 2;
 boundaryp14->connect_initial_vertex_to_polyline(boundaryp10, vertex_to_connect);

 triangles_pt[5] = boundaryp14;

 // ------------------------------------------------------------------------
 
 outer_boundaries_pt[1] = new TriangleMeshPolygon(triangles_pt);

 // ------------------------------------------------------------------------

 // ************************************************************************
 // Internal closed boundary
 // ************************************************************************

 Vector<TriangleMeshClosedCurve*> internal_boundaries_pt(2);

 // ************************************************************************
 // First internal boundary
 // ************************************************************************

 Vector<TriangleMeshCurveSection*> internal_circle1_pt(2);
 // ------------------------------------------------------------------------
 
 double x_centre = 1.5;
 double y_centre = 1.5;
 double r_circle = 0.2;
 Circle *circle1_pt = new Circle(x_centre, y_centre, r_circle);

 unsigned n_seg1 = 20;
 double z_start = 0.0;
 double z_end = MathematicalConstants::Pi;

 boundary_id = 15;
 TriangleMeshCurviLine *boundaryc15 = 
  new TriangleMeshCurviLine(circle1_pt,
                            z_start,
                            z_end,
                            n_seg1,
                            boundary_id);

 internal_circle1_pt[0] = boundaryc15;

 // ------------------------------------------------------------------------
 unsigned n_seg2 = 20;
 z_start = MathematicalConstants::Pi;
 z_end = 2.0*MathematicalConstants::Pi;

 boundary_id = 16;
 TriangleMeshCurviLine *boundaryc16 = 
  new TriangleMeshCurviLine(circle1_pt,
                            z_start,
                            z_end,
                            n_seg2,
                            boundary_id);

 internal_circle1_pt[1] = boundaryc16;

 // ------------------------------------------------------------------------
 Vector<double> hole1(2);
 hole1[0] = x_centre;
 hole1[1] = y_centre;

 internal_boundaries_pt[0] = 
  new TriangleMeshClosedCurve(internal_circle1_pt, hole1);

 // ************************************************************************
 // Second internal boundary
 // ************************************************************************

 Vector<TriangleMeshCurveSection*> internal_circle2_pt(2);
 // ------------------------------------------------------------------------
 
 x_centre = 4.0;
 y_centre = 1.4;
 r_circle = 0.2;
 Circle *circle2_pt = new Circle(x_centre, y_centre, r_circle);

 unsigned n_seg3 = 20;
 z_start = 0.0;
 z_end = MathematicalConstants::Pi;

 boundary_id = 17;
  TriangleMeshCurviLine *boundaryc17 = 
  new TriangleMeshCurviLine(circle2_pt,
                            z_start,
                            z_end,
                            n_seg3,
                            boundary_id);

 internal_circle2_pt[0] = boundaryc17;

 // ------------------------------------------------------------------------
 unsigned n_seg4 = 20;
 z_start = MathematicalConstants::Pi;
 z_end = 2.0*MathematicalConstants::Pi;

 boundary_id = 18;
 TriangleMeshCurviLine *boundaryc18 = 
  new TriangleMeshCurviLine(circle2_pt,
                            z_start,
                            z_end,
                            n_seg4,
                            boundary_id);

 internal_circle2_pt[1] = boundaryc18;

 // ------------------------------------------------------------------------
 Vector<double> hole2(2);
 hole2[0] = x_centre;
 hole2[1] = y_centre;
  
 internal_boundaries_pt[1] = 
  new TriangleMeshClosedCurve(internal_circle2_pt, hole2);

 // ************************************************************************
 // Internal open boundary
 // ************************************************************************

 // ------------------------------------------------------------------------
 Vector<TriangleMeshOpenCurve *> internal_open_boundaries_pt(3);
 // ------------------------------------------------------------------------
 
 // ************************************************************************
 // Internal open boundary 1
 // ************************************************************************

 // ------------------------------------------------------------------------
 Vector<TriangleMeshCurveSection*> internal_open_boundary1_pt(1);
 // ------------------------------------------------------------------------

 Vector<Vector<double> >verticesb19(2);
 for (unsigned i = 0; i < 2; i++)
  {
   verticesb19[i].resize(2);
  }

 verticesb19[0][0] = 0.5;
 verticesb19[0][1] = 1.0;

 verticesb19[1][0] = 0.5;
 verticesb19[1][1] = 0.0;

 boundary_id = 19;
 TriangleMeshPolyLine *boundaryp19 = 
  new TriangleMeshPolyLine(verticesb19, boundary_id);

 vertex_to_connect = 1;
 boundaryp19->
  connect_initial_vertex_to_polyline(boundaryp5, vertex_to_connect);

 boundaryp19->
  connect_final_vertex_to_polyline(boundaryp7, vertex_to_connect);
 
 internal_open_boundary1_pt[0] = boundaryp19;

 internal_open_boundaries_pt[0] = 
  new TriangleMeshOpenCurve(internal_open_boundary1_pt);
 
 // ************************************************************************
 // Internal open boundary 2
 // ************************************************************************

 // ------------------------------------------------------------------------
 Vector<TriangleMeshCurveSection*> internal_open_boundary2_pt(1);
 // ------------------------------------------------------------------------

 Vector<Vector<double> >verticesb20(3);
 for (unsigned i = 0; i < 3; i++)
  {
   verticesb20[i].resize(2);
  }

 verticesb20[0][0] = 3.5;
 verticesb20[0][1] = 2.5;

 verticesb20[1][0] = 4.0;
 verticesb20[1][1] = 2.5;

 verticesb20[2][0] = 4.5;
 verticesb20[2][1] = 2.5;

 boundary_id = 20;
 TriangleMeshPolyLine *boundaryp20 = 
  new TriangleMeshPolyLine(verticesb20, boundary_id);

 vertex_to_connect = 1;
 boundaryp20->
  connect_initial_vertex_to_polyline(boundaryp11, vertex_to_connect);

 boundaryp20->
  connect_final_vertex_to_polyline(boundaryp13, vertex_to_connect);

 internal_open_boundary2_pt[0] = boundaryp20;

 internal_open_boundaries_pt[1] = 
  new TriangleMeshOpenCurve(internal_open_boundary2_pt);

 // ************************************************************************
 // Internal open boundary 3
 // ************************************************************************

 // ------------------------------------------------------------------------
 Vector<TriangleMeshCurveSection*> internal_open_boundary3_pt(1);
 // ------------------------------------------------------------------------

 Vector<Vector<double> >verticesb21(2);
 for (unsigned i = 0; i < 2; i++)
  {
   verticesb21[i].resize(2);
  }

 verticesb21[0][0] = 4.0;
 verticesb21[0][1] = 3.0;

 verticesb21[1][0] = 4.0;
 verticesb21[1][1] = 2.5;

 boundary_id = 21;
 TriangleMeshPolyLine *boundaryp21 = 
  new TriangleMeshPolyLine(verticesb21, boundary_id);

 vertex_to_connect = 1;
 boundaryp21->
  connect_initial_vertex_to_polyline(boundaryp12, vertex_to_connect);

 boundaryp21->
  connect_final_vertex_to_polyline(boundaryp20, vertex_to_connect);
 
 internal_open_boundary3_pt[0] = boundaryp21;
 
 internal_open_boundaries_pt[2] = 
  new TriangleMeshOpenCurve(internal_open_boundary3_pt);

 // >> ====================================================================
 // >> Setup boundary (End)
 // >> ====================================================================

 // Setup parameters for exact tanh solution

 // Steepness of step
 TanhSolnForPoisson::Alpha=1.0;

 // Orientation of step
 TanhSolnForPoisson::Beta=1.0;

 // >> ====================================================================
 // >> Create mesh
 // >> ====================================================================

 double uniform_element_area=0.2;

 TriangleMeshParameters triangle_mesh_parameters(outer_boundaries_pt);
 // -----------------------------------------------------------------------
 triangle_mesh_parameters.internal_closed_curve_pt() = 
  internal_boundaries_pt;
 // -----------------------------------------------------------------------
 triangle_mesh_parameters.internal_open_curves_pt() = 
  internal_open_boundaries_pt;
 // -----------------------------------------------------------------------
 triangle_mesh_parameters.element_area() = uniform_element_area;
 // -----------------------------------------------------------------------
 // >> Regions
 Vector<double> region1(2);
 region1[0] = 0.25; region1[1] = 0.5;
 Vector<double> region2(2);
 region2[0] = 0.75; region2[1] = 0.5;
 Vector<double> region3(2);
 region3[0] = 1.5; region3[1] = 1.1;
 Vector<double> region4(2);
 region4[0] = 4.0; region4[1] = 1.1;
 Vector<double> region5(2);
 region5[0] = 4.0; region5[1] = 2.2;
 Vector<double> region6(2);
 region6[0] = 3.5; region6[1] = 2.8;
 Vector<double> region7(2);
 region7[0] = 4.5; region7[1] = 2.8;

 triangle_mesh_parameters.add_region_coordinates(10, region1);
 triangle_mesh_parameters.add_region_coordinates(20, region2);
 triangle_mesh_parameters.add_region_coordinates(30, region3);
 triangle_mesh_parameters.add_region_coordinates(40, region4);
 triangle_mesh_parameters.add_region_coordinates(50, region5);
 triangle_mesh_parameters.add_region_coordinates(60, region6);
 triangle_mesh_parameters.add_region_coordinates(70, region7);

 My_mesh_pt = new RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters);

 Problem::mesh_pt() = My_mesh_pt;
 
 // Set error estimator for bulk mesh
 error_estimator_pt = new Z2ErrorEstimator;

 My_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;
/*
 // Set targets for spatial adaptivity
 My_mesh_pt->max_permitted_error()=0.0001;
 My_mesh_pt->min_permitted_error()=0.001;
 My_mesh_pt->max_element_size()=0.2;
 My_mesh_pt->min_element_size()=0.001;
*/

 complete_problem_setup();

}

template<class ELEMENT>
void PoissonProblem<ELEMENT>::complete_problem_setup() {

  // Set the boundary conditions for this problem: All nodes are
  // free by default -- just pin the ones that have Dirichlet conditions
  // here.
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
    {
      mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
    }
  }

  // Complete the build of all elements so they are fully functional

  //Find number of elements in mesh
  unsigned n_element = mesh_pt()->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for(unsigned i=0;i<n_element;i++)
   {
    // Upcast from GeneralElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

    //Set the source function pointer
    el_pt->source_fct_pt() = Source_fct_pt;
   }

  // Setup equation numbering scheme
  cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;


   apply_boundary_conditions();

}

template<class ELEMENT>
void PoissonProblem<ELEMENT>::apply_boundary_conditions()
{

//Loop over the boundaries
unsigned num_bound = mesh_pt()->nboundary();
for(unsigned ibound=0;ibound<num_bound;ibound++)
 {
  // Loop over the nodes on boundary
  unsigned num_nod=mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
    double u;
    Vector<double> x(2);
    x[0]=nod_pt->x(0);
    x[1]=nod_pt->x(1);
    TanhSolnForPoisson::get_exact_u(x,u);
    nod_pt->set_value(0,u);
    //nod_pt->set_value(0,0.0);
   }
 }

}

//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info,
                                           ofstream &trace_file)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=2;

 // Output regions
 unsigned nreg=mesh_pt()->nregion();
 oomph_info << "nreg: " << nreg << std::endl;
 for (unsigned r=0;r<nreg;r++)
  {
   double att=mesh_pt()->region_attribute(r);
   unsigned region_id = (unsigned)att;
   oomph_info << "attribute in reg: " << att << std::endl;

   sprintf(filename,"%s/region%i_%i.dat",doc_info.directory().c_str(),
           region_id,doc_info.number());   
   some_file.open(filename);
   
   unsigned nel=mesh_pt()->nregion_element(region_id);
   oomph_info << "nel in reg: " << nel << std::endl;
  
   for (unsigned e=0;e<nel;e++)
    {     
     FiniteElement* el_pt=mesh_pt()->region_element_pt(region_id,e);
     el_pt->output(some_file,npts);
    }
   some_file.close();
  }

 // Check that all nodes on boundaries are boundary nodes
 unsigned nb=mesh_pt()->nboundary();
 for (unsigned b=0;b<nb;b++)
  {
   unsigned nnod=mesh_pt()->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt=mesh_pt()->boundary_node_pt(b,j);
     BoundaryNodeBase* bnod_pt=dynamic_cast<BoundaryNodeBase*>(nod_pt);
     if (bnod_pt==0)
      {
       oomph_info << "Cast failed\n";
      }
     else
      {
       //oomph_info << "Cast succeeded\n";
      }
    }
  }

 // Output boundaries
 //------------------
 sprintf(filename,"%s/boundaries.dat",doc_info.directory().c_str());
 some_file.open(filename);
 mesh_pt()->output_boundaries(some_file);
 some_file.close();


 // Output solution
 //----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();


 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
 some_file.close();


 // Doc error
 //----------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,TanhSolnForPoisson::get_exact_u,
                          error,norm); 
 some_file.close();
 cout << "error: " << sqrt(error) << std::endl; 
 cout << "norm : " << sqrt(norm) << std::endl << std::endl;

 // Output the norm
 trace_file << sqrt(norm) << std::endl;

 }

 



//========================================================================
/// Demonstrate how to solve Poisson problem
//========================================================================
int main(int argc, char* argv[])
{

 // initialise MPI
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Sample point container
 CommandLineArgs::specify_command_line_flag("--non_ref_bin");
 CommandLineArgs::specify_command_line_flag("--ref_bin");
#ifdef OOMPH_HAS_CGAL
 CommandLineArgs::specify_command_line_flag("--cgal");
#endif

 // Parse command line
 CommandLineArgs::parse_and_assign();
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 const unsigned max_adapt = 3;

 // Only set one!
 unsigned count=0;
 if (CommandLineArgs::command_line_flag_has_been_set("--non_ref_bin"))
  {
   count++;
   MeshAsGeomObject_Helper::Default_sample_point_container_version=
    UseNonRefineableBinArray;
  }
 if (CommandLineArgs::command_line_flag_has_been_set("--ref_bin"))
  {
   count++;
   MeshAsGeomObject_Helper::Default_sample_point_container_version=
    UseRefineableBinArray;
  }

#ifdef OOMPH_HAS_CGAL
 if (CommandLineArgs::command_line_flag_has_been_set("--cgal"))
  {
   count++;
   MeshAsGeomObject_Helper::Default_sample_point_container_version=
    UseCGALSamplePointContainer;
  }
#endif

 if (count>1)
  {
   std::ostringstream error_message;
   error_message
    << "Can only choose one of --non_ref_bin, --ref_bin or --cgal!";
   throw OomphLibError(error_message.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");

 // Do the problem with quadratic elements
 //---------------------------------------
 {
  cout << std::endl  << "Quadratic elements" << std::endl;
  cout <<               "===================" << std::endl << std::endl;

  //Set up the problem
  PoissonProblem<ProjectablePoissonElement<TPoissonElement<2,3> > >
   problem(&TanhSolnForPoisson::get_source);

  // Solve the problem
  problem.newton_solve(max_adapt);

  // Open trace file
  char trace_filename[100];
  sprintf(trace_filename,"%s/trace.dat", 
          doc_info.directory().c_str());

  /// Trace file to document norm of solution
  ofstream trace_file;
  trace_file.open(trace_filename);

  //Output solution
  problem.doc_solution(doc_info, trace_file);

  //Increment counter for solutions
  doc_info.number()++;

  // Close the trace file
  trace_file.close();
  
 }

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

}



