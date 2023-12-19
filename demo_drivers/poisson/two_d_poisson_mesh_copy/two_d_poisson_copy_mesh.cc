// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Driver for a simple 2D poisson problem

// Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/simple_rectangular_tri_mesh.h"
#include "meshes/triangle_mesh.h"

using namespace std;

using namespace oomph;

//===== start_of_namespace=============================================
/// Namespace for exact solution for Poisson equation with "sharp step"
//=====================================================================
namespace TanhSolnForPoisson
{
  /// Parameter for steepness of "step"
  double Alpha = 1.0;

  /// Parameter for angle Phi of "step"
  double TanPhi = 0.0;

  /// Exact solution as a Vector
  void get_exact_u(const Vector<double>& x, Vector<double>& u)
  {
    u[0] = tanh(1.0 - Alpha * (TanPhi * x[0] - x[1]));
  }

  /// Source function required to make the solution above an exact solution
  void source_function(const Vector<double>& x, double& source)
  {
    source = 2.0 * tanh(-1.0 + Alpha * (TanPhi * x[0] - x[1])) *
               (1.0 - pow(tanh(-1.0 + Alpha * (TanPhi * x[0] - x[1])), 2.0)) *
               Alpha * Alpha * TanPhi * TanPhi +
             2.0 * tanh(-1.0 + Alpha * (TanPhi * x[0] - x[1])) *
               (1.0 - pow(tanh(-1.0 + Alpha * (TanPhi * x[0] - x[1])), 2.0)) *
               Alpha * Alpha;
  }

} // namespace TanhSolnForPoisson

//====== start_of_problem_class=======================================
/// Refineable 2D Poisson problem on rectangular domain, discretised
/// with 2D QPoisson elements. The specific type of element is
/// specified via the template parameter.
//====================================================================
template<class ELEMENT>
class RefineablePoissonProblem : public Problem
{
private:
  enum
  {
    Bottom_boundary,
    Right_boundary,
    Top_boundary,
    Left_boundary,
  };
  RefineableTriangleMesh<ELEMENT>* My_mesh_pt;

  /// Pointer to source function
  PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

public:
  /// Constructor: Pass pointer to source function
  RefineablePoissonProblem(
    PoissonEquations<2>::PoissonSourceFctPt source_fct_pt,
    const unsigned& n_adapt);

  /// Destructor (empty)
  ~RefineablePoissonProblem() {}

  /// Update the problem specs before solve: Reset boundary conditions
  /// to the values from the exact solution.
  void actions_before_newton_solve();

  /// Update the problem after solve (empty)
  void actions_after_newton_solve() {}

  /// Doc the solution. DocInfo object stores flags/labels for where the
  /// output gets written to
  void doc_solution(DocInfo& doc_info);

  TriangleMesh<ELEMENT>* bulk_mesh_pt();
}; // end of problem class


//=====start_of_constructor===============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT>
RefineablePoissonProblem<ELEMENT>::RefineablePoissonProblem(
  PoissonEquations<2>::PoissonSourceFctPt source_fct_pt,
  const unsigned& n_adapt)
  : Source_fct_pt(source_fct_pt)
{
  // Setup mesh
  // Domain length in x-direction
  double l_x = 1.0;

  // Domain length in y-direction
  double l_y = 2.0;

  Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);

  // Each polyline only has two vertices -- provide storage for their
  // coordinates
  Vector<Vector<double>> vertex_coord(2);
  for (unsigned i = 0; i < 2; i++)
  {
    vertex_coord[i].resize(2);
  }

  // First polyline: Free_surface_boundary_id
  vertex_coord[0][0] = 0.0;
  vertex_coord[0][1] = 0.0;
  vertex_coord[1][0] = l_x;
  vertex_coord[1][1] = 0.0;

  // Build the 1st boundary polyline
  boundary_polyline_pt[0] =
    new TriangleMeshPolyLine(vertex_coord, Bottom_boundary);

  // Second boundary polyline: Outer wall with slip
  vertex_coord[0][0] = vertex_coord[1][0];
  vertex_coord[0][1] = vertex_coord[1][1];
  vertex_coord[1][0] = l_x;
  vertex_coord[1][1] = l_y;

  // Build the 2nd boundary polyline
  boundary_polyline_pt[1] =
    new TriangleMeshPolyLine(vertex_coord, Right_boundary);

  // Third boundary polyline: Outflow
  vertex_coord[0][0] = vertex_coord[1][0];
  vertex_coord[0][1] = vertex_coord[1][1];
  vertex_coord[1][0] = 0.0;
  vertex_coord[1][1] = l_y;

  // Build the 3rd boundary polyline
  boundary_polyline_pt[2] =
    new TriangleMeshPolyLine(vertex_coord, Top_boundary);

  // Fourth boundary polyline: Bottom wall
  vertex_coord[0][0] = vertex_coord[1][0];
  vertex_coord[0][1] = vertex_coord[1][1];
  vertex_coord[1][0] = 0.0;
  vertex_coord[1][1] = 0.0;

  // Build the 4th boundary polyline
  boundary_polyline_pt[3] =
    new TriangleMeshPolyLine(vertex_coord, Left_boundary);

  // Build and assign mesh
  TriangleMeshClosedCurve* outer_boundary =
    new TriangleMeshPolygon(boundary_polyline_pt);
  TriangleMeshParameters triangle_mesh_parameters(outer_boundary);
  My_mesh_pt = new RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters);

  My_mesh_pt->spatial_error_estimator_pt() = new Z2ErrorEstimator;
  add_sub_mesh(My_mesh_pt);
  build_global_mesh();

  for (unsigned n = 0; n < n_adapt; n++)
  {
    adapt();
  }

  // Set the boundary conditions for this problem: All nodes are
  // free by default -- only need to pin the ones that have Dirichlet
  // conditions here.
  unsigned n_bound = mesh_pt()->nboundary();
  for (unsigned i = 0; i < n_bound; i++)
  {
    unsigned n_node = mesh_pt()->nboundary_node(i);
    for (unsigned n = 0; n < n_node; n++)
    {
      mesh_pt()->boundary_node_pt(i, n)->pin(0);
    }
  }

  // Complete the build of all elements so they are fully functional

  // Loop over the elements to set up element-specific
  // things that cannot be handled by the (argument-free!) ELEMENT
  // constructor: Pass pointer to source function
  unsigned n_element = mesh_pt()->nelement();
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralsedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

    // Set the source function pointer
    el_pt->source_fct_pt() = Source_fct_pt;
  }


  // Setup equation numbering scheme
  cout << "Number of equations: " << assign_eqn_numbers() << std::endl;

} // end of constructor


//=================================start_of_actions_before_newton_solve===
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void RefineablePoissonProblem<ELEMENT>::actions_before_newton_solve()
{
  // How many boundaries are there?
  unsigned n_bound = mesh_pt()->nboundary();

  // Loop over the boundaries
  for (unsigned i = 0; i < n_bound; i++)
  {
    // How many nodes are there on this boundary?
    unsigned n_node = mesh_pt()->nboundary_node(i);

    // Loop over the nodes on boundary
    for (unsigned n = 0; n < n_node; n++)
    {
      // Get pointer to node
      Node* nod_pt = mesh_pt()->boundary_node_pt(i, n);

      // Extract nodal coordinates from node:
      Vector<double> x(2);
      x[0] = nod_pt->x(0);
      x[1] = nod_pt->x(1);

      // Compute the value of the exact solution at the nodal point
      Vector<double> u(1);
      TanhSolnForPoisson::get_exact_u(x, u);

      // Assign the value to the one (and only) nodal value at this node
      nod_pt->set_value(0, u[0]);
    }
  }
} // end of actions before solve


//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void RefineablePoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  ofstream some_file;
  char filename[100];

  // Number of plot points: npts x npts
  unsigned npts = 5;

  // Output solution
  //-----------------
  sprintf(
    filename, "%s/soln%i.dat", doc_info.directory().c_str(), doc_info.number());
  some_file.open(filename);
  mesh_pt()->output(some_file, npts);
  some_file.close();


  // Output exact solution
  //----------------------
  sprintf(filename,
          "%s/exact_soln%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  mesh_pt()->output_fct(some_file, npts, TanhSolnForPoisson::get_exact_u);
  some_file.close();

  // Doc error and return of the square of the L2 error
  //---------------------------------------------------
  double error, norm;
  sprintf(filename,
          "%s/error%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  mesh_pt()->compute_error(
    some_file, TanhSolnForPoisson::get_exact_u, error, norm);
  some_file.close();

  // Doc L2 error and norm of solution
  cout << "\nNorm of error   : " << sqrt(error) << std::endl;
  cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

} // end of doc

template<class ELEMENT>
TriangleMesh<ELEMENT>* RefineablePoissonProblem<ELEMENT>::bulk_mesh_pt()
{
  return dynamic_cast<TriangleMesh<ELEMENT>*>(My_mesh_pt);
}

//====== start_of_problem_class=======================================
/// 2D Poisson problem on rectangular domain, discretised with
/// 2D QPoisson elements. The specific type of element is
/// specified via the template parameter.
//====================================================================
template<class ELEMENT>
class MyPoissonProblem : public Problem
{
private:
  enum
  {
    Bottom_boundary,
    Right_boundary,
    Top_boundary,
    Left_boundary,
  };
  TriangleMesh<ELEMENT>* My_mesh_pt;

  /// Pointer to source function
  PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

public:
  /// Constructor: Pass pointer to source function
  template<class ORIGINAL_ELEMENT>
  MyPoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt,
                   TriangleMesh<ORIGINAL_ELEMENT>* passed_mesh_pt);

  /// Destructor (empty)
  ~MyPoissonProblem() {}

  /// Update the problem specs before solve: Reset boundary conditions
  /// to the values from the exact solution.
  void actions_before_newton_solve();

  /// Update the problem after solve (empty)
  void actions_after_newton_solve() {}

  /// Doc the solution. DocInfo object stores flags/labels for where the
  /// output gets written to
  void doc_solution(DocInfo& doc_info);
}; // end of problem class


//=====start_of_constructor===============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT>
template<class ORIGINAL_ELEMENT>
MyPoissonProblem<ELEMENT>::MyPoissonProblem(
  PoissonEquations<2>::PoissonSourceFctPt source_fct_pt,
  TriangleMesh<ORIGINAL_ELEMENT>* passed_mesh_pt)
  : Source_fct_pt(source_fct_pt)
{
  // Setup mesh
  if (passed_mesh_pt)
  {
    My_mesh_pt = new TriangleMesh<ELEMENT>();
    My_mesh_pt->build_from_another_mesh(passed_mesh_pt);
  }
  else
  {
    // Domain length in x-direction
    double l_x = 1.0;

    // Domain length in y-direction
    double l_y = 2.0;

    Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);

    // Each polyline only has two vertices -- provide storage for their
    // coordinates
    Vector<Vector<double>> vertex_coord(2);
    for (unsigned i = 0; i < 2; i++)
    {
      vertex_coord[i].resize(2);
    }

    // First polyline: Free_surface_boundary_id
    vertex_coord[0][0] = 0.0;
    vertex_coord[0][1] = 0.0;
    vertex_coord[1][0] = l_x;
    vertex_coord[1][1] = 0.0;

    // Build the 1st boundary polyline
    boundary_polyline_pt[0] =
      new TriangleMeshPolyLine(vertex_coord, Bottom_boundary);

    // Second boundary polyline: Outer wall with slip
    vertex_coord[0][0] = vertex_coord[1][0];
    vertex_coord[0][1] = vertex_coord[1][1];
    vertex_coord[1][0] = l_x;
    vertex_coord[1][1] = l_y;

    // Build the 2nd boundary polyline
    boundary_polyline_pt[1] =
      new TriangleMeshPolyLine(vertex_coord, Right_boundary);

    // Third boundary polyline: Outflow
    vertex_coord[0][0] = vertex_coord[1][0];
    vertex_coord[0][1] = vertex_coord[1][1];
    vertex_coord[1][0] = 0.0;
    vertex_coord[1][1] = l_y;

    // Build the 3rd boundary polyline
    boundary_polyline_pt[2] =
      new TriangleMeshPolyLine(vertex_coord, Top_boundary);

    // Fourth boundary polyline: Bottom wall
    vertex_coord[0][0] = vertex_coord[1][0];
    vertex_coord[0][1] = vertex_coord[1][1];
    vertex_coord[1][0] = 0.0;
    vertex_coord[1][1] = 0.0;

    // Build the 4th boundary polyline
    boundary_polyline_pt[3] =
      new TriangleMeshPolyLine(vertex_coord, Left_boundary);

    // Build and assign mesh
    TriangleMeshClosedCurve* outer_boundary =
      new TriangleMeshPolygon(boundary_polyline_pt);
    TriangleMeshParameters triangle_mesh_parameters(outer_boundary);
    My_mesh_pt = new TriangleMesh<ELEMENT>(triangle_mesh_parameters);
  }
  add_sub_mesh(My_mesh_pt);
  build_global_mesh();

  // Set the boundary conditions for this problem: All nodes are
  // free by default -- only need to pin the ones that have Dirichlet
  // conditions here.
  unsigned n_bound = mesh_pt()->nboundary();
  for (unsigned i = 0; i < n_bound; i++)
  {
    unsigned n_node = mesh_pt()->nboundary_node(i);
    for (unsigned n = 0; n < n_node; n++)
    {
      mesh_pt()->boundary_node_pt(i, n)->pin(0);
    }
  }

  // Complete the build of all elements so they are fully functional

  // Loop over the elements to set up element-specific
  // things that cannot be handled by the (argument-free!) ELEMENT
  // constructor: Pass pointer to source function
  unsigned n_element = mesh_pt()->nelement();
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralsedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

    // Set the source function pointer
    el_pt->source_fct_pt() = Source_fct_pt;
  }


  // Setup equation numbering scheme
  cout << "Number of equations: " << assign_eqn_numbers() << std::endl;

} // end of constructor


//========================================start_of_actions_before_newton_solve===
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void MyPoissonProblem<ELEMENT>::actions_before_newton_solve()
{
  // How many boundaries are there?
  unsigned n_bound = mesh_pt()->nboundary();

  // Loop over the boundaries
  for (unsigned i = 0; i < n_bound; i++)
  {
    // How many nodes are there on this boundary?
    unsigned n_node = mesh_pt()->nboundary_node(i);

    // Loop over the nodes on boundary
    for (unsigned n = 0; n < n_node; n++)
    {
      // Get pointer to node
      Node* nod_pt = mesh_pt()->boundary_node_pt(i, n);

      // Extract nodal coordinates from node:
      Vector<double> x(2);
      x[0] = nod_pt->x(0);
      x[1] = nod_pt->x(1);

      // Compute the value of the exact solution at the nodal point
      Vector<double> u(1);
      TanhSolnForPoisson::get_exact_u(x, u);

      // Assign the value to the one (and only) nodal value at this node
      nod_pt->set_value(0, u[0]);
    }
  }
} // end of actions before solve


//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void MyPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  ofstream some_file;
  char filename[100];

  // Number of plot points: npts x npts
  unsigned npts = 5;

  // Output solution
  //-----------------
  sprintf(
    filename, "%s/soln%i.dat", doc_info.directory().c_str(), doc_info.number());
  some_file.open(filename);
  mesh_pt()->output(some_file, npts);
  some_file.close();


  // Output exact solution
  //----------------------
  sprintf(filename,
          "%s/exact_soln%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  mesh_pt()->output_fct(some_file, npts, TanhSolnForPoisson::get_exact_u);
  some_file.close();

  // Doc error and return of the square of the L2 error
  //---------------------------------------------------
  double error, norm;
  sprintf(filename,
          "%s/error%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  mesh_pt()->compute_error(
    some_file, TanhSolnForPoisson::get_exact_u, error, norm);
  some_file.close();

  // Doc L2 error and norm of solution
  cout << "\nNorm of error   : " << sqrt(error) << std::endl;
  cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

} // end of doc


//===== start_of_main=====================================================
/// Driver code for 2D Poisson problem
//========================================================================
int main()
{
  // Set up parameters
  TanhSolnForPoisson::TanPhi = 1.0;
  TanhSolnForPoisson::Alpha = 1.0;

  // Create initial problem
  RefineablePoissonProblem<ProjectablePoissonElement<TPoissonElement<2, 3>>>
    problem(&TanhSolnForPoisson::source_function, 1);

  // Solve initial problem
  problem.newton_solve();

  // Document initial problem
  DocInfo doc_info;
  doc_info.set_directory("RESLT");
  doc_info.number() = 0;
  problem.doc_solution(doc_info);
  doc_info.number()++;

  TriangleMesh<ProjectablePoissonElement<TPoissonElement<2, 3>>>* old_mesh_pt =
    problem.bulk_mesh_pt();

  // Create new problem from initial problem's mesh
  MyPoissonProblem<TPoissonElement<2, 3>> new_problem(
    &TanhSolnForPoisson::source_function, old_mesh_pt);

  new_problem.doc_solution(doc_info);
  doc_info.number()++;

  // Solve new problem
  new_problem.newton_solve();

  // Document the result
  new_problem.doc_solution(doc_info);
  doc_info.number()++;

} // end of main
