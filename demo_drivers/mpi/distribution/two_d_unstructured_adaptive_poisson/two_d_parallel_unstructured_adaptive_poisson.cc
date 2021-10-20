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
#include <fenv.h> 

// Driver code for the solution of the possion equation in different
// domain geometries. The unstructured mesh is adapted in parallel.

// Check the documentation of the results and the name of the
// files. When working in parallel, each processor generates an output
// file with its corresponding "part" of the work

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

// ==================================================================
// Data to do the tests
// ==================================================================
namespace TestArguments
{
  unsigned Domain_configuration = 1; // Choose the domain geometry
  double Element_size = 1.0e-2; // same as default
  unsigned Max_adapt = 3;
  double Max_permitted_error = 1.0e-3; // same as default
  double Min_permitted_error = 1.0e-5; // same as default
  double Max_element_size = Element_size;
  double Min_element_size = 1.0e-14; // allow small elements if
                                     // required
  unsigned Load_balance = 0; // change this to enable load balance (or
                             // indicate it in the arguments when
                             // running this demo driver)
  
  /// Folder where is stored the distribution file indicating the
  /// elements that each processor is in charge
  std::string Folder_distribution_file="";
}

//====================================================================
/// Poisson problem
//====================================================================

// Poisson problem
template<class ELEMENT> 
class PoissonProblem : public Problem
{

public:


 ///  Constructor
 PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt,
                const double element_size, const unsigned domain_configuration);

 /// Destructor (empty)
 ~PoissonProblem(){}

 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve()
  {
   apply_boundary_conditions();
  }

 /// Update the problem specs before solve (empty)
 void actions_after_newton_solve() {}

 /// Apply actions before adapt (empty)
 void actions_before_adapt() {}
 
 /// Apply actions after adapt
 void actions_after_adapt()
  {
   complete_problem_setup();
  }
 
 /// Apply actions before distribute (same as actions before adapt)
 void actions_before_distribute()
  {
   actions_before_adapt();
  }
 
 /// Apply actions after distribute (same as actions after adapt)
 void actions_after_distribute()
  {
   // BEGIN --- Only required if calling load_balance() ----------------
    
   // After calling load_balance() we need to re-establish the
   // pointers of the unstructured meshes (ONLY). This is because the
   // meshes created by the build_mesh() called from the
   // load_balance() method (to replicate the behaviour of the
   // structured load balance case) were deleted as part of the load
   // balancing process.
   
    // Get the mesh pointer from mesh_pt() and set it on the
    // corresponding unstructured mesh pointer.
   Bulk_mesh_pt = mesh_pt();
   
   // Make the mesh_pt() from problem point to the corresponding
   // mesh. This is similar to what we do in the build_mesh() method
   Problem::mesh_pt() = Bulk_mesh_pt;
   
    // END --- Only required if calling load_balance() ----------------
   
   // Then we can safely call actions_after_adapt(). If we do not do
   // the previous two steps after calling load_balance() then the
   // pointer 'Bulk_mesh_pt' would be pointing to NULL, thus causing an
   // error.
   actions_after_adapt();
  }
 
 // Access funtion for the specific mesh
 RefineableTriangleMesh<ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RefineableTriangleMesh<ELEMENT>*>(Problem::mesh_pt());
  }
 
 /// Doc the solution (when working in parallel each processor creates
 /// output files, thus we must be careful with the naming of the
 /// output files)
 void doc_solution(DocInfo& doc_info, ofstream &trace_file);
 
 /// Required function ONLY to perform load balancing. However, it is
 /// a good practice to implement this function and build the mesh in
 /// there
 void build_mesh();
 
 /// Read the distribution data
 void read_custom_distribution_from_file(Vector<unsigned> &output_distribution);
  
 /// Saves the custom distribution to file
 void save_custom_distribution_to_file(Vector<unsigned> &input_distribution);
  
private:
 
 // Apply boundary conditions (used at the end of the constructor and
 // in actions after adapt() )
 void apply_boundary_conditions();
 
 // Wraps the setting of fully funtional elements in the mesh and the
 // application of boundary conditions (useful in actions after mesh
 // adaptations)
 void complete_problem_setup();
 
 // The adaptive unstructured mesh
 RefineableTriangleMesh<ELEMENT> *Bulk_mesh_pt;
 
 // Erros estimator for the unstructured adaptive mesh
 Z2ErrorEstimator* Error_estimator_pt;
 
 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;
 
 // The TriangleMeshParameters object is defined as a Problem object.
 // This is to get access to the domain configuration when calling
 // build_mesh() from load_balance()
 TriangleMeshParameters Triangle_mesh_parameters;
  
  // Set the element size
  double Element_size;
  
  // The domain configuration
  unsigned Domain_configuration;
  
  // The methods that define the domain's geometry
  void square_domain();
  void half_circle_domain();
  void half_circle_domain_with_internal_boundaries();
  void complex_domain_with_holes();
  
};

//========================================================================
/// Constructor for Poisson problem
//========================================================================

template<class ELEMENT>
PoissonProblem<ELEMENT>::
PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt,
               const double element_size, const unsigned domain_configuration)
  : Source_fct_pt(source_fct_pt), 
    Element_size(element_size),
    Domain_configuration(domain_configuration)
{
 
 // Setup parameters for exact tanh solution
 
 // Steepness of step
  TanhSolnForPoisson::Alpha=10.0;
 
 // Orientation of step
 TanhSolnForPoisson::Beta=1.0;
 
 if (Domain_configuration == 1)
   {
     square_domain();
   }
 else if (Domain_configuration == 2)
   {
     half_circle_domain();
   }
 else if (Domain_configuration == 3)
   {
     half_circle_domain_with_internal_boundaries();
   }
 else if (Domain_configuration == 4)
   {
     complex_domain_with_holes();
   }
 else
   {
     std::ostringstream error_message;
     error_message
       << "The choosen domain configuration is not implemented\n"
       << "You chose ("<<Domain_configuration<<")";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
   }
 
 // At this point the Triangle_mesh_parameters object stores the
 // domain geometry. We next specify the element size
 
 // Set the element size
 Triangle_mesh_parameters.element_area() = Element_size;
 
 // Create the error estimator object for the bulk mesh
 Error_estimator_pt = new Z2ErrorEstimator;
 
 // Build the mesh
 build_mesh();
 
}

//========================================================================
/// Required function to perform load balancing, create the mesh
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::build_mesh()
{
 // Create the mesh (using the geometry stored in the
 // Triangle_mesh_parameters object)
 Bulk_mesh_pt = new RefineableTriangleMesh<ELEMENT>(Triangle_mesh_parameters);
 
 // ... and set the just created mesh as the mesh in the problem object
 Problem::mesh_pt() = Bulk_mesh_pt;
 
 // Set the error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt() = Error_estimator_pt;
  
 // Establish the max. and min. permitted errors for adaptation
 Bulk_mesh_pt->max_permitted_error() = 
   TestArguments::Max_permitted_error; // refine if error is larger
                                       // than this
 
 Bulk_mesh_pt->min_permitted_error() = 
   TestArguments::Min_permitted_error; // unrefine if error is smaller
                                       // than this
 
 // Set the maximum element size
 Bulk_mesh_pt->max_element_size() = TestArguments::Element_size;
 
 // Set the minimum element size
 Bulk_mesh_pt->min_element_size() = TestArguments::Min_element_size;
 
 // Make all the elements fully funtional and set boundary conditions
 complete_problem_setup();
 
}

//========================================================================
/// Builds the square domain mesh
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::square_domain()
{
  // Square domain
  unsigned boundary_id;
  const double length = 1.0;
  
  // ---------------------------------------------------------------------
  // >> Boundary 0
  const unsigned num_vertices_b0 = 2;
  
  Vector<Vector <double> > vertices(num_vertices_b0);
  
  for (unsigned i = 0; i < num_vertices_b0; i++)
    {
      vertices[i].resize(2);
    }
  
  vertices[0][0] = 0;
  vertices[0][1] = 0;
  
  vertices[1][0] = 0;
  vertices[1][1] = length;
  
  boundary_id = 0;
  TriangleMeshPolyLine *boundary0_pt =
    new TriangleMeshPolyLine(vertices, boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Boundary 1
  const unsigned num_vertices_b1 = 2;
  
  vertices.resize(num_vertices_b1);
  for (unsigned i = 0; i < num_vertices_b1; i++)
    {
      vertices[i].resize(2);
    }
  
  vertices[0][0] = 0;
  vertices[0][1] = length;
  
  vertices[1][0] = length;
  vertices[1][1] = length;
  
  boundary_id = 1;
  TriangleMeshPolyLine *boundary1_pt =
    new TriangleMeshPolyLine(vertices, boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Boundary 2
  const unsigned num_vertices_b2 = 2;
  
  vertices.resize(num_vertices_b2);
  for (unsigned i = 0; i < num_vertices_b2; i++)
    {
      vertices[i].resize(2);
    }
  
  vertices[0][0] = length;
  vertices[0][1] = length;
  
  vertices[1][0] = length;
  vertices[1][1] = 0;
    
  boundary_id = 2;
  TriangleMeshPolyLine *boundary2_pt = 
    new TriangleMeshPolyLine(vertices, boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Boundary 3
  const unsigned num_vertices_b3 = 2;
  
  vertices.resize(num_vertices_b3);
  for (unsigned i = 0; i < num_vertices_b3; i++)
    {
      vertices[i].resize(2);
    }
  
  vertices[0][0] = length;
  vertices[0][1] = 0;
    
  vertices[1][0] = 0;
  vertices[1][1] = 0;
  
  boundary_id = 3;
  TriangleMeshPolyLine *boundary3_pt = 
    new TriangleMeshPolyLine(vertices, boundary_id);

  // ---------------------------------------------------------------------
  // >> Building the OUTER BOUNDARY
  
  // >> Setting up the domain with PolyLines
  
  Vector<TriangleMeshCurveSection*> outer_boundary_polylines_pt(4);
  
  outer_boundary_polylines_pt[0] = boundary0_pt;
  outer_boundary_polylines_pt[1] = boundary1_pt;
  outer_boundary_polylines_pt[2] = boundary2_pt;
  outer_boundary_polylines_pt[3] = boundary3_pt;

  // We have only one outer boundary (one outer polygon)
  Vector<TriangleMeshClosedCurve *> outer_boundary_pt(1);
  outer_boundary_pt[0] = 
    new TriangleMeshPolygon(outer_boundary_polylines_pt);
  
  // ---------------------------------------------------------------------
  // Store the domain configuration in the Triangle_mesh_parameters
  // object
  Triangle_mesh_parameters.outer_boundary_pt() = outer_boundary_pt; 
  
}

//========================================================================
/// Builds the half circle domain
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::half_circle_domain()
{
  unsigned boundary_id;
  
  // ---------------------------------------------------------------------
  // >> Boundary 0
  const unsigned num_vertices_b0 = 2;
  
  Vector<Vector <double> > vertices_b0(num_vertices_b0);
  
  for (unsigned i = 0; i < num_vertices_b0; i++)
    vertices_b0[i].resize(2);
  
  vertices_b0[0][0] = 0;
  vertices_b0[0][1] = 0;

  vertices_b0[1][0] = 2;
  vertices_b0[1][1] = 0;
  
  boundary_id = 0;
  TriangleMeshPolyLine *boundary0_pt =
    new TriangleMeshPolyLine(vertices_b0, boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Boundary 1
  const unsigned num_segments_b1 = 10;
  double x_centre_b1 = 1.0;
  double y_centre_b1 = 0.0;
  double r_circle_b1 = 1.0;
  
  Circle *boundary_circle_b1_pt = 
    new Circle(x_centre_b1, y_centre_b1, r_circle_b1);
  
  double z_start_b1 = 0.0;
  double z_end_b1 = MathematicalConstants::Pi;
  
  boundary_id = 1;
  TriangleMeshCurviLine *boundary1_pt = 
    new TriangleMeshCurviLine(boundary_circle_b1_pt,
                              z_start_b1,
                              z_end_b1,
                              num_segments_b1,
                              boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Building the OUTER BOUNDARY
  
  // >> Setting up the domain with PolyLines
  
  Vector<TriangleMeshCurveSection*> outer_boundary_curve_section_pt(2);

  outer_boundary_curve_section_pt[0] = boundary0_pt;
  outer_boundary_curve_section_pt[1] = boundary1_pt;
  
  // We have only one outer boundary (one outer polygon)
  Vector<TriangleMeshClosedCurve *> outer_boundary_pt(1);
  outer_boundary_pt[0] = 
    new TriangleMeshClosedCurve(outer_boundary_curve_section_pt);
  
  // ---------------------------------------------------------------------
  // Store the domain configuration in the Triangle_mesh_parameters
  // object
  Triangle_mesh_parameters.outer_boundary_pt() = outer_boundary_pt;
  
}

//========================================================================
/// Builds the half circle with straight internal boundaries
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::half_circle_domain_with_internal_boundaries()
{
  unsigned boundary_id;
  
  // ---------------------------------------------------------------------
  // >> Boundary 0
  const unsigned num_vertices_b0 = 6;
  
  Vector<Vector <double> > vertices_b0(num_vertices_b0);
  
  for (unsigned i = 0; i < num_vertices_b0; i++)
    vertices_b0[i].resize(2);
  
  vertices_b0[0][0] = 0;
  vertices_b0[0][1] = 0;
  
  vertices_b0[1][0] = 7.5;
  vertices_b0[1][1] = 0;
  
  vertices_b0[2][0] = 7.5;
  vertices_b0[2][1] = 5;
  
  vertices_b0[3][0] = 12.5;
  vertices_b0[3][1] = 5;
 
  vertices_b0[4][0] = 12.5;
  vertices_b0[4][1] = 0;
  
  vertices_b0[5][0] = 20;
  vertices_b0[5][1] = 0;
  
  boundary_id = 0;
  TriangleMeshPolyLine *boundary0_pt =
    new TriangleMeshPolyLine(vertices_b0, boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Boundary 1
  const unsigned num_segments = 10;
  const double x_centre = 10.0;
  const double y_centre = 0.0;
  const double r_circle = 10.0;
  
  Circle *boundary_circle_pt = new Circle(x_centre, y_centre, r_circle);
 
  const double z_start = 0.0;
  const double z_end = MathematicalConstants::Pi;
  
  boundary_id = 1;
  TriangleMeshCurviLine *boundary1_pt = 
    new TriangleMeshCurviLine(boundary_circle_pt,
                              z_start,
                              z_end,
                              num_segments,
                              boundary_id);
  
  // ---------------------------------------------------------------------
  // ---------------------------------------------------------------------
  // Internal boundaries (open curves)
  // ---------------------------------------------------------------------
  // >> Boundary 2
  const unsigned num_vertices_b2 = 3;
  Vector<Vector<double> > verticesb2(num_vertices_b2);
  
  for (unsigned i = 0; i < num_vertices_b2; i++)
    verticesb2[i].resize(2);
  
  verticesb2[0][0] = 6;
  verticesb2[0][1] = 6;
  
  verticesb2[1][0] = 7.5;
  verticesb2[1][1] = 6;
  
  verticesb2[2][0] = 9;
  verticesb2[2][1] = 6;
  
  boundary_id = 2;
  TriangleMeshPolyLine *boundary2_pt =
    new TriangleMeshPolyLine(verticesb2, boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Boundary 3
  const unsigned num_vertices_b3 = 3;
  Vector<Vector<double> > verticesb3(num_vertices_b3);
  
  for (unsigned i = 0; i < num_vertices_b3; i++)
  verticesb3[i].resize(2);
  
  verticesb3[0][0] = 9;
  verticesb3[0][1] = 6;
  
  verticesb3[1][0] = 10;
  verticesb3[1][1] = 6;
 
  verticesb3[2][0] = 11;
  verticesb3[2][1] = 6;
  
  boundary_id = 3;
  TriangleMeshPolyLine *boundary3_pt =
    new TriangleMeshPolyLine(verticesb3, boundary_id);
 
  // ---------------------------------------------------------------------
  // >> Boundary 4
  const unsigned num_vertices_b4 = 3;
  Vector<Vector<double> > verticesb4(num_vertices_b4);
  
  for (unsigned i = 0; i < num_vertices_b4; i++)
    verticesb4[i].resize(2);
  
  verticesb4[0][0] = 11;
  verticesb4[0][1] = 6;
  
  verticesb4[1][0] = 12.5;
  verticesb4[1][1] = 6;
  
  verticesb4[2][0] = 14;
  verticesb4[2][1] = 6;
  
  boundary_id = 4;
  TriangleMeshPolyLine *boundary4_pt =
    new TriangleMeshPolyLine(verticesb4, boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Boundary 5
  const unsigned num_vertices_b5 = 3;
  Vector<Vector<double> > verticesb5(num_vertices_b5);
  
  for (unsigned i = 0; i < num_vertices_b5; i++)
    verticesb5[i].resize(2);
  
  verticesb5[0][0] = 7;
  verticesb5[0][1] = 7;
  
  verticesb5[1][0] = 10;
  verticesb5[1][1] = 7;
  
  verticesb5[2][0] = 13;
  verticesb5[2][1] = 7;
  
  boundary_id = 5;
  TriangleMeshPolyLine *boundary5_pt =
    new TriangleMeshPolyLine(verticesb5, boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Building the OUTER BOUNDARY
  Vector<TriangleMeshCurveSection*> outer_boundary_polylines_pt(2);
  
  outer_boundary_polylines_pt[0] = boundary0_pt;
  outer_boundary_polylines_pt[1] = boundary1_pt;
  
  // We have only one outer boundary (one outer polygon)
  Vector<TriangleMeshClosedCurve *> outer_boundary_pt(1);
  outer_boundary_pt[0] = 
    new TriangleMeshClosedCurve(outer_boundary_polylines_pt);
 
  // ---------------------------------------------------------------------
  // >> Building the OPEN BOUNDARIES
  
  // We have two open curves
  Vector<TriangleMeshOpenCurve *> open_curve_pt(2);
  
  // The first open curve is formed by the following boundaries
  Vector<TriangleMeshCurveSection*> open_curve_one_pt(3);
  open_curve_one_pt[0] = boundary2_pt;
  open_curve_one_pt[1] = boundary3_pt;
  open_curve_one_pt[2] = boundary4_pt;
 
  open_curve_pt[0] = new TriangleMeshOpenCurve(open_curve_one_pt);
  
  // The second open curve is formed by only one boundary
  Vector<TriangleMeshCurveSection*> open_curve_two_pt(1);
  open_curve_two_pt[0] = boundary5_pt;
  
  open_curve_pt[1] = new TriangleMeshOpenCurve(open_curve_two_pt);
  
  // ---------------------------------------------------------------------
  // Store the geometry of the domain in the Triangle_mesh_parameters
  // object
  
  // >> Create the TriangleMeshParameters object
  Triangle_mesh_parameters.outer_boundary_pt() = outer_boundary_pt;
  // Set the pointer to internal open boundaries
  Triangle_mesh_parameters.internal_open_curves_pt() = open_curve_pt;
  
}

//========================================================================
/// Builds a complex domain with holes
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::complex_domain_with_holes()
{
  unsigned boundary_id;
  
  // ---------------------------------------------------------------------
  // >> Boundary 0
  const unsigned num_vertices_b0 = 6;
  
  Vector<Vector <double> > vertices_b0(num_vertices_b0);
  
  for (unsigned i = 0; i < num_vertices_b0; i++)
    vertices_b0[i].resize(2);
  
  vertices_b0[0][0] = 0;
  vertices_b0[0][1] = 0;
  
  vertices_b0[1][0] = 7.5;
  vertices_b0[1][1] = 0;
  
  vertices_b0[2][0] = 7.5;
  vertices_b0[2][1] = 5;
  
  vertices_b0[3][0] = 12.5;
  vertices_b0[3][1] = 5;
  
  vertices_b0[4][0] = 12.5;
  vertices_b0[4][1] = 0;
  
  vertices_b0[5][0] = 20;
  vertices_b0[5][1] = 0;
  
  boundary_id = 0;
  TriangleMeshPolyLine *boundary0_pt =
    new TriangleMeshPolyLine(vertices_b0, boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Boundary 1
  const unsigned num_segments_b1 = 10;
  double x_centre_b1 = 10.0;
  double y_centre_b1 = 0.0;
  double r_circle_b1 = 10.0;
  
  Circle *boundary_circle_b1_pt = 
    new Circle(x_centre_b1, y_centre_b1, r_circle_b1);
  
  double z_start_b1 = 0.0;
  double z_end_b1 = MathematicalConstants::Pi;
  
  boundary_id = 1;
  TriangleMeshCurviLine *boundary1_pt = 
    new TriangleMeshCurviLine(boundary_circle_b1_pt,
                              z_start_b1,
                              z_end_b1,
                              num_segments_b1,
                              boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Building the OUTER BOUNDARY
  
  // >> Setting up the domain with PolyLines
  
  Vector<TriangleMeshCurveSection*> outer_boundary_curve_section_pt(2);

  outer_boundary_curve_section_pt[0] = boundary0_pt;
  outer_boundary_curve_section_pt[1] = boundary1_pt;
  
  // We have only one outer boundary (one outer polygon)
  Vector<TriangleMeshClosedCurve *> outer_boundary_pt(1);
  outer_boundary_pt[0] = 
    new TriangleMeshClosedCurve(outer_boundary_curve_section_pt);
  
  // Now create an internal polygon that will represent a hole
  // ---------------------------------------------------------------------
  // This hole has a similar geometry as the outer domain but rotated 90'
  
  // >> Boundary 2
  const unsigned num_vertices_b2 = 6;
  
  Vector<Vector <double> > vertices_b2(num_vertices_b2);
  
  for (unsigned i = 0; i < num_vertices_b2; i++)
    vertices_b2[i].resize(2);
  
  vertices_b2[0][0] = 14;
  vertices_b2[0][1] = 7;
  
  vertices_b2[1][0] = 14;
  vertices_b2[1][1] = 5;
  
  vertices_b2[2][0] = 15;
  vertices_b2[2][1] = 5;
  
  vertices_b2[3][0] = 15;
  vertices_b2[3][1] = 3;
  
  vertices_b2[4][0] = 14;
  vertices_b2[4][1] = 3;
  
  vertices_b2[5][0] = 14;
  vertices_b2[5][1] = 1;
  
  boundary_id = 2;
  TriangleMeshPolyLine *boundary2_pt =
    new TriangleMeshPolyLine(vertices_b2, boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Boundary 3
  const unsigned num_segments_b3 = 10;
  double x_centre_b3 = 14.0;
  double y_centre_b3 = 4.0;
  double r_circle_b3 = 3.0;
  
  Circle *boundary_circle_b3_pt=new Circle(x_centre_b3,y_centre_b3,r_circle_b3);
  
  double z_start_b3 = MathematicalConstants::Pi * (3.0 / 2.0);
  double z_end_b3 = MathematicalConstants::Pi * (5.0 / 2.0);
  
  boundary_id = 3;
  TriangleMeshCurviLine *boundary3_pt = 
    new TriangleMeshCurviLine(boundary_circle_b3_pt,
                              z_start_b3,
                              z_end_b3,
                              num_segments_b3,
                              boundary_id);
  
  // ---------------------------------------------------------------------
  // Another internal polygon (hole). A hole shaped as a square
  // >> Boundary 4
  const unsigned num_vertices_b4 = 3;
  
  Vector<Vector <double> > vertices_b4(num_vertices_b4);
  
  for (unsigned i = 0; i < num_vertices_b4; i++)
    vertices_b4[i].resize(2);
  
  vertices_b4[0][0] = 3;
  vertices_b4[0][1] = 1;
  
  vertices_b4[1][0] = 3;
  vertices_b4[1][1] = 2.5;
  
  vertices_b4[2][0] = 3;
  vertices_b4[2][1] = 4;
  
  boundary_id = 4;
  TriangleMeshPolyLine *boundary4_pt =
    new TriangleMeshPolyLine(vertices_b4, boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Boundary 5
  const unsigned num_vertices_b5 = 3;
  
  Vector<Vector <double> > vertices_b5(num_vertices_b5);
  
  for (unsigned i = 0; i < num_vertices_b5; i++)
    vertices_b5[i].resize(2);
  
  vertices_b5[0][0] = 3;
  vertices_b5[0][1] = 4;
  
  vertices_b5[1][0] = 4.5;
  vertices_b5[1][1] = 4;
  
  vertices_b5[2][0] = 6;
  vertices_b5[2][1] = 4;
  
  boundary_id = 5;
  TriangleMeshPolyLine *boundary5_pt =
    new TriangleMeshPolyLine(vertices_b5, boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Boundary 6
  const unsigned num_vertices_b6 = 3;
  
  Vector<Vector <double> > vertices_b6(num_vertices_b6);
  
  for (unsigned i = 0; i < num_vertices_b6; i++)
    vertices_b6[i].resize(2);
  
  vertices_b6[0][0] = 6;
  vertices_b6[0][1] = 4;
  
  vertices_b6[1][0] = 6;
  vertices_b6[1][1] = 2.5;
  
  vertices_b6[2][0] = 6;
  vertices_b6[2][1] = 1;
  
  boundary_id = 6;
  TriangleMeshPolyLine *boundary6_pt =
    new TriangleMeshPolyLine(vertices_b6, boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Boundary 7
  const unsigned num_vertices_b7 = 3;
  
  Vector<Vector <double> > vertices_b7(num_vertices_b7);
  
  for (unsigned i = 0; i < num_vertices_b7; i++)
    vertices_b7[i].resize(2);
  
  vertices_b7[0][0] = 6;
  vertices_b7[0][1] = 1;
  
  vertices_b7[1][0] = 4.5;
  vertices_b7[1][1] = 1;
  
  vertices_b7[2][0] = 3;
  vertices_b7[2][1] = 1;
  
  boundary_id = 7;
  TriangleMeshPolyLine *boundary7_pt =
    new TriangleMeshPolyLine(vertices_b7, boundary_id);
  
  // ---------------------------------------------------------------------
  // Another hole (internal polygon), created as a circle
  // >> Boundary 8
  const unsigned num_segments_b8 = 10;
  double x_centre_b8 = 10.0;
  double y_centre_b8 = 7.5;
  double r_circle_b8 = 1.0;
  
  Circle *boundary_circle_b8_pt=new Circle(x_centre_b8,y_centre_b8,r_circle_b8);
  
  double z_start_b8 = 0.0;
  double z_end_b8 = MathematicalConstants::Pi;
  
  boundary_id = 8;
  TriangleMeshCurviLine *boundary8_pt = 
    new TriangleMeshCurviLine(boundary_circle_b8_pt,
                              z_start_b8,
                              z_end_b8,
                              num_segments_b8,
                              boundary_id);
  
  // ---------------------------------------------------------------------
  // >> Boundary 9
  const unsigned num_segments_b9 = 10;
  double x_centre_b9 = 10.0;
  double y_centre_b9 = 7.5;
  double r_circle_b9 = 1.0;
  
  Circle *boundary_circle_b9_pt=new Circle(x_centre_b9,y_centre_b9,r_circle_b9);
  
  double z_start_b9 = MathematicalConstants::Pi;
  double z_end_b9 = 2.0 * MathematicalConstants::Pi;
  
  boundary_id = 9;
  TriangleMeshCurviLine *boundary9_pt = 
    new TriangleMeshCurviLine(boundary_circle_b9_pt,
                              z_start_b9,
                              z_end_b9,
                              num_segments_b9,
                              boundary_id);
  
  // ---------------------------------------------------------------------
  // ---------------------------------------------------------------------
  // >> Building the INTERNAL BOUNDARIES
  
  // Set the holes information (indicate the number of holes)
  const unsigned nholes = 3;
  
  Vector<TriangleMeshClosedCurve*> inner_boundaries_pt(nholes);
  Vector<Vector<TriangleMeshCurveSection*> > 
    inner_boundary_curve_section_pt(nholes);
  
  // -----------------
  // 1st hole
  inner_boundary_curve_section_pt[0].resize(2);
  
  inner_boundary_curve_section_pt[0][0] = boundary2_pt;
  inner_boundary_curve_section_pt[0][1] = boundary3_pt;
  
  // Set the coordinates to mark the holes (1st hole)
  Vector<double> hole1(2);
  hole1[0] = 16.0;
  hole1[1] = 4.0;
  
  // -----------------
  // 2nd hole
  inner_boundary_curve_section_pt[1].resize(4);
 
  inner_boundary_curve_section_pt[1][0] = boundary4_pt;
  inner_boundary_curve_section_pt[1][1] = boundary5_pt;
  inner_boundary_curve_section_pt[1][2] = boundary6_pt;
  inner_boundary_curve_section_pt[1][3] = boundary7_pt;
  
  // Set the coordinates to mark the holes (1st hole)
  Vector<double> hole2(2);
  hole2[0] = 4.5;
  hole2[1] = 2.5;
  
  // -----------------
  // 3rd hole
  inner_boundary_curve_section_pt[2].resize(2);
  
  inner_boundary_curve_section_pt[2][0] = boundary8_pt;
  inner_boundary_curve_section_pt[2][1] = boundary9_pt;
  
  // Set the coordinates to mark the holes (1st hole)
  Vector<double> hole3(2);
  hole3[0] = 10.0;
  hole3[1] = 7.5;
  
  // Store the holes geometries
  inner_boundaries_pt[0] = 
    new TriangleMeshClosedCurve(inner_boundary_curve_section_pt[0],hole1);
  
  inner_boundaries_pt[1] = 
    new TriangleMeshClosedCurve(inner_boundary_curve_section_pt[1],hole2);
  
  inner_boundaries_pt[2] = 
    new TriangleMeshClosedCurve(inner_boundary_curve_section_pt[2],hole3);
 
  // ---------------------------------------------------------------------
  // Store the geometry of the domain in the Triangle_mesh_parameters
  // object
  
  // ---------------------------------------------------------------------
  // Create the TriangleMeshParameters object and set the outer
  // boundary
  Triangle_mesh_parameters.outer_boundary_pt() = outer_boundary_pt;
  // Set the pointer to internal boundaries (holes)
  Triangle_mesh_parameters.internal_closed_curve_pt() = inner_boundaries_pt;
  
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
 
 // Disable the output of halo element so that in the output we only
 // see those elements that the processor is in charge
 mesh_pt()->disable_output_of_halo_elements();
 
 // Output solution (be careful with the naming of the output files,
 // each file can be identified by the processor id
 //  ----------------
 sprintf(filename,"%s/soln%i_proc%i.dat",
         doc_info.directory().c_str(),
         doc_info.number(),
         this->communicator_pt()->my_rank());
 
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
 // Output exact solution (output files in parallel)
 //----------------
 sprintf(filename,"%s/exact_soln%i_proc%i.dat",
         doc_info.directory().c_str(),
         doc_info.number(),
         this->communicator_pt()->my_rank());
 
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u);
 some_file.close();
 
 // Doc error (local to the processor, only consider the elements that
 // the processor is in charge)
 // -----------
 double error,norm;
 mesh_pt()->compute_error(some_file,TanhSolnForPoisson::get_exact_u,
                          error,norm);
 cout << "error: " << sqrt(error) << std::endl;
 cout << "norm : " << sqrt(norm) << std::endl << std::endl;
 
 // Because we want the global norm we need to send/receive
 // information to/from the other processors
 double norm_soln=0.0;
 mesh_pt()->compute_norm(norm_soln);  
#ifdef OOMPH_HAS_MPI
 if (mesh_pt()->is_mesh_distributed())
  {
   double norm_reduced = 0.0;
   MPI_Allreduce(&norm_soln, &norm_reduced, 1, MPI_DOUBLE, MPI_SUM, 
                 this->communicator_pt()->mpi_comm());
   
   cout << "norm reduced: " << sqrt(norm_reduced) << std::endl << std::endl;
   
   // Output the global norm which considers all processors
   trace_file << sqrt(norm_reduced) << std::endl;
  }
#endif
 
}

//========================================================================
/// Complete problem setup
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::complete_problem_setup()
{
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.
 
 const unsigned num_bound = mesh_pt()->nboundary();
 
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
 const unsigned n_element = mesh_pt()->nelement();
 
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

//========================================================================
/// Apply boundary conditions
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::apply_boundary_conditions()
{
 const unsigned num_bound = mesh_pt()->nboundary();
 
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary
   const unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
     double u;
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);
     TanhSolnForPoisson::get_exact_u(x,u);
     nod_pt->set_value(0,u);
    }
  }
}

//========================================================================
/// Read the distribution data
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::
read_custom_distribution_from_file(Vector<unsigned> &output_distribution)
{
 // Get the rank of the processor
 const unsigned my_rank = this->communicator_pt()->my_rank();
 
 // Pointer for the file
 FILE *file_pt;
 // Store the file name
 char file_name[500];
 sprintf(file_name, "%s/input_distribution_%i.dat", 
         TestArguments::Folder_distribution_file.c_str(),
         my_rank);
 
 oomph_info << "Read custom distribution from file: " << file_name 
            << std::endl;
 
 // Try to open the file
 file_pt = fopen(file_name, "r");
 if (file_pt == 0)
  {
   // Error, the file could not be opened
   std::ostringstream error_stream;
   error_stream << "Input file could not be opened: " 
                << file_name << std::endl;
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // Read the number of data
 unsigned n_data = 0;
 int iret = fscanf(file_pt, "NDATA:(%i)\n", &n_data);
 // Increase iret to avoid warning
 iret++;
 // Read the number of meshes
 unsigned tmp_n_mesh = 0;
 iret = fscanf(file_pt, "NMESHES:(%i)\n", &tmp_n_mesh);
 
 // Get the number of meshes in the current problem
 const unsigned n_sub_mesh = nsub_mesh();
 if (tmp_n_mesh != n_sub_mesh)
  {
   std::ostringstream error_message;
   error_message
     << "The number of sub-meshes in the input file is different from the "
     << "number of sub-meshes in the current problem\n"
     << "N.submeshes in problem: (" << n_sub_mesh << ")\n"
     << "N.submeshes in input file (" << tmp_n_mesh << ")";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
  }
 
 // Clear the output data and resize according to the number of data
 // to read
 output_distribution.clear();
 output_distribution.resize(n_data);
 
 // Finally, after all test have been successful loop and read over
 // the data
 for (unsigned e = 0; e < n_data; e++)
  {
   unsigned counter = 0;
   iret = fscanf(file_pt, "ELEMENT(%i):PROCESSOR(%i)\n",
          &counter, &output_distribution[e]);
  }
 
 // Close the file
 fclose(file_pt);
 
 oomph_info << "Read custom distribution from file: " << file_name 
            << " [DONE]" << std::endl;
  
}

//========================================================================
/// Saves the custom (input) distribution to file
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::
save_custom_distribution_to_file(Vector<unsigned> &input_distribution)
{
 // Get the rank of the processor
 const unsigned my_rank = this->communicator_pt()->my_rank();
 
 // Get the number of elements in the custom distribution
 const unsigned n_ele_custom_distribution = input_distribution.size();
 
 // Get the number of meshes in the current problem
 const unsigned n_sub_meshes = nsub_mesh();
 
 // Char for the output file name
 char file_name[500];
 // Set the name of the file
 sprintf(file_name, "DISTRIBUTION/input_distribution_%i.dat", my_rank);
 
 oomph_info << "Save custom distribution to file: " << file_name 
            << std::endl;
 
 // Open the file
 FILE* file_pt = fopen(file_name, "w");
 if (file_pt != NULL)
  {
   // Store the number of elements (DATA) in the mesh(es)
   fprintf(file_pt, "NDATA:(%i)\n", n_ele_custom_distribution);
   // Store the number of meshes
   fprintf(file_pt, "NMESHES:(%i)\n", n_sub_meshes);
   
   // Loop over the elements in the custom distribution and save it to
   // file
   for (unsigned e = 0; e < n_ele_custom_distribution; e++)
    {
     fprintf(file_pt, "ELEMENT(%i):PROCESSOR(%i)\n", e, input_distribution[e]);
    } // for (i < n_ele_custom_distribution)
   
  }
 else
   {
    // Error, the file could not be opened
    std::ostringstream error_stream;
    error_stream << "Output file could not be opened: " 
                 << file_name << std::endl;
    throw OomphLibError(error_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
 
 // Close the output file
 fclose(file_pt);
 
 oomph_info << "Save custom distribution to file: " << file_name 
            << " [DONE]" << std::endl;
 
}

//========================================================================
/// Demonstrate how to solve Poisson problem using parallel
//unstructured mesh adaptation
//========================================================================
int main(int argc, char* argv[])
{ 

 // initialise MPI
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif
 
 // Get the processor rank (we use this to name the trace file)
 const unsigned my_rank = MPI_Helpers::communicator_pt()->my_rank();
 
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Domain configurations: 1 Square domain
 //                        2 Half circle domain
 //                        3 Half circle domain with internal boundaries
 //                        4 Complex domain with holes
 // Domain configuration
 CommandLineArgs::specify_command_line_flag("--domain_configuration", 
                                            &TestArguments::Domain_configuration);

 // Element size
 CommandLineArgs::specify_command_line_flag("--element_size",
                                            &TestArguments::Element_size);

 // Sample point container
 CommandLineArgs::specify_command_line_flag("--non_ref_bin");
 CommandLineArgs::specify_command_line_flag("--ref_bin");
#ifdef OOMPH_HAS_CGAL
 CommandLineArgs::specify_command_line_flag("--cgal");
#endif

 // Max adaptations
 CommandLineArgs::specify_command_line_flag("--max_adapt",
                                            &TestArguments::Max_adapt);
 
 // Max. permitted error
 CommandLineArgs::specify_command_line_flag("--max_permitted_error", 
                                            &TestArguments::Max_permitted_error);
 
 // Min. permitted error
 CommandLineArgs::specify_command_line_flag("--min_permitted_error", 
                                            &TestArguments::Min_permitted_error);
 
 // Max. element size
 CommandLineArgs::specify_command_line_flag("--max_element_size", 
                                            &TestArguments::Max_element_size);
 
 // Min. element size
 CommandLineArgs::specify_command_line_flag("--min_element_size", 
                                            &TestArguments::Min_element_size);
 
 // Do load balance?
 CommandLineArgs::specify_command_line_flag("--load_balance", 
                                            &TestArguments::Load_balance);
 
 // Max adaptations
 CommandLineArgs::specify_command_line_flag("--folder_distribution_file",
                                            &TestArguments::Folder_distribution_file);
 
 // Parse command line
 CommandLineArgs::parse_and_assign();
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
  
 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
  
 // Do not monitor memory usage
 MemoryUsage::Bypass_all_memory_usage_monitoring=true;
 
 // Swith timings on
 Global_timings::Doc_comprehensive_timings = true;
 
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


 // Use quadratic elements to solve the problem
 //----------------------------------------------
 
 //Set up the problem
 PoissonProblem<ProjectablePoissonElement<TPoissonElement<2,3> > >
   problem(&TanhSolnForPoisson::get_source, TestArguments::Element_size, TestArguments::Domain_configuration);
 
 // Open the trace file (the name of the trace file is different for
 // each processor)
 char trace_filename[100];
 sprintf(trace_filename,"%s/trace_proc%i.dat", 
         doc_info.directory().c_str(), my_rank);
 
 /// Trace file to document error and norm of solution
 ofstream trace_file;
 trace_file.open(trace_filename);
 
 // Doc the initial mesh
 problem.doc_solution(doc_info, trace_file);
 
 //Increment counter for solutions
 doc_info.number()++;
 
 // Once the mesh has been created, suppress the generation of points
 // along the boundaries by Triangle during mesh adaptation
 problem.mesh_pt()->disable_automatic_creation_of_vertices_on_boundaries();
 
 // Set the number of bins for area transfer
 problem.mesh_pt()->nbin_x_for_area_transfer() = 100; // default values
 problem.mesh_pt()->nbin_y_for_area_transfer() = 100; // default values
 
#ifdef OOMPH_HAS_MPI
 // Store the distribution of the elements
 Vector<unsigned> distributed_elements;
 // Should we read the distribution from file
 if(CommandLineArgs::command_line_flag_has_been_set("--folder_distribution_file"))
  {
   // Read a custom distribution from file
   Vector<unsigned> input_distribution;
   problem.read_custom_distribution_from_file(input_distribution);
   // Distribute the problem and store the distribution of the elements
   distributed_elements = problem.distribute(input_distribution);
  }
 else
  {
   // Distribute problem
   distributed_elements = problem.distribute();
   // Uncomment the next lines if you want to save the distribution to
   // file. Create a file with the custom distribution
   //problem.save_custom_distribution_to_file(distributed_elements);
  }
 
 // Output the initial distributed mesh
 char file_initial_distributed_mesh[100];
 // This show us the elements assigned to this processor (remember to
 // disable the output of halo elements)
 problem.mesh_pt()->disable_output_of_halo_elements();
 sprintf(file_initial_distributed_mesh, 
         "%s/output_initial_distributed_mesh_%i.dat", 
         doc_info.directory().c_str(), my_rank);
 
 problem.mesh_pt()->output(file_initial_distributed_mesh, 2);
#endif // OOMPH_HAS_MPI
  
 // Solve the problem doing (parallel unstructured) mesh adaptation
 problem.newton_solve(TestArguments::Max_adapt);
  
 // Document the solution
 problem.doc_solution(doc_info, trace_file);
 
 //Increment counter for solutions
 doc_info.number()++;
 
 // Do load balance?
 if (TestArguments::Load_balance > 0)
  {
   // Perform load balance
   problem.load_balance();
  } // if (load_balance > 0)
 
 //Output solution (balanced solution)
 problem.doc_solution(doc_info, trace_file);
 
 //Increment counter for solutions
 doc_info.number()++;
 
 // Close the trace file
 trace_file.close();
 
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif
 
 // Return exit successful
 return 0;
 
}
