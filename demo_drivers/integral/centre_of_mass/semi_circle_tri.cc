//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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

#include <iostream>

#include "generic.h"
#include "meshes.h"
#include "integral.h"

using namespace std;
using namespace oomph;

void x_moment_fct(const double t, const Vector<double>& x, double& f)
{
  f = x[0];
}

void y_moment_fct(const double t, const Vector<double>& x, double& f)
{
  f = x[1];
}

namespace problem_parameter
{
  enum
  {
    LOWER_BOUNDARY,
    CURVED_BOUNDARY,
  };
}; // namespace problem_parameter

template<class ELEMENT>
class IntegralProblem : public Problem
{
public:
  /// Constructor
  IntegralProblem();

  /// Destructor
  ~IntegralProblem(){};

  /// Document the solution
  void doc_solution(DocInfo& doc_info);

private:
  /// Generate mesh
  void generate_mesh();

  /// Upcast elements and finalise setup
  void upcast_and_finalise_elements();

  /// Pointer to the outer boundary polyline
  TriangleMeshPolygon* Rect_boundary_polyline_pt;

  /// Pointer to the "bulk" integral mesh
  Mesh* Integral_mesh_pt;

  /// Pointer to the info mesh
  Mesh* Info_mesh_pt;
};

template<class ELEMENT>
IntegralProblem<ELEMENT>::IntegralProblem()
{
  cout << "Generate mesh" << endl;
  generate_mesh();

  cout << "Upcast and finalise elements" << endl;
  upcast_and_finalise_elements();

  // Setup equation numbering scheme
  cout << "Number of equations: " << endl;
  cout << assign_eqn_numbers() << endl;
  cout << "Number of unknowns: " << endl;
  cout << ndof() << endl;
}

template<class ELEMENT>
void IntegralProblem<ELEMENT>::generate_mesh()
{
  double l_x = 2.0;

  cout << "Generate boundary" << endl;
  Vector<TriangleMeshCurveSection*> boundary_polyline_pt;

  Vector<Vector<double>> vertex_coord(2);
  for (unsigned i = 0; i < 2; i++)
  {
    vertex_coord[i].resize(2);
  }

  /// First vertex
  vertex_coord[0][0] = 0;
  vertex_coord[0][1] = 0;
  /// Second vertex
  vertex_coord[1][0] = l_x;
  vertex_coord[1][1] = 0;

  boundary_polyline_pt.push_back(
    new TriangleMeshPolyLine(vertex_coord, problem_parameter::LOWER_BOUNDARY));

  const unsigned npoints = 32;
  const double zeta_step = MathematicalConstants::Pi / double(npoints - 1);
  // Intrinsic coordinate along GeomObject defining the bubble
  Vector<double> zeta(1);
  const double x_c = 0.5 * l_x;
  const double y_c = 0.0;
  const double r = 0.5 * l_x;
  Circle circle(x_c, y_c, r);

  /// Remove first vertex
  vertex_coord.erase(vertex_coord.begin());
  /// Fill in the semi circle
  for (unsigned ipoint = 1; ipoint < npoints; ipoint++)
  {
    // Get the coordinates
    zeta[0] = zeta_step * double(ipoint);
    Vector<double> current_vertex(2);
    circle.position(zeta, current_vertex);
    vertex_coord.push_back(current_vertex);
  }

  boundary_polyline_pt.push_back(
    new TriangleMeshPolyLine(vertex_coord, problem_parameter::CURVED_BOUNDARY));

  // Create the triangle mesh polygon for rectangle boundary
  Rect_boundary_polyline_pt = new TriangleMeshPolygon(boundary_polyline_pt);

  // Convert to "closed curve" objects
  TriangleMeshClosedCurve* rect_closed_curve_pt = Rect_boundary_polyline_pt;

  // Generate mesh parameters for external mesh generator "Triangle"
  TriangleMeshParameters triangle_mesh_parameters(rect_closed_curve_pt);

  const double maximum_default_element_area = 1e-2;
  triangle_mesh_parameters.element_area() = maximum_default_element_area;

  // Call external mesh generator
  Integral_mesh_pt = new TriangleMesh<ELEMENT>(triangle_mesh_parameters);

  this->Info_mesh_pt = new Mesh;
  this->Info_mesh_pt->add_element_pt(new InfoElement(3));

  add_sub_mesh(Integral_mesh_pt);
  add_sub_mesh(Info_mesh_pt);

  build_global_mesh();
}

template<class ELEMENT>
void IntegralProblem<ELEMENT>::upcast_and_finalise_elements()
{
  /// Bulk mesh
  // Find number of elements in mesh
  unsigned n_element = Integral_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Integral_mesh_pt->element_pt(i));

    el_pt->setup_integrand(
      this->Info_mesh_pt->element_pt(0)->internal_data_pt(0));
    el_pt->setup_integrand(
      this->Info_mesh_pt->element_pt(0)->internal_data_pt(1), &x_moment_fct);
    el_pt->setup_integrand(
      this->Info_mesh_pt->element_pt(0)->internal_data_pt(2), &y_moment_fct);
  }
}

template<class ELEMENT>
void IntegralProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  string data_directory = doc_info.directory();
  string filename = data_directory + "integral.dat";
  ofstream output_stream;
  output_stream.open(filename.c_str());
  const double area =
    this->Info_mesh_pt->element_pt(0)->internal_data_pt(0)->value(0);
  const double x_moment =
    this->Info_mesh_pt->element_pt(0)->internal_data_pt(1)->value(0);
  const double y_moment =
    this->Info_mesh_pt->element_pt(0)->internal_data_pt(2)->value(0);
  output_stream << "Area: " << area << endl;
  output_stream << "X moment: " << x_moment << endl;
  output_stream << "Y moment: " << y_moment << endl;
  output_stream << "Com: " << x_moment / area << ", " << y_moment / area
                << endl;
  output_stream.close();

  filename = data_directory + "domain.dat";
  output_stream.open(filename.c_str());
  Integral_mesh_pt->output(output_stream);
  output_stream.close();
}

int main(int argc, char* argv[])
{
  /// Store command line arguments
  CommandLineArgs::setup(argc, argv);

  /// Create a DocInfo object for output processing
  DocInfo doc_info;
  doc_info.set_directory("RESLT/");

  IntegralProblem<TIntegralElement<2,3>> problem;

  /// Call problem solve
  problem.newton_solve();

  /// Document solution
  problem.doc_solution(doc_info);

  return 0;
}
