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

  /// Pointer to the "bulk" integral mesh
  Mesh* Integral_mesh_pt;

  /// Pointer to the info mesh
  Mesh* Info_mesh_pt;
};

template<class ELEMENT>
IntegralProblem<ELEMENT>::IntegralProblem()
{
  cout << "Generate mesh" << endl;
  this->generate_mesh();

  cout << "Upcast and finalise elements" << endl;
  this->upcast_and_finalise_elements();

  // Setup equation numbering scheme
  cout << "Number of equations: " << endl;
  cout << this->assign_eqn_numbers() << endl;
  cout << "Number of unknowns: " << endl;
  cout << this->ndof() << endl;
}

template<class ELEMENT>
void IntegralProblem<ELEMENT>::generate_mesh()
{
  const unsigned nx = 10;
  const unsigned ny = 10;
  const unsigned nx_cut_out = 5;
  const unsigned ny_cut_out = 5;
  const double lx = 2.0;
  const double ly = 1.0;
  this->Integral_mesh_pt =
    new BackwardStepQuadMesh<ELEMENT>(nx, ny, nx_cut_out, ny_cut_out, lx, ly);

  this->Info_mesh_pt = new Mesh;
  this->Info_mesh_pt->add_element_pt(new InfoElement(3));

  this->add_sub_mesh(this->Integral_mesh_pt);
  this->add_sub_mesh(this->Info_mesh_pt);

  this->build_global_mesh();
}

template<class ELEMENT>
void IntegralProblem<ELEMENT>::upcast_and_finalise_elements()
{
  /// Bulk mesh
  // Find number of elements in mesh
  unsigned n_element = this->Integral_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    ELEMENT* el_pt =
      dynamic_cast<ELEMENT*>(this->Integral_mesh_pt->element_pt(i));

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

  IntegralProblem<QIntegralElement<2,3>> problem;

  /// Call problem solve
  problem.newton_solve();

  /// Document solution
  problem.doc_solution(doc_info);

  return 0;
}
