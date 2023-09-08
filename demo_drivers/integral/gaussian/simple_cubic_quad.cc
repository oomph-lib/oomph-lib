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
#include <iostream>

#include "generic.h"
#include "meshes.h"
#include "integral.h"

using namespace std;
using namespace oomph;

// Multivariate Gaussian function
void gaussian_fct(const double t, const Vector<double>& x, double& f)
{
  const double a = 1;
  const double b = 1;
  const double c = 1;

  const unsigned N = x.size();
  double sum = 0;
  for (unsigned n = 0; n < N; n++)
  {
    sum += pow(x[n] - b, 2.0);
  }
  f = a * exp(-sum / pow(c, 2.0));
}

// Specific Problem class
template<class ELEMENT>
class IntegralProblem : public Problem
{
public:
  // Constructor
  IntegralProblem();

  // Destructor
  ~IntegralProblem(){};

  // Document the solution
  void doc_solution(DocInfo& doc_info);

private:
  // Generate mesh
  void generate_mesh();

  // Upcast elements and finalise setup
  void upcast_and_finalise_elements();

  // Pointer to the "bulk" integral mesh
  Mesh* Integral_mesh_pt;

  // Pointer to the info mesh
  Mesh* Info_mesh_pt;
};

// Constructor
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
}

template<class ELEMENT>
void IntegralProblem<ELEMENT>::generate_mesh()
{
  // Number of elements per size
  const unsigned nx = 10;
  const unsigned ny = 10;
  const unsigned nz = 5;
  // Length of each size
  const double lx = 2.0;
  const double ly = 1.0;
  const double lz = 0.5;

  // Create the domain mesh
  this->Integral_mesh_pt = new SimpleCubicMesh<ELEMENT>(nx, ny, nz, lx, ly, lz);

  // Create the info element (we only need one internal data) and its mesh
  this->Info_mesh_pt = new Mesh;
  const unsigned n_data = 1;
  this->Info_mesh_pt->add_element_pt(new InfoElement(n_data));

  // Add sub meshes to the global mesh
  this->add_sub_mesh(this->Integral_mesh_pt);
  this->add_sub_mesh(this->Info_mesh_pt);

  // Build the global mesh
  this->build_global_mesh();
}

// Upcast the elements and setup the integrand function pointers
template<class ELEMENT>
void IntegralProblem<ELEMENT>::upcast_and_finalise_elements()
{
  // Bulk mesh
  // Find number of elements in mesh
  unsigned n_element = this->Integral_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    ELEMENT* el_pt =
      dynamic_cast<ELEMENT*>(this->Integral_mesh_pt->element_pt(i));

    // Setup the integrand function and where the integral is to be stored.
    el_pt->setup_integrand(
      this->Info_mesh_pt->element_pt(0)->internal_data_pt(0), &gaussian_fct);
  }
}

// Document the solution and integral
template<class ELEMENT>
void IntegralProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  string data_directory = doc_info.directory();

  // Output the integral
  string filename = data_directory + "integral.dat";
  ofstream output_stream;
  output_stream.open(filename.c_str());
  double my_integral =
    this->Info_mesh_pt->element_pt(0)->internal_data_pt(0)->value(0);
  output_stream << "Integral: " << my_integral << endl;
  output_stream.close();

  // Output the domain
  filename = data_directory + "domain.dat";
  output_stream.open(filename.c_str());
  Integral_mesh_pt->output(output_stream);
  output_stream.close();
}

int main(int argc, char* argv[])
{
  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  // Create a DocInfo object for output processing
  DocInfo doc_info;
  doc_info.set_directory("RESLT/");

  // Create the integral problem
  IntegralProblem<QIntegralElement<3, 3>> problem;

  // Call problem solve
  problem.newton_solve();

  // Document solution
  problem.doc_solution(doc_info);

  return 0;
}
