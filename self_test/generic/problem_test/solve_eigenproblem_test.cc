// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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

// Oomph-lib includes
#include "generic.h"
#include <random>

using namespace std;
using namespace oomph;

/// Base eigenproblem element class used to generate the Jacobian and mass
/// matrix which correspond to the eigenproblem (J - lambda M) v = 0
class BaseEigenElement : public GeneralisedElement
{
public:
  BaseEigenElement() {}

  void set_size(const unsigned& n)
  {
    N_value = n;

    Data_index = add_internal_data(new Data(N_value));
  }

  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix) = 0;

protected:
  unsigned N_value;
  unsigned Data_index;
};

/// RandomAsymmetricEigenElement generates an eigenproblem whose Jacobian and
/// mass matrices have random integer elements between -128 and 127. Seed = 0.
class RandomAsymmetricEigenElement : public BaseEigenElement
{
public:
  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    unsigned seed = 0;
    mt19937 generator(seed);
    for (unsigned i = 0; i < N_value; i++)
    {
      unsigned local_eqn = internal_local_eqn(Data_index, i);
      for (unsigned j = 0; j < N_value; j++)
      {
        unsigned local_unknown = internal_local_eqn(Data_index, j);
        jacobian(local_eqn, local_unknown) += generator() % 256 - 128;
        mass_matrix(local_eqn, local_unknown) += generator() % 256 - 128;
      }
    }
  }
};

/// Eigenproblem class. Creates a mesh with a single eigenelement and assigns
/// the appropriate equation numbers.
template<class ELEMENT>
class Eigenproblem : public Problem
{
public:
  Eigenproblem(const unsigned& size)
  {
    this->mesh_pt() = new Mesh;

    ELEMENT* el_pt = new ELEMENT;

    el_pt->set_size(size);

    this->mesh_pt()->add_element_pt(el_pt);

    assign_eqn_numbers();
  }

  ~Eigenproblem()
  {
    delete this->mesh_pt();
  }
};

void doc_solution(string filename,
                  Vector<std::complex<double>>& eigenvalue,
                  Vector<DoubleVector>& eigenvector_real,
                  Vector<DoubleVector>& eigenvector_imag)
{
  ofstream output_stream;
  output_stream.open(filename);
  unsigned N = eigenvalue.size();
  for (unsigned i = 0; i < N; i++)
  {
    output_stream << eigenvalue[i].real() << " , " << eigenvalue[i].imag()
                  << endl;
  }
  output_stream << endl;

  for (unsigned i = 0; i < N; i++)
  {
    for (unsigned j = 0; j < N; j++)
    {
      output_stream << eigenvector_real[j][i] << " , " << eigenvector_imag[j][i]
                    << endl;
    }
    output_stream << endl;
  }
  output_stream.close();
}

void doc_solution_legacy(string filename,
                         Vector<std::complex<double>>& eigenvalue,
                         Vector<DoubleVector>& eigenvector)
{
  ofstream output_stream;
  output_stream.open(filename);
  unsigned N = eigenvalue.size();
  for (unsigned i = 0; i < N; i++)
  {
    output_stream << eigenvalue[i].real() << " , " << eigenvalue[i].imag()
                  << endl;
  }
  output_stream << endl;

  for (unsigned i = 0; i < N; i++)
  {
    for (unsigned j = 0; j < N; j++)
    {
      output_stream << eigenvector[j][i] << endl;
    }
    output_stream << endl;
  }
  output_stream.close();
}

/// Main function. Apply solver tests to each eigensolver.
int main()
{
  DocInfo doc_info;
  doc_info.set_directory("RESLT/");

  // Matrix dimensions
  const unsigned N = 32;

  // Create eigenproblem
  Eigenproblem<RandomAsymmetricEigenElement> problem(N);

  // problem.eigen_solver_pt() = new ANASAZI;
  // Anasazi::Use_temporary_code_for_andrew_legacy_version = true;

  // Set up additional arguments
  // Output all eigenvalues
  unsigned n_eval = N;

  // Store outputs
  Vector<complex<double>> eigenvalue(N);
  Vector<DoubleVector> eigenvector_real(N);
  Vector<DoubleVector> eigenvector_imag(N);

  /// Test solve_eigenproblem
  problem.solve_eigenproblem(
    n_eval, eigenvalue, eigenvector_real, eigenvector_imag);
  doc_solution(doc_info.directory() + "solve_eigenproblem_test.dat",
               eigenvalue,
               eigenvector_real,
               eigenvector_imag);

  /// ----------------------------------------
  ///              Legacy tests
  /// ----------------------------------------

  Vector<DoubleVector> eigenvector(N);

  /// Test solve_eigenproblem_legacy
  problem.solve_eigenproblem_legacy(n_eval, eigenvalue, eigenvector);
  doc_solution_legacy(doc_info.directory() +
                        "solve_eigenproblem_legacy_test.dat",
                      eigenvalue,
                      eigenvector);


  return (EXIT_SUCCESS);
}
