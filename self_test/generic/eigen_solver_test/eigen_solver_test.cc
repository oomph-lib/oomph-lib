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

/// IdentityEigenElement generates an eigenproblem whose Jacobian and mass
/// matrices are equal to the identiy matrix
class IdentityEigenElement : public BaseEigenElement
{
public:
  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    for (unsigned i = 0; i < N_value; i++)
    {
      unsigned local_eqn = internal_local_eqn(Data_index, i);
      for (unsigned j = 0; j < N_value; j++)
      {
        unsigned local_unknown = internal_local_eqn(Data_index, j);
        if (i == j)
        {
          jacobian(local_eqn, local_unknown) += 1;
          mass_matrix(local_eqn, local_unknown) += 1;
        }
      }
    }
  }
};

/// AsymmetricEigenElement generates an eigenproblem whose Jacobian and mass
/// matrices are simple asymmetric examples
class AsymmetricEigenElement : public BaseEigenElement
{
public:
  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    for (unsigned i = 0; i < N_value; i++)
    {
      unsigned local_eqn = internal_local_eqn(Data_index, i);
      for (unsigned j = 0; j < N_value; j++)
      {
        unsigned local_unknown = internal_local_eqn(Data_index, j);
        if (i >= j)
        {
          if (i == j)
          {
            jacobian(local_eqn, local_unknown) += i;
          }
          else
          {
            jacobian(local_eqn, local_unknown) += 1;
          }
          mass_matrix(local_eqn, local_unknown) += 1;
        }
      }
    }
  }
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
    srand(seed);
    for (unsigned i = 0; i < N_value; i++)
    {
      unsigned local_eqn = internal_local_eqn(Data_index, i);
      for (unsigned j = 0; j < N_value; j++)
      {
        unsigned local_unknown = internal_local_eqn(Data_index, j);
        jacobian(local_eqn, local_unknown) += rand() % 256 - 128;
        mass_matrix(local_eqn, local_unknown) += rand() % 256 - 128;
      }
    }
  }
};

/// RosserSymmetricEigenElement generates the classic Rosser eigenproblem
/// matrices. Size 8x8.
class RosserSymmetricEigenElement : public BaseEigenElement
{
public:
  void set_size(const unsigned& n)
  {
    /// Override the input argument as the Rosser matrix is fixed at 8x8
    N_value = 8;

    Data_index = add_internal_data(new Data(N_value));
  }

  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    int A[8][8] = {{611, 196, -192, 407, -8, -52, -49, 29},
                   {196, 899, 113, -192, -71, -43, -8, -44},
                   {-192, 113, 899, 196, 61, 49, 8, 52},
                   {407, -192, 196, 611, 8, 44, 59, -23},
                   {-8, -71, 61, 8, 411, -599, 208, 208},
                   {-52, -43, 49, 44, -599, 411, 208, 208},
                   {-49, -8, 8, 59, 208, 208, 99, -911},
                   {29, -44, 52, -23, 208, 208, -911, 99}};


    for (unsigned i = 0; i < N_value; i++)
    {
      unsigned local_eqn = internal_local_eqn(Data_index, i);
      for (unsigned j = 0; j < N_value; j++)
      {
        unsigned local_unknown = internal_local_eqn(Data_index, j);
        jacobian(local_eqn, local_unknown) += double(A[i][j]);
        if (i == j)
        {
          mass_matrix(local_eqn, local_unknown) += 1;
        }
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

/// TestProblem class. Apply a series of tests upon a given eigensolver using
/// the eigenelement provided as a template.
template<class ELEMENT>
class TestProblem
{
public:
  TestProblem(EigenSolver* const& eigen_solver_pt,
              const unsigned& N,
              const unsigned& n_timing_loops)
    : Eigen_solver_pt(eigen_solver_pt),
      Matrix_size(N),
      N_timing_loops(n_timing_loops)
  {
    // Create eigenproblem
    Problem_pt = new Eigenproblem<ELEMENT>(Matrix_size);

    test_solve_eigenproblem_legacy();

    test_solve_eigenproblem();

    test_solve_adjoint_eigenproblem_legacy();

    test_solve_adjoint_eigenproblem();
  }

  ~TestProblem()
  {
    delete Problem_pt;
  }

  void test_solve_eigenproblem_legacy()
  {
    // Set up additional arguments
    // Output all eigenvalues
    unsigned n_eval = Matrix_size;

    // Store outputs
    Vector<complex<double>> eval(Matrix_size);
    Vector<DoubleVector> evec(Matrix_size);

    // Start clock
    clock_t t_start = clock();
    for (unsigned i = 0; i < N_timing_loops; i++)
    {
      try
      {
        // Call solve_eigenproblem_legacy
        Eigen_solver_pt->solve_eigenproblem_legacy(
          Problem_pt, n_eval, eval, evec);
      }
      catch (...)
      {
        oomph_info << "Failed to run test_solve_eigenproblem_legacy " << endl;
        return;
      }
    }
    // Stop clock
    clock_t t_end = clock();

    // Document duration
    double t_length = (double)(t_end - t_start) / CLOCKS_PER_SEC;
    ofstream timing_stream;
    timing_stream.open("timing.dat", ios_base::app);
    timing_stream << "time: " << t_length / double(N_timing_loops) << endl;
    timing_stream.close();

    // Document solution
    doc_solution(eval);
  }

  void test_solve_eigenproblem()
  {
    // Set up additional arguments
    // Output all eigenvalues
    unsigned n_eval = Matrix_size;

    // Store outputs
    Vector<complex<double>> eval(Matrix_size);
    Vector<DoubleVector> eigenvector_real(Matrix_size);
    Vector<DoubleVector> eigenvector_imag(Matrix_size);

    // Start clock
    clock_t t_start = clock();
    for (unsigned i = 0; i < N_timing_loops; i++)
    {
      try
      {
        // Call solve_eigenproblem_legacy
        Eigen_solver_pt->solve_eigenproblem(
          Problem_pt, n_eval, eval, eigenvector_real, eigenvector_imag);
      }
      catch (...)
      {
        oomph_info << "Failed to run test_solve_eigenproblem " << endl;
        return;
      }
    }
    // Stop clock
    clock_t t_end = clock();

    // Document duration
    double t_length = (double)(t_end - t_start) / CLOCKS_PER_SEC;
    ofstream timing_stream;
    timing_stream.open("timing.dat", ios_base::app);
    timing_stream << "time: " << t_length / double(N_timing_loops) << endl;
    timing_stream.close();

    // Document solution
    doc_solution(eval);
  }

  void test_solve_adjoint_eigenproblem_legacy()
  {
    // Set up additional arguments
    // Output all eigenvalues
    unsigned n_eval = Matrix_size;

    // Store outputs
    Vector<complex<double>> eval(Matrix_size);
    Vector<DoubleVector> evec(Matrix_size);

    // Start clock
    clock_t t_start = clock();
    for (unsigned i = 0; i < N_timing_loops; i++)
    {
      try
      {
        // Call solve_eigenproblem_legacy
        Eigen_solver_pt->solve_adjoint_eigenproblem_legacy(
          Problem_pt, n_eval, eval, evec);
      }
      catch (...)
      {
        oomph_info << "Failed to run test_solve_adjoint_eigenproblem_legacy "
                   << endl;
        return;
      }
    }
    // Stop clock
    clock_t t_end = clock();

    // Document duration
    double t_length = (double)(t_end - t_start) / CLOCKS_PER_SEC;
    ofstream timing_stream;
    timing_stream.open("timing.dat", ios_base::app);
    timing_stream << "time: " << t_length / double(N_timing_loops) << endl;
    timing_stream.close();

    // Document solution
    doc_solution(eval);
  }


  void test_solve_adjoint_eigenproblem()
  {
    // Set up additional arguments
    // Output all eigenvalues
    unsigned n_eval = Matrix_size;

    // Store outputs
    Vector<complex<double>> eval(Matrix_size);
    Vector<DoubleVector> eigenvector_real(Matrix_size);
    Vector<DoubleVector> eigenvector_imag(Matrix_size);

    // Start clock
    clock_t t_start = clock();
    for (unsigned i = 0; i < N_timing_loops; i++)
    {
      try
      {
        // Call solve_eigenproblem_legacy
        Eigen_solver_pt->solve_adjoint_eigenproblem(
          Problem_pt, n_eval, eval, eigenvector_real, eigenvector_imag);
      }
      catch (...)
      {
        oomph_info << "Failed to run test_solve_adjoint_eigenproblem " << endl;
        return;
      }
    }
    // Stop clock
    clock_t t_end = clock();

    // Document duration
    double t_length = (double)(t_end - t_start) / CLOCKS_PER_SEC;
    ofstream timing_stream;
    timing_stream.open("timing.dat", ios_base::app);
    timing_stream << "time: " << t_length / double(N_timing_loops) << endl;
    timing_stream.close();

    // Document solution
    doc_solution(eval);
  }

  void doc_solution(Vector<complex<double>> eval)
  {
    for (unsigned i = 0; i < Matrix_size; i++)
    {
      oomph_info << eval[i] << ", ";
    }
    oomph_info << endl;
  }

private:
  EigenSolver* Eigen_solver_pt;
  unsigned Matrix_size;
  unsigned N_timing_loops;
  Problem* Problem_pt;
};

/// TestSolver class. Tests the templated eigensolver against a series of
/// eigenelements.
template<class SOLVER>
class TestSolver
{
public:
  TestSolver(const unsigned N, const unsigned n_timing_loops)
    : Matrix_size(N), N_timing_loops(n_timing_loops)
  {
    EigenSolver* eigen_solver_pt = new SOLVER;

    TestProblem<IdentityEigenElement>(
      eigen_solver_pt, Matrix_size, N_timing_loops);

    TestProblem<AsymmetricEigenElement>(
      eigen_solver_pt, Matrix_size, N_timing_loops);

    TestProblem<RosserSymmetricEigenElement>(
      eigen_solver_pt, 8, N_timing_loops);

    TestProblem<RandomAsymmetricEigenElement>(
      eigen_solver_pt, Matrix_size, N_timing_loops);

    delete eigen_solver_pt;
  }

private:
  unsigned Matrix_size;
  unsigned N_timing_loops;
};

/// Main function. Apply solver tests to each eigensolver.
int main()
{
  /// Temporary (minimum working example)
  // mwe();

  // Number of times to repeat the operation for better timings
  const unsigned n_timing_loops = 2;

  // Matrix dimensions
  const unsigned N = 32;

  TestSolver<LAPACK_QZ>(N, n_timing_loops);

  TestSolver<ANASAZI>(N, n_timing_loops);

  return (EXIT_SUCCESS);
}
