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

  void fill_in_contribution_to_residuals(oomph::Vector<double>& residuals) {}

protected:
  unsigned N_value;
  unsigned Data_index;
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
    cout << "N_value: " << N_value << endl;
    for (unsigned i = 0; i < N_value; i++)
    {
      unsigned local_eqn = internal_local_eqn(Data_index, i);
      for (unsigned j = 0; j < N_value; j++)
      {
        unsigned local_unknown = internal_local_eqn(Data_index, j);
        cout << i << ", " << local_eqn << ", " << j << ", " << local_unknown
             << endl;

        if (i == j)
        {
          jacobian(local_eqn, local_unknown) += 1;
          mass_matrix(local_eqn, local_unknown) += 1;
        }

        // if (i >= j)
        //{
        //  if (i == j)
        //  {
        //    jacobian(local_eqn, local_unknown) += i;
        //  }
        //  else
        //  {
        //    jacobian(local_eqn, local_unknown) += 1;
        //  }
        //  mass_matrix(local_eqn, local_unknown) += 1;
        //}
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

    //   el_pt = new ELEMENT;

    //   el_pt->set_size(size);

    //   this->mesh_pt()->add_element_pt(el_pt);

    //   el_pt = new ELEMENT;

    //   el_pt->set_size(size);

    //   this->mesh_pt()->add_element_pt(el_pt);

    // build_global_mesh();

    cout << "Number of equations:" << assign_eqn_numbers() << endl;
  }

  ~Eigenproblem()
  {
    delete this->mesh_pt();
  }
};

/// Minimum working example to demonstrate failing
/// ANASAZI.solve_eigenproblem_legacy
void mwe()
{
  const unsigned N = 64;
  const unsigned n_eval = 1;
  const bool do_adjoint_problem = false;

  EigenSolver* Eigen_solver_pt = 0;
#ifdef OOMPH_HAS_TRILINOS
  Eigen_solver_pt = new ANASAZI;
#else
  return;
#endif
  Problem* Problem_pt = new Eigenproblem<AsymmetricEigenElement>(N);

  Vector<complex<double>> eval(N);
  Vector<DoubleVector> evec(N);
  Eigen_solver_pt->solve_eigenproblem_legacy(
    Problem_pt, n_eval, eval, evec, do_adjoint_problem);
  // Vector<DoubleVector> evecI(N);
  // Eigen_solver_pt->solve_eigenproblem(
  //  Problem_pt, n_eval, eval, evec, evecI, do_adjoint_problem);

  for (unsigned i = 0; i < N; i++)
  {
    oomph_info << eval[i] << ", ";
  }
  oomph_info << endl;
}


/// Main function.
int main(int argc, char** argv)
{
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::init(argc, argv);
#endif

#ifdef OOMPH_HAS_TRILINOS
  Anasazi::Use_temporary_code_for_andrew_legacy_version = true;
  mwe();
#endif

#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

  return (EXIT_SUCCESS);
}
