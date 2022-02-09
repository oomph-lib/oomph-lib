// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
// Driver to solve the harmonic equation with homogeneous Dirichlet boundary
// conditions.

// Generic oomph-lib routines
#include "generic.h"

// Include the mesh
#include "meshes/one_d_mesh.h"

using namespace std;

using namespace oomph;

/// Base eigenproblem element class used to generate the Jacobian and mass
/// matrix which correspond to the eigenproblem (J - lambda M) v = 0
class BaseEigenElement : public FiniteElement
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

  void shape(const Vector<double>& s, Shape& psi) const
  {
    cout << "Shape" << endl;
  };

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


//=================================================================
/// A class for all elements that solve the simple one-dimensional
/// eigenvalue problem
/// \f[
/// \frac{\partial^2 u}{\partial x_i^2}  + \lambda u = 0
/// \f]
/// These elements are very closely related to the Poisson
/// elements and could inherit from them. They are here developed
/// from scratch for pedagogical purposes.
/// This class  contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//================================================================
class HarmonicEquations : public virtual FiniteElement
{
public:
  /// Empty Constructor
  HarmonicEquations() {}

  /// Access function: Eigenfunction value at local node n
  /// Note that solving the eigenproblem does not assign values
  /// to this storage space. It is used for output purposes only.
  virtual inline double u(const unsigned& n) const
  {
    return nodal_value(n, 0);
  }

  /// Output the eigenfunction with default number of plot points
  void output(ostream& outfile)
  {
    unsigned nplot = 5;
    output(outfile, nplot);
  }

  /// Output FE representation of soln: x,y,u or x,y,z,u at
  /// Nplot  plot points
  void output(ostream& outfile, const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(1);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);
      // Output the coordinate and the eigenfunction
      outfile << interpolated_x(s, 0) << " " << interpolated_u(s) << std::endl;
    }
    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }

  /// Assemble the contributions to the jacobian and mass
  /// matrices
  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Set up memory for the shape functions and their derivatives
    Shape psi(n_node);
    DShape dpsidx(n_node, 1);

    // Set the number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Integers to store the local equation and unknown numbers
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_eulerian_at_knot(ipt, psi, dpsidx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Assemble the contributions to the mass matrix
      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the local equation number
        local_eqn = u_local_eqn(l);
        int global_eqn = eqn_number(local_eqn);
        /*IF it's not a boundary condition*/
        if (local_eqn >= 0)
        {
          // Loop over the shape functions
          for (unsigned l2 = 0; l2 < n_node; l2++)
          {
            local_unknown = u_local_eqn(l2);
            int global_unknown = eqn_number(local_unknown);
            // If at a non-zero degree of freedom add in the entry
            if (local_unknown >= 0)
            {
              // if (global_unknown == global_eqn)
              // {
              //   jacobian(global_eqn, global_unknown) += global_eqn;
              // }
              // else
              // {
              //   jacobian(global_eqn, global_unknown) += 1;
              // }
              // mass_matrix(global_eqn, global_unknown) += 1;

              // if (local_unknown == local_eqn)
              // {
              //   jacobian(local_eqn, local_unknown) += local_eqn;
              // }
              // else
              // {
              //   jacobian(local_eqn, local_unknown) += 1;
              // }
              // mass_matrix(local_eqn, local_unknown) += 1;

              jacobian(local_eqn, local_unknown) +=
                dpsidx(l, 0) * dpsidx(l2, 0) * W;
              mass_matrix(local_eqn, local_unknown) += psi(l) * psi(l2) * W;
            }
          }
        }
      }
    }
  } // end_of_fill_in_contribution_to_jacobian_and_mass_matrix

  /// Return FE representation of function value u(s) at local coordinate s
  inline double interpolated_u(const Vector<double>& s) const
  {
    unsigned n_node = nnode();

    // Local shape functions
    Shape psi(n_node);

    // Find values of basis function
    this->shape(s, psi);

    // Initialise value of u
    double interpolated_u = 0.0;

    // Loop over the local nodes and sum
    for (unsigned l = 0; l < n_node; l++)
    {
      interpolated_u += u(l) * psi[l];
    }

    // Return the interpolated value of the eigenfunction
    return (interpolated_u);
  }

protected:
  /// Shape/test functions and derivs w.r.t. to global coords at
  /// local coord. s; return  Jacobian of mapping
  virtual double dshape_eulerian(const Vector<double>& s,
                                 Shape& psi,
                                 DShape& dpsidx) const = 0;

  /// Shape/test functions and derivs w.r.t. to global coords at
  /// integration point ipt; return  Jacobian of mapping
  virtual double dshape_eulerian_at_knot(const unsigned& ipt,
                                         Shape& psi,
                                         DShape& dpsidx) const = 0;

  /// Access function that returns the local equation number
  /// of the unknown in the problem. Default is to assume that it is the
  /// first (only) value stored at the node.
  virtual inline int u_local_eqn(const unsigned& n)
  {
    return nodal_local_eqn(n, 0);
  }

private:
};


//======================================================================
/// QHarmonicElement<NNODE_1D> elements are 1D  Elements with
/// NNODE_1D nodal points that are used to solve the Harmonic eigenvalue
/// Problem described by HarmonicEquations.
//======================================================================
template<unsigned NNODE_1D>
class QHarmonicElement : public virtual QElement<1, NNODE_1D>,
                         public HarmonicEquations
{
public:
  /// Constructor: Call constructors for QElement and
  /// Poisson equations
  QHarmonicElement() : QElement<1, NNODE_1D>(), HarmonicEquations() {}

  ///  Required  # of `values' (pinned or dofs)
  /// at node n
  inline unsigned required_nvalue(const unsigned& n) const
  {
    return 1;
  }

  /// Output function overloaded from HarmonicEquations
  void output(ostream& outfile)
  {
    HarmonicEquations::output(outfile);
  }

  ///  Output function overloaded from HarmonicEquations
  void output(ostream& outfile, const unsigned& Nplot)
  {
    HarmonicEquations::output(outfile, Nplot);
  }


protected:
  /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
  inline double dshape_eulerian(const Vector<double>& s,
                                Shape& psi,
                                DShape& dpsidx) const
  {
    return QElement<1, NNODE_1D>::dshape_eulerian(s, psi, dpsidx);
  }


  /// Shape, test functions & derivs. w.r.t. to global coords. at
  /// integration point ipt. Return Jacobian.
  inline double dshape_eulerian_at_knot(const unsigned& ipt,
                                        Shape& psi,
                                        DShape& dpsidx) const
  {
    return QElement<1, NNODE_1D>::dshape_eulerian_at_knot(ipt, psi, dpsidx);
  }

}; // end_of_QHarmonic_class_definition


//==start_of_problem_class============================================
/// 1D Harmonic problem in unit interval.
//====================================================================
template<class ELEMENT, class EIGEN_SOLVER>
class HarmonicProblem : public Problem
{
public:
  /// Constructor: Pass number of elements and pointer to source function
  HarmonicProblem(const unsigned& n_element);

  /// Destructor (empty)
  ~HarmonicProblem()
  {
    delete this->mesh_pt();
    delete this->eigen_solver_pt();
  }

  /// Solve the problem
  void solve(const unsigned& label);

  /// Compare given Jacobian to a numerical one
  void debug_jacobian()
  {
    CRDoubleMatrix jacobian;
    CRDoubleMatrix mass_matrix;

    cout << "Get Jacobian" << endl;
    get_eigenproblem_matrices(mass_matrix, jacobian);

    cout << ndof() << jacobian.nrow() << jacobian.ncol() << endl;
    for (unsigned i = 0; i < ndof(); i++)
    {
      for (unsigned j = 0; j < ndof(); j++)
      {
        printf("i: %4u, j: %4u", i, j);
        printf(", %10.5g", jacobian(j, i));
        cout << endl;
      }
    }
    cout << "End of Jacobian" << endl;
  }
}; // end of problem class


//=====start_of_constructor===============================================
/// Constructor for 1D Harmonic problem in unit interval.
/// Discretise the 1D domain with n_element elements of type ELEMENT.
/// Specify function pointer to source function.
//========================================================================
template<class ELEMENT, class EIGEN_SOLVER>
HarmonicProblem<ELEMENT, EIGEN_SOLVER>::HarmonicProblem(
  const unsigned& n_element)
{
  // Create the eigen solver
  this->eigen_solver_pt() = new EIGEN_SOLVER;

#ifdef OOMPH_HAS_TRILINOS
  // hierher Temporary work-around to keep the legacy version working
  Anasazi::Use_temporary_code_for_andrew_legacy_version = true;
#endif

  // //Get the positive eigenvalues, shift is zero by default
  // static_cast<EIGEN_SOLVER*>(eigen_solver_pt())
  //  ->get_eigenvalues_right_of_shift();

  {
    // Problem::mesh_pt() = new Mesh;
    // AsymmetricEigenElement* el_pt = new AsymmetricEigenElement;
    // el_pt->set_size(n_element);
    // mesh_pt()->add_element_pt(el_pt);
  }

  {
    // Set domain length
    double L = 1.0;
    Problem::mesh_pt() = new OneDMesh<ELEMENT>(n_element, L);
    mesh_pt()->boundary_node_pt(0, 0)->pin(0);
    mesh_pt()->boundary_node_pt(1, 0)->pin(0);
  }

  // Setup equation numbering scheme
  cout << "Number of DoF: " << assign_eqn_numbers() << endl;
} // end of constructor

//=======================start_of_solve==============================
/// Solve the eigenproblem
//===================================================================
template<class ELEMENT, class EIGEN_SOLVER>
void HarmonicProblem<ELEMENT, EIGEN_SOLVER>::solve(const unsigned& label)
{
  // Set external storage for the eigenvalues
  Vector<complex<double>> eigenvalues;
  // Set external storage for the eigenvectors
  Vector<DoubleVector> eigenvectors;
  // Desired number eigenvalues
  unsigned n_eval = 4;

  // Solve the eigenproblem
  this->solve_eigenproblem_legacy(n_eval, eigenvalues, eigenvectors);

  char filename[100];
  sprintf(filename, "eigenvalues%i.dat", label);

  // Open an output file for the sorted eigenvalues
  ofstream evalues(filename);
  for (unsigned i = 0; i < n_eval; i++)
  {
    // Print to screen
    cout << eigenvalues[i].real() << " " << eigenvalues[i].imag() << std::endl;
    // Send to file
    evalues << eigenvalues[i].real() << " " << eigenvalues[i].imag()
            << std::endl;
  }

  evalues.close();
} // end_of_solve


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//======start_of_main==================================================
/// Driver for 1D Poisson problem
//=====================================================================
int main(int argc, char** argv)
{
// Want to test Trilinos if we have it, so we must initialise MPI
// if we have compiled with it
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::init(argc, argv);
#endif

  // Set up the problem:
  unsigned n_element = 8; // Number of elements

#ifdef OOMPH_HAS_TRILINOS
  // Solve with Anasazi
  HarmonicProblem<QHarmonicElement<3>, ANASAZI> problem(n_element);
  //problem.debug_jacobian();
  problem.solve(3);
#endif

#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

} // end of main
