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
// OOOMPH-LIB includes
#include "assembly_handler.h"
#include "elements.h"
#include "problem.h"
#include "mesh.h"

namespace oomph
{
  ///////////////////////////////////////////////////////////////////////
  // Non-inline functions for the AssemblyHandler class
  //////////////////////////////////////////////////////////////////////

  //===================================================================
  /// Get the number of elemental degrees of freedom. Direct call
  /// to the function in the element.
  //===================================================================
  unsigned AssemblyHandler::ndof(GeneralisedElement* const& elem_pt)
  {
    return elem_pt->ndof();
  }

  //===================================================================
  /// Get the vector of dofs in the element elem_pt at time level t
  /// Direct call to the function in the element.
  //===================================================================
  void AssemblyHandler::dof_vector(GeneralisedElement* const& elem_pt,
                                   const unsigned& t,
                                   Vector<double>& dof)
  {
    return elem_pt->dof_vector(t, dof);
  }

  //===================================================================
  /// Get the vector of pointers to dofs in the element elem_pt
  /// Direct call to the function in the element.
  //===================================================================
  void AssemblyHandler::dof_pt_vector(GeneralisedElement* const& elem_pt,
                                      Vector<double*>& dof_pt)
  {
    return elem_pt->dof_pt_vector(dof_pt);
  }

  /// Return the t-th level of storage associated with the i-th
  /// (local) dof stored in the problem
  double& AssemblyHandler::local_problem_dof(Problem* const& problem_pt,
                                             const unsigned& t,
                                             const unsigned& i)
  {
    return *(problem_pt->dof_pt(i) + t);
  }


  //==================================================================
  /// Get the global equation number of the local unknown. Direct call
  /// to the function in the element.
  //==================================================================
  unsigned long AssemblyHandler::eqn_number(GeneralisedElement* const& elem_pt,
                                            const unsigned& ieqn_local)
  {
    return elem_pt->eqn_number(ieqn_local);
  }

  //==================================================================
  /// Get the residuals by calling the underlying element's residuals
  /// directly.
  //=================================================================
  void AssemblyHandler::get_residuals(GeneralisedElement* const& elem_pt,
                                      Vector<double>& residuals)
  {
    elem_pt->get_residuals(residuals);
  }

  //=======================================================================
  /// Calculate the elemental Jacobian matrix "d equation / d variable" by
  /// calling the element's get_jacobian function.
  //======================================================================
  void AssemblyHandler::get_jacobian(GeneralisedElement* const& elem_pt,
                                     Vector<double>& residuals,
                                     DenseMatrix<double>& jacobian)
  {
    elem_pt->get_jacobian(residuals, jacobian);
  }

  //=======================================================================
  /// Calculate all desired vectors and matrices that are required by
  /// the element by calling the get_jacobian function
  //=======================================================================
  void AssemblyHandler::get_all_vectors_and_matrices(
    GeneralisedElement* const& elem_pt,
    Vector<Vector<double>>& vec,
    Vector<DenseMatrix<double>>& matrix)
  {
    get_jacobian(elem_pt, vec[0], matrix[0]);
  }

  //=======================================================================
  /// Calculate the derivative of the residuals with respect to
  /// a parameter, by calling the elemental function
  //======================================================================
  void AssemblyHandler::get_dresiduals_dparameter(
    GeneralisedElement* const& elem_pt,
    double* const& parameter_pt,
    Vector<double>& dres_dparam)
  {
    elem_pt->get_dresiduals_dparameter(parameter_pt, dres_dparam);
  }


  //=====================================================================
  /// Calculate the derivative of the residuals and jacobian
  /// with respect to a parameter by calling the elemental function
  //========================================================================
  void AssemblyHandler::get_djacobian_dparameter(
    GeneralisedElement* const& elem_pt,
    double* const& parameter_pt,
    Vector<double>& dres_dparam,
    DenseMatrix<double>& djac_dparam)
  {
    elem_pt->get_djacobian_dparameter(parameter_pt, dres_dparam, djac_dparam);
  }

  /// Calculate the product of the Hessian (derivative of Jacobian with
  /// respect to all variables) an eigenvector, Y, and
  /// other specified vectors, C
  /// (d(J_{ij})/d u_{k}) Y_{j} C_{k}
  /// At the moment the dof pointer is passed in to enable
  /// easy calculation of finite difference default
  void AssemblyHandler::get_hessian_vector_products(
    GeneralisedElement* const& elem_pt,
    Vector<double> const& Y,
    DenseMatrix<double> const& C,
    DenseMatrix<double>& product)
  {
    elem_pt->get_hessian_vector_products(Y, C, product);
  }


  //=======================================================================
  /// Return the eigenfunction(s) associated with the bifurcation that
  /// has been detected in bifurcation tracking problems. Default
  /// Broken implementation
  //=======================================================================
  double* AssemblyHandler::bifurcation_parameter_pt() const
  {
    std::ostringstream error_stream;
    error_stream << "There is no bifurcation parameter associated with the "
                    "current assembly handler.\n"
                 << "Eigenfunction are only calculated by the Fold, PitchFork "
                    "and Hopf Handlers"
                 << "\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    return 0;
  }


  //=======================================================================
  /// Return the eigenfunction(s) associated with the bifurcation that
  /// has been detected in bifurcation tracking problems. Default
  /// Broken implementation
  //=======================================================================
  void AssemblyHandler::get_eigenfunction(Vector<DoubleVector>& eigenfunction)
  {
    std::ostringstream error_stream;
    error_stream << "There is no eigenfunction associated with the current "
                    "assembly handler.\n"
                 << "Eigenfunction are only calculated by the Fold, PitchFork "
                    "and Hopf Handlers"
                 << "\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  ///==========================================================================
  /// Compute the inner products of the given vector of pairs of
  /// history values over the element. The values of the index in the pair
  /// may be the same.
  //===========================================================================
  void AssemblyHandler::get_inner_products(
    GeneralisedElement* const& elem_pt,
    Vector<std::pair<unsigned, unsigned>> const& history_index,
    Vector<double>& inner_product)
  {
    elem_pt->get_inner_products(history_index, inner_product);
  }

  //==========================================================================
  /// Compute the vectors that when taken as a dot product with
  /// other history values give the inner product over the element.
  /// In other words if we call get_inner_product_vectors({0,1},output)
  /// output[0] will be a vector such that dofs.output[0] gives the inner
  /// product of current dofs with themselves.
  //==========================================================================
  void AssemblyHandler::get_inner_product_vectors(
    GeneralisedElement* const& elem_pt,
    Vector<unsigned> const& history_index,
    Vector<Vector<double>>& inner_product_vector)
  {
    elem_pt->get_inner_product_vectors(history_index, inner_product_vector);
  }


  ///////////////////////////////////////////////////////////////////////
  // Non-inline functions for the ExplicitTimeStepHandler class
  //////////////////////////////////////////////////////////////////////


  //===================================================================
  /// Get the number of elemental degrees of freedom. Direct call
  /// to the function in the element.
  //===================================================================
  unsigned ExplicitTimeStepHandler::ndof(GeneralisedElement* const& elem_pt)
  {
    return elem_pt->ndof();
  }

  //==================================================================
  /// Get the global equation number of the local unknown. Direct call
  /// to the function in the element.
  //==================================================================
  unsigned long ExplicitTimeStepHandler::eqn_number(
    GeneralisedElement* const& elem_pt, const unsigned& ieqn_local)
  {
    return elem_pt->eqn_number(ieqn_local);
  }

  //==================================================================
  /// Call the element's residuals
  //=================================================================
  void ExplicitTimeStepHandler::get_residuals(
    GeneralisedElement* const& elem_pt, Vector<double>& residuals)
  {
    elem_pt->get_residuals(residuals);
  }

  //=======================================================================
  /// Replace get jacobian with the call to get the mass matrix
  //======================================================================
  void ExplicitTimeStepHandler::get_jacobian(GeneralisedElement* const& elem_pt,
                                             Vector<double>& residuals,
                                             DenseMatrix<double>& jacobian)
  {
    elem_pt->get_mass_matrix(residuals, jacobian);
  }


  //=======================================================================
  /// Calculate all desired vectors and matrices that are required by
  /// the problem  by calling those of the underlying element.
  //=======================================================================
  void ExplicitTimeStepHandler::get_all_vectors_and_matrices(
    GeneralisedElement* const& elem_pt,
    Vector<Vector<double>>& vec,
    Vector<DenseMatrix<double>>& matrix)
  {
#ifdef PARANOID
    // Check dimension
    if (matrix.size() != 1)
    {
      throw OomphLibError("ExplicitTimeSteps should return one matrix",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif
    // Get the residuals and <mass matrix
    elem_pt->get_mass_matrix(vec[0], matrix[0]);
  }


  ///////////////////////////////////////////////////////////////////////
  // Non-inline functions for the EigenProblemHandler class
  //////////////////////////////////////////////////////////////////////


  //===================================================================
  /// Get the number of elemental degrees of freedom. Direct call
  /// to the function in the element.
  //===================================================================
  unsigned EigenProblemHandler::ndof(GeneralisedElement* const& elem_pt)
  {
    return elem_pt->ndof();
  }

  //==================================================================
  /// Get the global equation number of the local unknown. Direct call
  /// to the function in the element.
  //==================================================================
  unsigned long EigenProblemHandler::eqn_number(
    GeneralisedElement* const& elem_pt, const unsigned& ieqn_local)
  {
    return elem_pt->eqn_number(ieqn_local);
  }

  //==================================================================
  /// Cannot call get_residuals for an eigenproblem, so throw an error
  //=================================================================
  void EigenProblemHandler::get_residuals(GeneralisedElement* const& elem_pt,
                                          Vector<double>& residuals)
  {
    throw OomphLibError(
      "An eigenproblem does not have a get_residuals function",
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
  }

  //=======================================================================
  /// Cannot call get_jacobian for an eigenproblem, so throw an error
  //======================================================================
  void EigenProblemHandler::get_jacobian(GeneralisedElement* const& elem_pt,
                                         Vector<double>& residuals,
                                         DenseMatrix<double>& jacobian)
  {
    throw OomphLibError("An eigenproblem does not have a get_jacobian function",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }


  //=======================================================================
  /// Calculate all desired vectors and matrices that are required by
  /// the problem  by calling those of the underlying element.
  //=======================================================================
  void EigenProblemHandler::get_all_vectors_and_matrices(
    GeneralisedElement* const& elem_pt,
    Vector<Vector<double>>& vec,
    Vector<DenseMatrix<double>>& matrix)
  {
#ifdef PARANOID
    // Check dimension
    if (matrix.size() != 2)
    {
      throw OomphLibError("EigenProblems should return two matrices",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif
    // Find the number of variables
    unsigned n_var = elem_pt->ndof();
    // Assign a dummy residuals vector
    Vector<double> dummy(n_var);
    // Get the jacobian and mass matrices
    elem_pt->get_jacobian_and_mass_matrix(dummy, matrix[0], matrix[1]);

    // If we have a non-zero shift, then shift the A matrix
    if (Sigma_real != 0.0)
    {
      // Set the shifted matrix
      for (unsigned i = 0; i < n_var; i++)
      {
        for (unsigned j = 0; j < n_var; j++)
        {
          matrix[0](i, j) -= Sigma_real * matrix[1](i, j);
        }
      }
    }
  }


  //======================================================================
  /// Clean up the memory that may have been allocated by the solver
  //=====================================================================
  AugmentedBlockFoldLinearSolver::~AugmentedBlockFoldLinearSolver()
  {
    if (Alpha_pt != 0)
    {
      delete Alpha_pt;
    }
    if (E_pt != 0)
    {
      delete E_pt;
    }
  }

  //===================================================================
  /// Use a block factorisation to solve the augmented system
  /// associated with a PitchFork bifurcation.
  //===================================================================
  void AugmentedBlockFoldLinearSolver::solve(Problem* const& problem_pt,
                                             DoubleVector& result)
  {
    // if the result is setup then it should not be distributed
#ifdef PARANOID
    if (result.built())
    {
      if (result.distributed())
      {
        throw OomphLibError("The result vector must not be distributed",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Locally cache the pointer to the handler.
    FoldHandler* handler_pt =
      static_cast<FoldHandler*>(problem_pt->assembly_handler_pt());

    // Switch the handler to "block solver" mode
    handler_pt->solve_augmented_block_system();

    // We need to find out the number of dofs in the problem
    unsigned n_dof = problem_pt->ndof();

    // create the linear algebra distribution for this solver
    // currently only global (non-distributed) distributions are allowed
    LinearAlgebraDistribution dist(problem_pt->communicator_pt(), n_dof, false);
    this->build_distribution(dist);

    // if the result vector is not setup then rebuild with distribution = global
    if (!result.built())
    {
      result.build(this->distribution_pt(), 0.0);
    }

    // Setup storage for temporary vectors
    DoubleVector a(this->distribution_pt(), 0.0),
      b(this->distribution_pt(), 0.0);

    // Allocate storage for Alpha which can be used in the resolve
    if (Alpha_pt != 0)
    {
      delete Alpha_pt;
    }
    Alpha_pt = new DoubleVector(this->distribution_pt(), 0.0);

    // We are going to do resolves using the underlying linear solver
    Linear_solver_pt->enable_resolve();

    // Solve the first system Aa = R
    Linear_solver_pt->solve(problem_pt, a);

    // The vector in the top right-hand block is the jacobian multiplied
    // by the null vector

    // Get the present null vector from the handler
    DoubleVector y(this->distribution_pt(), 0.0);
    for (unsigned n = 0; n < (n_dof - 1); ++n)
    {
      y[n] = handler_pt->Y[n];
    }
    // For simplicity later, add a zero at the end
    y[n_dof - 1] = 0.0;

    // Loop over the elements and assemble the product
    // Local storage for the product terms
    DoubleVector Jy(this->distribution_pt(), 0.0);

    // Calculate the product of the jacobian matrices, etc
    unsigned long n_element = problem_pt->mesh_pt()->nelement();
    for (unsigned long e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the ndofs in each element
      unsigned n_var = elem_pt->ndof();
      // Get the jacobian matrix
      DenseMatrix<double> jac(n_var);
      // storage for the residual
      Vector<double> res(n_var);
      // Get unperturbed raw jacobian
      elem_pt->get_jacobian(res, jac);

      // Multiply the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = elem_pt->eqn_number(n);
        for (unsigned m = 0; m < n_var; m++)
        {
          unsigned unknown = elem_pt->eqn_number(m);
          Jy[eqn_number] += jac(n, m) * y[unknown];
        }
      }
    }
    // The final entry of the vector will be zero

    // Now resolve to find alpha
    Linear_solver_pt->resolve(Jy, *Alpha_pt);

    // The vector that multiplie the product matrix is actually y - alpha
    DoubleVector y_minus_alpha(this->distribution_pt(), 0.0);
    for (unsigned n = 0; n < n_dof; n++)
    {
      y_minus_alpha[n] = y[n] - (*Alpha_pt)[n];
    }

    // We can now construct our multipliers
    // Prepare to scale
    double dof_length = 0.0, a_length = 0.0, alpha_length = 0.0;
    for (unsigned n = 0; n < n_dof; n++)
    {
      if (std::fabs(problem_pt->dof(n)) > dof_length)
      {
        dof_length = std::fabs(problem_pt->dof(n));
      }
      if (std::fabs(a[n]) > a_length)
      {
        a_length = std::fabs(a[n]);
      }
      if (std::fabs(y_minus_alpha[n]) > alpha_length)
      {
        alpha_length = std::fabs(y_minus_alpha[n]);
      }
    }

    double a_mult = dof_length / a_length;
    double alpha_mult = dof_length / alpha_length;
    const double FD_step = 1.0e-8;
    a_mult += FD_step;
    alpha_mult += FD_step;
    a_mult *= FD_step;
    alpha_mult *= FD_step;

    // Local storage for the product terms
    DoubleVector Jprod_a(this->distribution_pt(), 0.0),
      Jprod_alpha(this->distribution_pt(), 0.0);

    // Calculate the product of the jacobian matrices, etc
    for (unsigned long e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the ndofs in each element
      unsigned n_var = handler_pt->ndof(elem_pt);
      // Get the jacobian matrices
      DenseMatrix<double> jac(n_var), jac_a(n_var), jac_alpha(n_var);
      // elemental residual storage
      Vector<double> res_elemental(n_var);
      // Get unperturbed jacobian
      handler_pt->get_jacobian(elem_pt, res_elemental, jac);

      // Backup the dofs
      Vector<double> dof_bac(n_var);
      // Perturb the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        dof_bac[n] = problem_pt->dof(eqn_number);
        // Pertub by vector a
        problem_pt->dof(eqn_number) += a_mult * a[eqn_number];
      }

      problem_pt->actions_after_change_in_bifurcation_parameter();

      // Now get the new jacobian
      handler_pt->get_jacobian(elem_pt, res_elemental, jac_a);

      // Perturb the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        problem_pt->dof(eqn_number) = dof_bac[n];
        // Pertub by vector a
        problem_pt->dof(eqn_number) += alpha_mult * y_minus_alpha[eqn_number];
      }

      problem_pt->actions_after_change_in_bifurcation_parameter();

      // Now get the new jacobian
      handler_pt->get_jacobian(elem_pt, res_elemental, jac_alpha);

      // Reset the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        problem_pt->dof(eqn_number) = dof_bac[n];
      }

      problem_pt->actions_after_change_in_bifurcation_parameter();

      // OK, now work out the products
      // Note the (n_var-1), we are only interested in the non-augmented
      // jacobian
      for (unsigned n = 0; n < (n_var - 1); n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        double prod_a = 0.0, prod_alpha = 0.0;
        for (unsigned m = 0; m < (n_var - 1); m++)
        {
          unsigned unknown = handler_pt->eqn_number(elem_pt, m);
          prod_a += (jac_a(n, m) - jac(n, m)) * y[unknown];
          prod_alpha += (jac_alpha(n, m) - jac(n, m)) * y[unknown];
        }
        Jprod_a[eqn_number] += prod_a / a_mult;
        Jprod_alpha[eqn_number] += prod_alpha / alpha_mult;
      }
    }

    Jprod_alpha[n_dof - 1] = 0.0;
    Jprod_a[n_dof - 1] = 0.0;

    // OK, now we can formulate the next vectors
    // The assumption here is that the result has been set to the
    // residuals.
    for (unsigned n = 0; n < n_dof - 1; n++)
    {
      b[n] = result[n_dof + n] - Jprod_a[n];
    }
    // The final residual is the entry corresponding to the original parameter
    b[n_dof - 1] = result[n_dof - 1];

    // Allocate storage for E which can be used in the resolve
    if (E_pt != 0)
    {
      delete E_pt;
    }
    E_pt = new DoubleVector(this->distribution_pt(), 0.0);
    DoubleVector f(this->distribution_pt(), 0.0);
    Linear_solver_pt->resolve(b, f);
    Linear_solver_pt->resolve(Jprod_alpha, *E_pt);

    // Calculate the final entry in the vector e
    const double e_final = (*E_pt)[n_dof - 1];
    // Calculate the final entry in the vector d
    const double d_final = f[n_dof - 1] / e_final;
    // Assemble the final corrections
    for (unsigned n = 0; n < n_dof - 1; n++)
    {
      result[n] = a[n] - (*Alpha_pt)[n] * d_final + d_final * y[n];
      result[n_dof + n] = f[n] - (*E_pt)[n] * d_final;
    }
    // The result corresponding to the parameter
    result[n_dof - 1] = a[n_dof - 1] - (*Alpha_pt)[n_dof - 1] * d_final;

    // The sign of the jacobian is the sign of the final entry in e
    problem_pt->sign_of_jacobian() =
      static_cast<int>(std::fabs(e_final) / e_final);

    // Switch things to our block solver
    handler_pt->solve_full_system();

    // If we are not storing things, clear the results
    if (!Enable_resolve)
    {
      // We no longer need to store the matrix
      Linear_solver_pt->disable_resolve();
      delete Alpha_pt;
      Alpha_pt = 0;
      delete E_pt;
      E_pt = 0;
    }
    // Otherwise, also store the problem pointer
    else
    {
      Problem_pt = problem_pt;
    }
  }


  //======================================================================
  // Hack the re-solve to use the block factorisation
  //======================================================================
  void AugmentedBlockFoldLinearSolver::resolve(const DoubleVector& rhs,
                                               DoubleVector& result)
  {
    // Check that the factors have been stored
    if (Alpha_pt == 0)
    {
      throw OomphLibError("The required vectors have not been stored",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Get the pointer to the problem
    Problem* const problem_pt = Problem_pt;

    FoldHandler* handler_pt =
      static_cast<FoldHandler*>(problem_pt->assembly_handler_pt());

    // Switch things to our block solver
    handler_pt->solve_augmented_block_system();

    // We need to find out the number of dofs in the problem
    unsigned n_dof = problem_pt->ndof();

#ifdef PARANOID
    // if the result is setup then it should not be distributed
    if (result.built())
    {
      if (result.distributed())
      {
        throw OomphLibError("The result vector must not be distributed",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
    // the rhs must be setup
    if (!rhs.built())
    {
      throw OomphLibError("The rhs vector must be setup",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // if the result vector is not setup then rebuild with distribution = global
    if (!result.built())
    {
      result.build(rhs.distribution_pt(), 0.0);
    }

    // Setup storage
    DoubleVector a(this->distribution_pt(), 0.0),
      b(this->distribution_pt(), 0.0);

    // Set the values of the a vector
    for (unsigned n = 0; n < (n_dof - 1); n++)
    {
      a[n] = rhs[n];
    }
    // The entry associated with the additional parameter is zero
    a[n_dof - 1] = 0.0;

    Linear_solver_pt->enable_resolve();

    // Copy rhs vector into local storage so it doesn't get overwritten
    // if the linear solver decides to initialise the solution vector, say,
    // which it's quite entitled to do!
    DoubleVector input_a(a);

    Linear_solver_pt->resolve(input_a, a);

    // We can now construct our multipliers
    // Prepare to scale
    double dof_length = 0.0, a_length = 0.0;
    for (unsigned n = 0; n < n_dof; n++)
    {
      if (std::fabs(problem_pt->dof(n)) > dof_length)
      {
        dof_length = std::fabs(problem_pt->dof(n));
      }

      if (std::fabs(a[n]) > a_length)
      {
        a_length = std::fabs(a[n]);
      }
    }
    double a_mult = dof_length / a_length;
    const double FD_step = 1.0e-8;
    a_mult += FD_step;
    a_mult *= FD_step;

    DoubleVector Jprod_a(this->distribution_pt(), 0.0);

    unsigned long n_element = problem_pt->mesh_pt()->nelement();
    for (unsigned long e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the ndofs in each element
      unsigned n_var = handler_pt->ndof(elem_pt);
      // Get some jacobian matrices
      DenseMatrix<double> jac(n_var), jac_a(n_var);
      // the elemental residual
      Vector<double> res_elemental(n_var);
      // Get unperturbed jacobian
      handler_pt->get_jacobian(elem_pt, res_elemental, jac);

      // Backup the dofs
      Vector<double> dof_bac(n_var);
      // Perturb the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        dof_bac[n] = problem_pt->dof(eqn_number);
        // Pertub by vector a
        problem_pt->dof(eqn_number) += a_mult * a[eqn_number];
      }

      problem_pt->actions_after_change_in_bifurcation_parameter();

      // Now get the new jacobian
      handler_pt->get_jacobian(elem_pt, res_elemental, jac_a);

      // Reset the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        problem_pt->dof(eqn_number) = dof_bac[n];
      }

      problem_pt->actions_after_change_in_bifurcation_parameter();

      // OK, now work out the products
      for (unsigned n = 0; n < (n_var - 1); n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        double prod_a = 0.0;
        for (unsigned m = 0; m < (n_var - 1); m++)
        {
          unsigned unknown = handler_pt->eqn_number(elem_pt, m);
          prod_a += (jac_a(n, m) - jac(n, m)) * handler_pt->Y[unknown];
        }
        Jprod_a[eqn_number] += prod_a / a_mult;
      }
    }

    Jprod_a[n_dof - 1] = 0.0;

    // OK, now we can formulate the next vectors
    for (unsigned n = 0; n < n_dof - 1; n++)
    {
      b[n] = rhs[n_dof + n] - Jprod_a[n];
    }
    // The residuals for the final entry should be those associated
    // with the parameter
    b[n_dof - 1] = rhs[n_dof - 1];

    DoubleVector f(this->distribution_pt(), 0.0);

    Linear_solver_pt->resolve(b, f);

    // Calculate the final entry in the vector d
    const double d_final = f[n_dof - 1] / (*E_pt)[n_dof - 1];
    // Assemble the final corrections
    for (unsigned n = 0; n < n_dof - 1; n++)
    {
      result[n] = a[n] - (*Alpha_pt)[n] * d_final + d_final * handler_pt->Y[n];
      result[n_dof + n] = f[n] - (*E_pt)[n] * d_final;
    }
    // The result corresponding to the paramater
    result[n_dof - 1] = a[n_dof - 1] - (*Alpha_pt)[n_dof - 1] * d_final;

    Linear_solver_pt->disable_resolve();

    // Switch things to our block solver
    handler_pt->solve_full_system();
  }


  ///////////////////////////////////////////////////////////////////////
  // Non-inline functions for the FoldHandler class
  //////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Constructor: Initialise the fold handler,
  /// by setting initial guesses for Y, Phi and calculating Count.
  /// If the system changes, a new  handler must be constructed.
  //========================================================================
  FoldHandler::FoldHandler(Problem* const& problem_pt,
                           double* const& parameter_pt)
    : Solve_which_system(Full_augmented), Parameter_pt(parameter_pt)
  {
    // Set the problem pointer
    Problem_pt = problem_pt;
    // Set the number of degrees of freedom
    Ndof = problem_pt->ndof();

    // create the linear algebra distribution for this solver
    // currently only global (non-distributed) distributions are allowed
    LinearAlgebraDistribution* dist_pt =
      new LinearAlgebraDistribution(problem_pt->communicator_pt(), Ndof, false);

    // Resize the vectors of additional dofs and constants
    Phi.resize(Ndof);
    Y.resize(Ndof);
    Count.resize(Ndof, 0);

    // Loop over all the elements in the problem
    unsigned n_element = problem_pt->mesh_pt()->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the local freedoms in the element
      unsigned n_var = elem_pt->ndof();
      for (unsigned n = 0; n < n_var; n++)
      {
        // Increase the associated global equation number counter
        ++Count[elem_pt->eqn_number(n)];
      }
    }

    // Calculate the value Phi by
    // solving the system JPhi = dF/dlambda

    // Locally cache the linear solver
    LinearSolver* const linear_solver_pt = problem_pt->linear_solver_pt();

    // Save the status before entry to this routine
    bool enable_resolve = linear_solver_pt->is_resolve_enabled();

    // We need to do a resolve
    linear_solver_pt->enable_resolve();

    // Storage for the solution
    DoubleVector x(dist_pt, 0.0);

    // Solve the standard problem, we only want to make sure that
    // we factorise the matrix, if it has not been factorised. We shall
    // ignore the return value of x.
    linear_solver_pt->solve(problem_pt, x);

    // Get the vector dresiduals/dparameter
    problem_pt->get_derivative_wrt_global_parameter(parameter_pt, x);

    // Copy rhs vector into local storage so it doesn't get overwritten
    // if the linear solver decides to initialise the solution vector, say,
    // which it's quite entitled to do!
    DoubleVector input_x(x);

    // Now resolve the system with the new RHS and overwrite the solution
    linear_solver_pt->resolve(input_x, x);

    // Restore the storage status of the linear solver
    if (enable_resolve)
    {
      linear_solver_pt->enable_resolve();
    }
    else
    {
      linear_solver_pt->disable_resolve();
    }

    // Add the global parameter as an unknown to the problem
    problem_pt->Dof_pt.push_back(parameter_pt);


    // Normalise the initial guesses for phi
    double length = 0.0;
    for (unsigned n = 0; n < Ndof; n++)
    {
      length += x[n] * x[n];
    }
    length = sqrt(length);

    // Now add the null space components to the problem unknowns
    // and initialise them and Phi to the same normalised values
    for (unsigned n = 0; n < Ndof; n++)
    {
      problem_pt->Dof_pt.push_back(&Y[n]);
      Y[n] = Phi[n] = -x[n] / length;
    }

    // delete the dist_pt
    problem_pt->Dof_distribution_pt->build(
      problem_pt->communicator_pt(), Ndof * 2 + 1, true);
    // Remove all previous sparse storage used during Jacobian assembly
    Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

    delete dist_pt;
  }


  /// Constructor in which the eigenvector is passed as an initial
  /// guess
  FoldHandler::FoldHandler(Problem* const& problem_pt,
                           double* const& parameter_pt,
                           const DoubleVector& eigenvector)
    : Solve_which_system(Full_augmented), Parameter_pt(parameter_pt)
  {
    // Set the problem pointer
    Problem_pt = problem_pt;
    // Set the number of degrees of freedom
    Ndof = problem_pt->ndof();

    // create the linear algebra distribution for this solver
    // currently only global (non-distributed) distributions are allowed
    LinearAlgebraDistribution* dist_pt =
      new LinearAlgebraDistribution(problem_pt->communicator_pt(), Ndof, false);

    // Resize the vectors of additional dofs and constants
    Phi.resize(Ndof);
    Y.resize(Ndof);
    Count.resize(Ndof, 0);

    // Loop over all the elements in the problem
    unsigned n_element = problem_pt->mesh_pt()->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the local freedoms in the element
      unsigned n_var = elem_pt->ndof();
      for (unsigned n = 0; n < n_var; n++)
      {
        // Increase the associated global equation number counter
        ++Count[elem_pt->eqn_number(n)];
      }
    }

    // Add the global parameter as an unknown to the problem
    problem_pt->Dof_pt.push_back(parameter_pt);


    // Normalise the initial guesses for the eigenvecor
    double length = 0.0;
    for (unsigned n = 0; n < Ndof; n++)
    {
      length += eigenvector[n] * eigenvector[n];
    }
    length = sqrt(length);

    // Now add the null space components to the problem unknowns
    // and initialise them and Phi to the same normalised values
    for (unsigned n = 0; n < Ndof; n++)
    {
      problem_pt->Dof_pt.push_back(&Y[n]);
      Y[n] = Phi[n] = eigenvector[n] / length;
    }

    // delete the dist_pt
    problem_pt->Dof_distribution_pt->build(
      problem_pt->communicator_pt(), Ndof * 2 + 1, true);
    // Remove all previous sparse storage used during Jacobian assembly
    Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

    delete dist_pt;
  }


  /// Constructor in which the eigenvector and normalisation
  /// vector are  passed as an initial guess
  FoldHandler::FoldHandler(Problem* const& problem_pt,
                           double* const& parameter_pt,
                           const DoubleVector& eigenvector,
                           const DoubleVector& normalisation)
    : Solve_which_system(Full_augmented), Parameter_pt(parameter_pt)
  {
    // Set the problem pointer
    Problem_pt = problem_pt;
    // Set the number of degrees of freedom
    Ndof = problem_pt->ndof();

    // create the linear algebra distribution for this solver
    // currently only global (non-distributed) distributions are allowed
    LinearAlgebraDistribution* dist_pt =
      new LinearAlgebraDistribution(problem_pt->communicator_pt(), Ndof, false);

    // Resize the vectors of additional dofs and constants
    Phi.resize(Ndof);
    Y.resize(Ndof);
    Count.resize(Ndof, 0);

    // Loop over all the elements in the problem
    unsigned n_element = problem_pt->mesh_pt()->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the local freedoms in the element
      unsigned n_var = elem_pt->ndof();
      for (unsigned n = 0; n < n_var; n++)
      {
        // Increase the associated global equation number counter
        ++Count[elem_pt->eqn_number(n)];
      }
    }

    // Add the global parameter as an unknown to the problem
    problem_pt->Dof_pt.push_back(parameter_pt);


    // Normalise the initial guesses for the eigenvecor
    double length = 0.0;
    for (unsigned n = 0; n < Ndof; n++)
    {
      length += eigenvector[n] * normalisation[n];
    }
    length = sqrt(length);

    // Now add the null space components to the problem unknowns
    // and initialise them and Phi to the same normalised values
    for (unsigned n = 0; n < Ndof; n++)
    {
      problem_pt->Dof_pt.push_back(&Y[n]);
      Y[n] = eigenvector[n] / length;
      Phi[n] = normalisation[n];
    }

    // delete the dist_pt
    problem_pt->Dof_distribution_pt->build(
      problem_pt->communicator_pt(), Ndof * 2 + 1, true);
    // Remove all previous sparse storage used during Jacobian assembly
    Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

    delete dist_pt;
  }


  //=================================================================
  /// The number of unknowns is 2n+1
  //================================================================
  unsigned FoldHandler::ndof(GeneralisedElement* const& elem_pt)
  {
    unsigned raw_ndof = elem_pt->ndof();
    // Return different values depending on the type of block decomposition
    switch (Solve_which_system)
    {
      case Full_augmented:
        return (2 * raw_ndof + 1);
        break;

      case Block_augmented_J:
        return (raw_ndof + 1);
        break;

      case Block_J:
        return raw_ndof;
        break;

      default:
        std::ostringstream error_stream;
        error_stream
          << "The Solve_which_system flag can only take values 0, 1, 2"
          << " not " << Solve_which_system << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  //=====================================================================
  /// Return the global equation number associated with local equation
  /// number ieqn_local. We choose to number the unknowns according
  /// to the augmented system.
  //=======================================================================
  unsigned long FoldHandler::eqn_number(GeneralisedElement* const& elem_pt,
                                        const unsigned& ieqn_local)
  {
    // Find the number of non-augmented dofs in the element
    unsigned raw_ndof = elem_pt->ndof();
    // Storage for the global eqn number
    unsigned long global_eqn = 0;
    // If we are a "standard" unknown, just return the standard equation number
    if (ieqn_local < raw_ndof)
    {
      global_eqn = elem_pt->eqn_number(ieqn_local);
    }
    // Otherwise if we are at an unknown corresponding to the bifurcation
    // parameter return the global equation number of the parameter
    else if (ieqn_local == raw_ndof)
    {
      global_eqn = Ndof;
    }
    // Otherwise we are in the unknown corresponding to a null vector
    // return the global unknown Ndof + 1 + global unknown of "standard"
    // unknown.
    else
    {
      global_eqn = Ndof + 1 + elem_pt->eqn_number(ieqn_local - 1 - raw_ndof);
    }

    // Return the global equation number
    return global_eqn;
  }

  //====================================================================
  /// Formulate the augmented system
  //====================================================================
  void FoldHandler::get_residuals(GeneralisedElement* const& elem_pt,
                                  Vector<double>& residuals)
  {
    // Need to get raw residuals and jacobian
    unsigned raw_ndof = elem_pt->ndof();

    // Find out which system we are solving
    switch (Solve_which_system)
    {
        // If we are solving the standard system
      case Block_J:
      {
        // Get the basic residuals
        elem_pt->get_residuals(residuals);
      }
      break;

      // If we are solving the augmented-by-one system
      case Block_augmented_J:
      {
        // Get the basic residuals
        elem_pt->get_residuals(residuals);

        // Zero the final residual
        residuals[raw_ndof] = 0.0;
      }
      break;

      // If we are solving the full augmented system
      case Full_augmented:
      {
        DenseMatrix<double> jacobian(raw_ndof);
        // Get the basic residuals and jacobian initially
        elem_pt->get_jacobian(residuals, jacobian);

        // The normalisation equation must be initialised to -1.0/number of
        // elements
        residuals[raw_ndof] = -1.0 / Problem_pt->mesh_pt()->nelement();

        // Now assemble the equations Jy = 0
        for (unsigned i = 0; i < raw_ndof; i++)
        {
          residuals[raw_ndof + 1 + i] = 0.0;
          for (unsigned j = 0; j < raw_ndof; j++)
          {
            residuals[raw_ndof + 1 + i] +=
              jacobian(i, j) * Y[elem_pt->eqn_number(j)];
          }
          // Add the contribution to phi.y=1
          // Need to divide by the number of elements that contribute to this
          // unknown so that we do get phi.y exactly.
          unsigned global_eqn = elem_pt->eqn_number(i);
          residuals[raw_ndof] +=
            (Phi[global_eqn] * Y[global_eqn]) / Count[global_eqn];
        }
      }
      break;

      default:
        std::ostringstream error_stream;
        error_stream
          << "The Solve_which_system flag can only take values 0, 1, 2"
          << " not " << Solve_which_system << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //=================================================================
  /// Calculate the elemental Jacobian matrix "d equation
  /// / d variable" for the augmented system
  //=================================================================
  void FoldHandler::get_jacobian(GeneralisedElement* const& elem_pt,
                                 Vector<double>& residuals,
                                 DenseMatrix<double>& jacobian)
  {
    // Find the number of augmented dofs
    unsigned augmented_ndof = ndof(elem_pt);
    // Find the non-augmented dofs
    unsigned raw_ndof = elem_pt->ndof();

    // Which system are we solving
    switch (Solve_which_system)
    {
        // If we are solving the original system
      case Block_J:
      {
        // Just get the raw jacobian and residuals
        elem_pt->get_jacobian(residuals, jacobian);
      }
      break;

      // If we are solving the augmented-by-one system
      case Block_augmented_J:
      {
        // Get the full residuals, we need them
        get_residuals(elem_pt, residuals);

        // Need to get the raw jacobian (and raw residuals)
        Vector<double> newres(augmented_ndof);
        elem_pt->get_jacobian(newres, jacobian);

        // Now do finite differencing stuff
        const double FD_step = 1.0e-8;
        // Fill in the first lot of finite differences
        {
          // increase the global parameter
          double* unknown_pt = Problem_pt->Dof_pt[Ndof];
          double init = *unknown_pt;
          *unknown_pt += FD_step;

          // Now do any possible updates
          Problem_pt->actions_after_change_in_bifurcation_parameter();

          // Get the new (modified) residuals
          get_residuals(elem_pt, newres);

          // The final column  is given by the difference
          // between the residuals
          for (unsigned n = 0; n < raw_ndof; n++)
          {
            jacobian(n, augmented_ndof - 1) =
              (newres[n] - residuals[n]) / FD_step;
          }
          // Reset the global parameter
          *unknown_pt = init;

          // Now do any possible updates
          Problem_pt->actions_after_change_in_bifurcation_parameter();
        }

        // Fill in the bottom row
        for (unsigned n = 0; n < raw_ndof; n++)
        {
          unsigned local_eqn = elem_pt->eqn_number(n);
          jacobian(augmented_ndof - 1, n) = Phi[local_eqn] / Count[local_eqn];
        }
      }
      break;

        // Otherwise solving the full system
      case Full_augmented:
      {
        // Get the residuals for the augmented system
        get_residuals(elem_pt, residuals);

        // Need to get the raw residuals and jacobian
        Vector<double> newres(raw_ndof);
        DenseMatrix<double> newjac(raw_ndof);
        elem_pt->get_jacobian(newres, jacobian);

        // Fill in the jacobian on the diagonal sub-block of
        // the null-space equations
        for (unsigned n = 0; n < raw_ndof; n++)
        {
          for (unsigned m = 0; m < raw_ndof; m++)
          {
            jacobian(raw_ndof + 1 + n, raw_ndof + 1 + m) = jacobian(n, m);
          }
        }

        // Now finite difference wrt the global unknown
        const double FD_step = 1.0e-8;
        // Fill in the first lot of finite differences
        {
          // increase the global parameter
          double* unknown_pt = Problem_pt->Dof_pt[Ndof];
          double init = *unknown_pt;
          *unknown_pt += FD_step;
          // Need to update the function
          Problem_pt->actions_after_change_in_bifurcation_parameter();

          // Get the new raw residuals and jacobian
          elem_pt->get_jacobian(newres, newjac);

          // The end of the first row is given by the difference
          // between the residuals
          for (unsigned n = 0; n < raw_ndof; n++)
          {
            jacobian(n, raw_ndof) = (newres[n] - residuals[n]) / FD_step;
            // The end of the second row is given by the difference multiplied
            // by the product
            for (unsigned l = 0; l < raw_ndof; l++)
            {
              jacobian(raw_ndof + 1 + n, raw_ndof) +=
                (newjac(n, l) - jacobian(n, l)) * Y[elem_pt->eqn_number(l)] /
                FD_step;
            }
          }
          // Reset the global parameter
          *unknown_pt = init;

          // Need to update the function
          Problem_pt->actions_after_change_in_bifurcation_parameter();
        }

        // Now fill in the first column of the second rows
        {
          for (unsigned n = 0; n < raw_ndof; n++)
          {
            unsigned long global_eqn = eqn_number(elem_pt, n);
            // Increase the first lot
            double* unknown_pt = Problem_pt->Dof_pt[global_eqn];
            double init = *unknown_pt;
            *unknown_pt += FD_step;
            Problem_pt->actions_before_newton_convergence_check(); /// ALICE

            // Get the new jacobian
            elem_pt->get_jacobian(newres, newjac);

            // Work out the differences
            for (unsigned k = 0; k < raw_ndof; k++)
            {
              // jacobian(raw_ndof+k,n) = 0.0;
              for (unsigned l = 0; l < raw_ndof; l++)
              {
                jacobian(raw_ndof + 1 + k, n) +=
                  (newjac(k, l) - jacobian(k, l)) * Y[elem_pt->eqn_number(l)] /
                  FD_step;
              }
            }
            *unknown_pt = init;
            Problem_pt->actions_before_newton_convergence_check(); /// ALICE
          }
        }

        // Fill in the row corresponding to the parameter
        for (unsigned n = 0; n < raw_ndof; n++)
        {
          unsigned global_eqn = elem_pt->eqn_number(n);
          jacobian(raw_ndof, raw_ndof + 1 + n) =
            Phi[global_eqn] / Count[global_eqn];
        }

        //   //Loop over the dofs
        //     for(unsigned n=0;n<augmented_ndof;n++)
        //      {
        //       unsigned long global_eqn = eqn_number(elem_pt,n);
        //       double* unknown_pt = Problem_pt->Dof_pt[global_eqn];
        //       double init = *unknown_pt;
        //       *unknown_pt += FD_step;

        //       //Get the new residuals
        //       get_residuals(elem_pt,newres);

        //       for(unsigned m=0;m<augmented_ndof;m++)
        //        {
        //         jacobian(m,n) = (newres[m] - residuals[m])/FD_step;
        //        }
        //       //Reset the unknown
        //       *unknown_pt = init;
        //      }
      }
      break;

      default:
        std::ostringstream error_stream;
        error_stream
          << "The Solve_which_system flag can only take values 0, 1, 2"
          << " not " << Solve_which_system << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //====================================================================
  /// Formulate the derivatives of the augmented system with respect
  /// to a parameter
  //====================================================================
  void FoldHandler::get_dresiduals_dparameter(
    GeneralisedElement* const& elem_pt,
    double* const& parameter_pt,
    Vector<double>& dres_dparam)
  {
    // Need to get raw residuals and jacobian
    unsigned raw_ndof = elem_pt->ndof();

    // Find out which system we are solving
    switch (Solve_which_system)
    {
        // If we are solving the standard system
      case Block_J:
      {
        // Get the basic residual derivatives
        elem_pt->get_dresiduals_dparameter(parameter_pt, dres_dparam);
      }
      break;

      // If we are solving the augmented-by-one system
      case Block_augmented_J:
      {
        // Get the basic residual derivatives
        elem_pt->get_dresiduals_dparameter(parameter_pt, dres_dparam);

        // Zero the final derivative
        dres_dparam[raw_ndof] = 0.0;
      }
      break;

      // If we are solving the full augmented system
      case Full_augmented:
      {
        DenseMatrix<double> djac_dparam(raw_ndof);
        // Get the basic residuals and jacobian derivatives initially
        elem_pt->get_djacobian_dparameter(
          parameter_pt, dres_dparam, djac_dparam);

        // The normalisation equation does not depend on the parameter
        dres_dparam[raw_ndof] = 0.0;

        // Now assemble the equations dJy/dparameter = 0
        for (unsigned i = 0; i < raw_ndof; i++)
        {
          dres_dparam[raw_ndof + 1 + i] = 0.0;
          for (unsigned j = 0; j < raw_ndof; j++)
          {
            dres_dparam[raw_ndof + 1 + i] +=
              djac_dparam(i, j) * Y[elem_pt->eqn_number(j)];
          }
        }
      }
      break;

      default:
        std::ostringstream error_stream;
        error_stream
          << "The Solve_which_system flag can only take values 0, 1, 2"
          << " not " << Solve_which_system << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //========================================================================
  /// Overload the derivative of the residuals and jacobian
  /// with respect to a parameter so that it breaks because it should not
  /// be required
  //========================================================================
  void FoldHandler::get_djacobian_dparameter(GeneralisedElement* const& elem_pt,
                                             double* const& parameter_pt,
                                             Vector<double>& dres_dparam,
                                             DenseMatrix<double>& djac_dparam)
  {
    std::ostringstream error_stream;
    error_stream
      << "This function has not been implemented because it is not required\n";
    error_stream << "in standard problems.\n";
    error_stream
      << "If you find that you need it, you will have to implement it!\n\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


  //=====================================================================
  /// Overload the hessian vector product function so that
  /// it breaks because it should not be required
  //========================================================================
  void FoldHandler::get_hessian_vector_products(
    GeneralisedElement* const& elem_pt,
    Vector<double> const& Y,
    DenseMatrix<double> const& C,
    DenseMatrix<double>& product)
  {
    std::ostringstream error_stream;
    error_stream
      << "This function has not been implemented because it is not required\n";
    error_stream << "in standard problems.\n";
    error_stream
      << "If you find that you need it, you will have to implement it!\n\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


  //==========================================================================
  /// Return the eigenfunction(s) associated with the bifurcation that
  /// has been detected in bifurcation tracking problems
  //==========================================================================
  void FoldHandler::get_eigenfunction(Vector<DoubleVector>& eigenfunction)
  {
    // There is only one (real) null vector
    eigenfunction.resize(1);
    LinearAlgebraDistribution dist(Problem_pt->communicator_pt(), Ndof, false);
    // Rebuild the vector
    eigenfunction[0].build(&dist, 0.0);
    // Set the value to be the null vector
    for (unsigned n = 0; n < Ndof; n++)
    {
      eigenfunction[0][n] = Y[n];
    }
  }


  //=======================================================================
  /// Destructor return the problem to its original (non-augmented) state
  //=======================================================================
  FoldHandler::~FoldHandler()
  {
    // If we are using the block solver reset the problem's linear solver
    // to the original one
    AugmentedBlockFoldLinearSolver* block_fold_solver_pt =
      dynamic_cast<AugmentedBlockFoldLinearSolver*>(
        Problem_pt->linear_solver_pt());

    if (block_fold_solver_pt)
    {
      // Reset the problem's linear solver
      Problem_pt->linear_solver_pt() = block_fold_solver_pt->linear_solver_pt();
      // Delete the block solver
      delete block_fold_solver_pt;
    }

    // Resize the number of dofs
    Problem_pt->Dof_pt.resize(Ndof);
    Problem_pt->Dof_distribution_pt->build(
      Problem_pt->communicator_pt(), Ndof, false);
    // Remove all previous sparse storage used during Jacobian assembly
    Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);
  }

  //====================================================================
  /// Set to solve the augmented-by-one block system.
  //===================================================================
  void FoldHandler::solve_augmented_block_system()
  {
    // Only bother to do anything if we haven't already set the flag
    if (Solve_which_system != Block_augmented_J)
    {
      // If we were solving the system with the original jacobian add the
      // parameter
      if (Solve_which_system == Block_J)
      {
        Problem_pt->Dof_pt.push_back(Parameter_pt);
      }

      // Restrict the problem to the standard variables and
      // the bifurcation parameter only
      Problem_pt->Dof_pt.resize(Ndof + 1);
      Problem_pt->Dof_distribution_pt->build(
        Problem_pt->communicator_pt(), Ndof + 1, false);
      // Remove all previous sparse storage used during Jacobian assembly
      Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);
      // Set the solve flag
      Solve_which_system = Block_augmented_J;
    }
  }


  //====================================================================
  /// Set to solve the block system. The boolean flag specifies
  /// whether the block decomposition should use exactly the same jacobian
  //===================================================================
  void FoldHandler::solve_block_system()
  {
    // Only bother to do anything if we haven't already set the flag
    if (Solve_which_system != Block_J)
    {
      // Restrict the problem to the standard variables
      Problem_pt->Dof_pt.resize(Ndof);
      Problem_pt->Dof_distribution_pt->build(
        Problem_pt->communicator_pt(), Ndof, false);
      // Remove all previous sparse storage used during Jacobian assembly
      Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

      // Set the solve flag
      Solve_which_system = Block_J;
    }
  }

  //=================================================================
  /// Set to Solve non-block system
  //=================================================================
  void FoldHandler::solve_full_system()
  {
    // Only do something if we are not solving the full system
    if (Solve_which_system != Full_augmented)
    {
      // If we were solving the problem without any augmentation,
      // add the parameter again
      if (Solve_which_system == Block_J)
      {
        Problem_pt->Dof_pt.push_back(Parameter_pt);
      }

      // Always add the additional unknowns back into the problem
      for (unsigned n = 0; n < Ndof; n++)
      {
        Problem_pt->Dof_pt.push_back(&Y[n]);
      }

      // update the Dof distribution pt
      Problem_pt->Dof_distribution_pt->build(
        Problem_pt->communicator_pt(), Ndof * 2 + 1, false);
      // Remove all previous sparse storage used during Jacobian assembly
      Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

      Solve_which_system = Full_augmented;
    }
  }


  //======================================================================
  /// Clean up the memory that may have been allocated by the solver
  //=====================================================================
  BlockPitchForkLinearSolver::~BlockPitchForkLinearSolver()
  {
    if (B_pt != 0)
    {
      delete B_pt;
    }
    if (C_pt != 0)
    {
      delete C_pt;
    }
    if (D_pt != 0)
    {
      delete D_pt;
    }
    if (dJy_dparam_pt != 0)
    {
      delete dJy_dparam_pt;
    }
  }

  //===================================================================
  /// Use a block factorisation to solve the augmented system
  /// associated with a PitchFork bifurcation.
  //===================================================================
  void BlockPitchForkLinearSolver::solve(Problem* const& problem_pt,
                                         DoubleVector& result)
  {
    std::cout << "Block pitchfork solve" << std::endl;
    // Locally cache the pointer to the handler.
    PitchForkHandler* handler_pt =
      static_cast<PitchForkHandler*>(problem_pt->assembly_handler_pt());

    // Get the augmented distribution from the handler and use it as
    // the distribution of the linear solver
    LinearAlgebraDistribution aug_dist =
      handler_pt->Augmented_dof_distribution_pt;
    this->build_distribution(aug_dist);

    // If the result vector is not setup then die
    if (!result.built())
    {
      throw OomphLibError("Result vector must be built\n",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Store the distribution of the result, which is probably not in
    // the natural distribution of the augmented solver
    LinearAlgebraDistribution result_dist(result.distribution_pt());

    // Locally cache a pointer to the parameter
    double* const parameter_pt = handler_pt->Parameter_pt;

    // Firstly get the derivatives of the full augmented system with
    // respect to the parameter

    // Allocate storage for the derivatives of the residuals with respect
    // to the global parameter using the augmented distribution
    DoubleVector dRdparam;
    // Then get the appropriate derivatives which will come back in
    // some distribution
    problem_pt->get_derivative_wrt_global_parameter(parameter_pt, dRdparam);
    // Redistribute into the augmented distribution
    dRdparam.redistribute(&aug_dist);

    // Switch the handler to "block solver" mode (sort out distribution)
    handler_pt->solve_block_system();

    // Temporary vector for the result (I shouldn't have to set this up)
    DoubleVector x1;

    // We are going to do resolves using the underlying linear solver
    Linear_solver_pt->enable_resolve();
    // Solve the first (standard) system Jx1 = R
    Linear_solver_pt->solve(problem_pt, x1);

    // Allocate storage for B, C and D which can be used in the resolve
    // and the derivatives of the jacobian/eigenvector product with
    // respect to the parameter
    if (B_pt != 0)
    {
      delete B_pt;
    }
    B_pt = new DoubleVector(Linear_solver_pt->distribution_pt(), 0.0);
    if (C_pt != 0)
    {
      delete C_pt;
    }
    C_pt = new DoubleVector(Linear_solver_pt->distribution_pt(), 0.0);
    if (D_pt != 0)
    {
      delete D_pt;
    }
    D_pt = new DoubleVector(Linear_solver_pt->distribution_pt(), 0.0);
    // Need this to be in the distribution of the dofs
    if (dJy_dparam_pt != 0)
    {
      delete dJy_dparam_pt;
    }
    dJy_dparam_pt = new DoubleVector(handler_pt->Dof_distribution_pt, 0.0);

    // Get the symmetry vector from the handler should have the
    // Dof distribution
    DoubleVector psi = handler_pt->Psi;

    // Temporary vector for the rhs that is dR/dparam (augmented distribution)
    // DoubleVector F(Linear_solver_pt->distribution_pt(),0.0);
    DoubleVector F(handler_pt->Dof_distribution_pt);
    const unsigned n_dof_local = F.nrow_local();

    // f.nrow local copied from dRdparam
    for (unsigned n = 0; n < n_dof_local; n++)
    {
      F[n] = dRdparam[n];
    }
    // Fill in the rhs that is dJy/dparam  //dJy_dparam nrow local

    // In standard cases the offset here will be one
    unsigned offset = 1;
#ifdef OOMPH_HAS_MPI
    // If we are distributed and not on the first processor
    // then there is no offset and the eigenfunction is
    // directly after the standard dofs
    int my_rank = problem_pt->communicator_pt()->my_rank();
    if ((problem_pt->distributed()) && (my_rank != 0))
    {
      offset = 0;
    }
#endif
    for (unsigned n = 0; n < n_dof_local; n++)
    {
      (*dJy_dparam_pt)[n] = dRdparam[n_dof_local + offset + n];
    }

    // Now resolve to find c and d
    // First, redistribute F and psi into the linear algebra distribution
    F.redistribute(Linear_solver_pt->distribution_pt());
    psi.redistribute(Linear_solver_pt->distribution_pt());

    Linear_solver_pt->resolve(F, *C_pt);
    Linear_solver_pt->resolve(psi, *D_pt);

    // We can now construct various dot products
    double psi_d = psi.dot(*D_pt);
    double psi_c = psi.dot(*C_pt);
    double psi_x1 = psi.dot(x1);

    // Calculate another intermediate constant
    const double Psi = psi_d / psi_c;

    // Construct the second intermediate value,
    // assumption is that result has been
    // set to the current value of the residuals
    // First, redistribute into the Natural distribution of the augmented system
    result.redistribute(&aug_dist);
    // The required parameter is that corresponding to the dof, which
    // is only stored on the root processor
    double parameter_residual = result[n_dof_local];
#ifdef OOMPH_HAS_MPI
    // Broadcast to all others, if we have a distributed problem
    if (problem_pt->distributed())
    {
      MPI_Bcast(&parameter_residual,
                1,
                MPI_DOUBLE,
                0,
                problem_pt->communicator_pt()->mpi_comm());
    }
#endif
    // Now construct the value
    double x2 = (psi_x1 - parameter_residual) / psi_c;

    // Now construct the vectors that multiply the jacobian terms
    Vector<DoubleVectorWithHaloEntries> D_and_X1(2);
    D_and_X1[0].build(Linear_solver_pt->distribution_pt(), 0.0);
    D_and_X1[1].build(Linear_solver_pt->distribution_pt(), 0.0);
    // Fill in the appropriate terms
    // Get the number of local dofs from the Linear_solver_pt distribution
    const unsigned n_dof_local_linear_solver =
      Linear_solver_pt->distribution_pt()->nrow_local();
    for (unsigned n = 0; n < n_dof_local_linear_solver; n++)
    {
      const double C_ = (*C_pt)[n];
      D_and_X1[0][n] = (*D_pt)[n] - Psi * C_;
      D_and_X1[1][n] = x1[n] - x2 * C_;
    }

    // Local storage for the result of the product terms
    Vector<DoubleVectorWithHaloEntries> Jprod_D_and_X1(2);

    // Redistribute the inputs to have the Dof distribution pt
    D_and_X1[0].redistribute(handler_pt->Dof_distribution_pt);
    D_and_X1[1].redistribute(handler_pt->Dof_distribution_pt);

    // Set up the halo scheme
#ifdef OOMPH_HAS_MPI
    D_and_X1[0].build_halo_scheme(problem_pt->Halo_scheme_pt);
    D_and_X1[1].build_halo_scheme(problem_pt->Halo_scheme_pt);
#endif

    // Get the products from the new problem function
    problem_pt->get_hessian_vector_products(
      handler_pt->Y, D_and_X1, Jprod_D_and_X1);

    // OK, now we can formulate the next vectors
    //(again assuming result contains residuals)
    // Need to redistribute F to the dof distribution
    // F.redistribute(handler_pt->Dof_distribution_pt);
    DoubleVector G(handler_pt->Dof_distribution_pt);

    for (unsigned n = 0; n < n_dof_local; n++)
    {
      G[n] = result[n_dof_local + offset + n] - Jprod_D_and_X1[1][n] -
             x2 * dRdparam[n_dof_local + offset + n];
      Jprod_D_and_X1[0][n] *= -1.0;
      Jprod_D_and_X1[0][n] -= Psi * dRdparam[n_dof_local + offset + n];
    }


    // Then redistribute back to the linear solver distribution
    G.redistribute(Linear_solver_pt->distribution_pt());
    Jprod_D_and_X1[0].redistribute(Linear_solver_pt->distribution_pt());
    Jprod_D_and_X1[1].redistribute(Linear_solver_pt->distribution_pt());

    // Linear solve to get B
    Linear_solver_pt->resolve(Jprod_D_and_X1[0], *B_pt);
    // Liner solve to get x3
    DoubleVector x3(Linear_solver_pt->distribution_pt(), 0.0);
    Linear_solver_pt->resolve(G, x3);

    // Construst a couple of additional products
    double l_x3 = psi.dot(x3);
    double l_b = psi.dot(*B_pt);

    // get the last intermediate variable
    // The required parameter is that corresponding to the dof, which
    // is only stored on the root processor
    double sigma_residual = result[2 * (n_dof_local + offset) - 1];
#ifdef OOMPH_HAS_MPI
    // Broadcast to all others, if we have a distributed problem
    if (problem_pt->distributed())
    {
      MPI_Bcast(&sigma_residual,
                1,
                MPI_DOUBLE,
                0,
                problem_pt->communicator_pt()->mpi_comm());
    }
#endif


    const double delta_sigma = (l_x3 - sigma_residual) / l_b;
    const double delta_lambda = x2 - delta_sigma * Psi;

    // Need to do some more rearrangements here because result is global
    // but the other vectors are not!

    // Create temporary DoubleVectors to hold the results
    DoubleVector res1(Linear_solver_pt->distribution_pt());
    DoubleVector res2(Linear_solver_pt->distribution_pt());

    for (unsigned n = 0; n < n_dof_local_linear_solver; n++)
    {
      res1[n] = x1[n] - delta_lambda * (*C_pt)[n] - delta_sigma * (*D_pt)[n];
      res2[n] = x3[n] - delta_sigma * (*B_pt)[n];
    }

    // Now redistribute these into the Dof distribution pointer
    res1.redistribute(handler_pt->Dof_distribution_pt);
    res2.redistribute(handler_pt->Dof_distribution_pt);

    // Now can copy over into results into the result vector
    for (unsigned n = 0; n < n_dof_local; n++)
    {
      result[n] = res1[n];
      result[n_dof_local + offset + n] = res2[n];
    }

    // Add the final contributions to the residuals
    // only on the root processor if we have a distributed problem
#ifdef OOMPH_HAS_MPI
    if ((!problem_pt->distributed()) || (my_rank == 0))
#endif
    {
      result[n_dof_local] = delta_lambda;
      result[2 * n_dof_local + 1] = delta_sigma;
    }


    // The sign of the determinant is given by the sign of
    // the product psi_c and l_b
    // NOT CHECKED YET!
    problem_pt->sign_of_jacobian() =
      static_cast<int>(std::fabs(psi_c * l_b) / (psi_c * l_b));

    // Redistribute the result into its incoming distribution
    result.redistribute(&result_dist);

    // Switch things to our block solver
    handler_pt->solve_full_system();

    // If we are not storing things, clear the results
    if (!Enable_resolve)
    {
      // We no longer need to store the matrix
      Linear_solver_pt->disable_resolve();
      delete B_pt;
      B_pt = 0;
      delete C_pt;
      C_pt = 0;
      delete D_pt;
      D_pt = 0;
      delete dJy_dparam_pt;
      dJy_dparam_pt = 0;
    }
    // Otherwise also store the pointer to the problem
    else
    {
      Problem_pt = problem_pt;
    }
  }


  //==============================================================
  // Hack the re-solve to use the block factorisation
  //==============================================================
  void BlockPitchForkLinearSolver::resolve(const DoubleVector& rhs,
                                           DoubleVector& result)
  {
    std::cout << "Block pitchfork resolve" << std::endl;
    // Check that the factors have been stored
    if (B_pt == 0)
    {
      throw OomphLibError("The required vectors have not been stored",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }


    // Cache pointer to the problem
    Problem* const problem_pt = Problem_pt;

    // Locally cache pointer to the handler
    PitchForkHandler* handler_pt =
      static_cast<PitchForkHandler*>(problem_pt->assembly_handler_pt());

    // Get the augmented distribution from the handler
    LinearAlgebraDistribution aug_dist =
      handler_pt->Augmented_dof_distribution_pt;
    // this->build_distribution(aug_dist);

    // Find the number of dofs of the augmented system
    // const unsigned n_aug_dof = problem_pt->ndof();

    // Storage for the result distribution
    LinearAlgebraDistribution result_dist;

    // if the result vector is not setup then rebuild with distribution
    // = natural distribution of augmented solver
    if (!result.built())
    {
      result.build(&aug_dist, 0.0);
    }
    // Otherwise store the incoming distribution and redistribute
    else
    {
      result_dist.build(result.distribution_pt());
      result.redistribute(&aug_dist);
    }


    // Locally cache a pointer to the parameter
    // double* const parameter_pt = handler_pt->Parameter_pt;

    // Switch things to our block solver
    handler_pt->solve_block_system();
    // We need to find out the number of dofs
    // unsigned n_dof = problem_pt->ndof();

    // create the linear algebra distribution for this solver
    // currently only global (non-distributed) distributions are allowed
    // LinearAlgebraDistribution
    // dist(problem_pt->communicator_pt(),n_dof,false);
    // this->build_distribution(dist);

    // if the result vector is not setup then rebuild with distribution = global
    // if (!result.built())
    // {
    //  result.build(this->distribution_pt(),0.0);
    // }


    // Copy the rhs into local storage
    DoubleVector rhs_local = rhs;
    // and redistribute into the augmented distribution
    rhs_local.redistribute(&aug_dist);

    // Setup storage for output
    DoubleVector x1(Linear_solver_pt->distribution_pt(), 0.0);
    DoubleVector x3(Linear_solver_pt->distribution_pt(), 0.0);

    // Create input for RHS with the natural distribution
    DoubleVector input_x1(handler_pt->Dof_distribution_pt);
    const unsigned n_dof_local = input_x1.nrow_local();

    // Set the values of the a vector
    for (unsigned n = 0; n < n_dof_local; n++)
    {
      input_x1[n] = rhs_local[n];
    }
    // Need to redistribute into the linear algebra distribution
    input_x1.redistribute(Linear_solver_pt->distribution_pt());

    // We want to resolve, of course
    Linear_solver_pt->enable_resolve();
    // Now solve for the first vector
    Linear_solver_pt->resolve(input_x1, x1);

    // Get the symmetry vector from the handler
    DoubleVector psi = handler_pt->Psi;
    // redistribute local copy
    psi.redistribute(Linear_solver_pt->distribution_pt());

    // We can now construct various dot products
    double psi_d = psi.dot(*D_pt);
    double psi_c = psi.dot(*C_pt);
    double psi_x1 = psi.dot(x1);

    // Calculate another intermediate constant
    const double Psi = psi_d / psi_c;

    // Construct the second intermediate value,
    // assumption is that rhs has been set to the current value of the residuals
    double parameter_residual = rhs_local[n_dof_local];
#ifdef OOMPH_HAS_MPI
    // Broadcast to all others, if we have a distributed problem
    if (problem_pt->distributed())
    {
      MPI_Bcast(&parameter_residual,
                1,
                MPI_DOUBLE,
                0,
                problem_pt->communicator_pt()->mpi_comm());
    }
#endif

    double x2 = (psi_x1 - parameter_residual) / psi_c;

    // Now construct the vectors that multiply the jacobian terms
    // Vector<double> X1(n_dof/*+1*/);
    Vector<DoubleVectorWithHaloEntries> X1(1);
    X1[0].build(Linear_solver_pt->distribution_pt(), 0.0);

    const unsigned n_dof_local_linear_solver =
      Linear_solver_pt->distribution_pt()->nrow_local();
    for (unsigned n = 0; n < n_dof_local_linear_solver; n++)
    {
      X1[0][n] = x1[n] - x2 * (*C_pt)[n];
    }

    // Local storage for the product term
    Vector<DoubleVectorWithHaloEntries> Jprod_X1(1);

    // Redistribute the inputs to have the Dof distribution pt
    X1[0].redistribute(handler_pt->Dof_distribution_pt);

    // Set up the halo scheme
#ifdef OOMPH_HAS_MPI
    X1[0].build_halo_scheme(problem_pt->Halo_scheme_pt);
#endif

    // Get the product from the problem
    problem_pt->get_hessian_vector_products(handler_pt->Y, X1, Jprod_X1);

    // In standard cases the offset here will be one
    unsigned offset = 1;
#ifdef OOMPH_HAS_MPI
    // If we are distributed and not on the first processor
    // then there is no offset and the eigenfunction is
    // directly after the standard dofs
    int my_rank = problem_pt->communicator_pt()->my_rank();
    if ((problem_pt->distributed()) && (my_rank != 0))
    {
      offset = 0;
    }
#endif

    // OK, now we can formulate the next vectors
    //(again assuming result contains residuals)
    // Local storage for the product terms
    DoubleVector Mod_Jprod_X1(handler_pt->Dof_distribution_pt, 0.0);

    for (unsigned n = 0; n < n_dof_local; n++)
    {
      Mod_Jprod_X1[n] = rhs_local[n_dof_local + offset + n] - Jprod_X1[0][n] -
                        x2 * (*dJy_dparam_pt)[n];
    }

    // Redistribute back to the linear solver distribution
    Mod_Jprod_X1.redistribute(Linear_solver_pt->distribution_pt());

    // Liner solve to get x3
    Linear_solver_pt->resolve(Mod_Jprod_X1, x3);

    // Construst a couple of additional products
    double l_x3 = psi.dot(x3);
    double l_b = psi.dot(*B_pt);

    // get the last intermediate variable
    // The required parameter is that corresponding to the dof, which
    // is only stored on the root processor
    double sigma_residual = rhs_local[2 * (n_dof_local + offset) - 1];
#ifdef OOMPH_HAS_MPI
    // Broadcast to all others, if we have a distributed problem
    if (problem_pt->distributed())
    {
      MPI_Bcast(&sigma_residual,
                1,
                MPI_DOUBLE,
                0,
                problem_pt->communicator_pt()->mpi_comm());
    }
#endif

    // get the last intermediate variable
    const double delta_sigma = (l_x3 - sigma_residual) / l_b;
    const double delta_lambda = x2 - delta_sigma * Psi;

    // Create temporary DoubleVectors to hold the results
    DoubleVector res1(Linear_solver_pt->distribution_pt());
    DoubleVector res2(Linear_solver_pt->distribution_pt());

    for (unsigned n = 0; n < n_dof_local_linear_solver; n++)
    {
      res1[n] = x1[n] - delta_lambda * (*C_pt)[n] - delta_sigma * (*D_pt)[n];
      res2[n] = x3[n] - delta_sigma * (*B_pt)[n];
    }

    // Now redistribute these into the Dof distribution pointer
    res1.redistribute(handler_pt->Dof_distribution_pt);
    res2.redistribute(handler_pt->Dof_distribution_pt);

    // Now can copy over into results into the result vector
    for (unsigned n = 0; n < n_dof_local; n++)
    {
      result[n] = res1[n];
      result[n_dof_local + offset + n] = res2[n];
    }

    // Add the final contributions to the residuals
    // only on the root processor if we have a distributed problem
#ifdef OOMPH_HAS_MPI
    if ((!problem_pt->distributed()) || (my_rank == 0))
#endif
    {
      result[n_dof_local] = delta_lambda;
      result[2 * n_dof_local + 1] = delta_sigma;
    }

    // Redistribute the result into its incoming distribution (if it had one)
    if (result_dist.built())
    {
      result.redistribute(&result_dist);
    }

    Linear_solver_pt->disable_resolve();

    // Switch things to our block solver
    handler_pt->solve_full_system();
  }


  //--------------------------------------------------------------


  //======================================================================
  /// Clean up the memory that may have been allocated by the solver
  //=====================================================================
  AugmentedBlockPitchForkLinearSolver::~AugmentedBlockPitchForkLinearSolver()
  {
    if (Alpha_pt != 0)
    {
      delete Alpha_pt;
    }
    if (E_pt != 0)
    {
      delete E_pt;
    }
  }

  //===================================================================
  /// Use a block factorisation to solve the augmented system
  /// associated with a PitchFork bifurcation.
  //===================================================================
  void AugmentedBlockPitchForkLinearSolver::solve(Problem* const& problem_pt,
                                                  DoubleVector& result)
  {
    std::cout << "Augmented pitchfork solve" << std::endl;
    // if the result is setup then it should not be distributed
#ifdef PARANOID
    if (result.built())
    {
      if (result.distributed())
      {
        throw OomphLibError("The result vector must not be distributed",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif


    // Locally cache the pointer to the handler.
    PitchForkHandler* handler_pt =
      static_cast<PitchForkHandler*>(problem_pt->assembly_handler_pt());

    // Switch the handler to "block solver" mode
    handler_pt->solve_augmented_block_system();

    // We need to find out the number of dofs in the problem
    unsigned n_dof = problem_pt->ndof();

    // create the linear algebra distribution for this solver
    // currently only global (non-distributed) distributions are allowed
    LinearAlgebraDistribution dist(problem_pt->communicator_pt(), n_dof, false);
    this->build_distribution(dist);

    // if the result vector is not setup then rebuild with distribution = global
    if (!result.built())
    {
      result.build(this->distribution_pt(), 0.0);
    }

    // Setup storage for temporary vectors
    DoubleVector a(this->distribution_pt(), 0.0),
      b(this->distribution_pt(), 0.0);

    // Allocate storage for Alpha which can be used in the resolve
    if (Alpha_pt != 0)
    {
      delete Alpha_pt;
    }
    Alpha_pt = new DoubleVector(this->distribution_pt(), 0.0);

    // We are going to do resolves using the underlying linear solver
    Linear_solver_pt->enable_resolve();
    // Solve the first system Aa = R
    Linear_solver_pt->solve(problem_pt, a);

    // Get the symmetry vector from the handler
    DoubleVector psi(this->distribution_pt(), 0.0);
    for (unsigned n = 0; n < (n_dof - 1); ++n)
    {
      psi[n] = handler_pt->Psi[n];
    }
    // Set the final entry to zero
    psi[n_dof - 1] = 0.0;

    // Now resolve to find alpha
    Linear_solver_pt->resolve(psi, *Alpha_pt);

    // We can now construct our multipliers
    // Prepare to scale
    double dof_length = 0.0, a_length = 0.0, alpha_length = 0.0;
    for (unsigned n = 0; n < n_dof; n++)
    {
      if (std::fabs(problem_pt->dof(n)) > dof_length)
      {
        dof_length = std::fabs(problem_pt->dof(n));
      }
      if (std::fabs(a[n]) > a_length)
      {
        a_length = std::fabs(a[n]);
      }
      if (std::fabs((*Alpha_pt)[n]) > alpha_length)
      {
        alpha_length = std::fabs((*Alpha_pt)[n]);
      }
    }

    double a_mult = dof_length / a_length;
    double alpha_mult = dof_length / alpha_length;
    const double FD_step = 1.0e-8;
    a_mult += FD_step;
    alpha_mult += FD_step;
    a_mult *= FD_step;
    alpha_mult *= FD_step;

    // Local storage for the product terms
    DoubleVector Jprod_a(this->distribution_pt(), 0.0),
      Jprod_alpha(this->distribution_pt(), 0.0);

    // Calculate the product of the jacobian matrices, etc
    unsigned long n_element = problem_pt->mesh_pt()->nelement();
    for (unsigned long e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the ndofs in each element
      unsigned n_var = handler_pt->ndof(elem_pt);
      // Get the jacobian matrices
      DenseMatrix<double> jac(n_var), jac_a(n_var), jac_alpha(n_var);
      // the elemental residual
      Vector<double> res_elemental(n_var);
      // Get unperturbed jacobian
      handler_pt->get_jacobian(elem_pt, res_elemental, jac);

      // Backup the dofs
      Vector<double> dof_bac(n_var);
      // Perturb the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        dof_bac[n] = problem_pt->dof(eqn_number);
        // Pertub by vector a
        problem_pt->dof(eqn_number) += a_mult * a[eqn_number];
      }

      problem_pt->actions_after_change_in_bifurcation_parameter();

      // Now get the new jacobian
      handler_pt->get_jacobian(elem_pt, res_elemental, jac_a);

      // Perturb the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        problem_pt->dof(eqn_number) = dof_bac[n];
        // Pertub by vector a
        problem_pt->dof(eqn_number) += alpha_mult * (*Alpha_pt)[eqn_number];
      }

      problem_pt->actions_after_change_in_bifurcation_parameter();

      // Now get the new jacobian
      handler_pt->get_jacobian(elem_pt, res_elemental, jac_alpha);

      // Reset the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        problem_pt->dof(eqn_number) = dof_bac[n];
      }

      problem_pt->actions_after_change_in_bifurcation_parameter();

      // OK, now work out the products
      for (unsigned n = 0; n < (n_var - 1); n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        double prod_a = 0.0, prod_alpha = 0.0;
        for (unsigned m = 0; m < (n_var - 1); m++)
        {
          unsigned unknown = handler_pt->eqn_number(elem_pt, m);
          prod_a += (jac_a(n, m) - jac(n, m)) * handler_pt->Y[unknown];
          prod_alpha += (jac_alpha(n, m) - jac(n, m)) * handler_pt->Y[unknown];
        }
        Jprod_a[eqn_number] += prod_a / a_mult;
        Jprod_alpha[eqn_number] += prod_alpha / alpha_mult;
      }
    }

    Jprod_alpha[n_dof - 1] = 0.0;
    Jprod_a[n_dof - 1] = 0.0;

    // OK, now we can formulate the next vectors
    // The assumption here is that the result has been set to the
    // residuals.
    for (unsigned n = 0; n < n_dof - 1; n++)
    {
      b[n] = result[n_dof + n] - Jprod_a[n];
    }
    b[n_dof - 1] = result[2 * n_dof - 1];


    // Allocate storage for E which can be used in the resolve
    if (E_pt != 0)
    {
      delete E_pt;
    }
    E_pt = new DoubleVector(this->distribution_pt(), 0.0);

    DoubleVector f(this->distribution_pt(), 0.0);

    Linear_solver_pt->resolve(b, f);
    Linear_solver_pt->resolve(Jprod_alpha, *E_pt);

    // Calculate the final entry in the vector e
    const double e_final = (*E_pt)[n_dof - 1];
    // Calculate the final entry in the vector d
    const double d_final = -f[n_dof - 1] / e_final;
    // Assemble the final corrections
    for (unsigned n = 0; n < n_dof - 1; n++)
    {
      result[n] = a[n] - (*Alpha_pt)[n] * d_final;
      result[n_dof + n] = f[n] + (*E_pt)[n] * d_final;
    }

    result[n_dof - 1] = a[n_dof - 1] - (*Alpha_pt)[n_dof - 1] * d_final;
    result[2 * n_dof - 1] = d_final;

    // The sign of the jacobian is the sign of the final entry in e
    problem_pt->sign_of_jacobian() =
      static_cast<int>(std::fabs(e_final) / e_final);


    // Switch things to our block solver
    handler_pt->solve_full_system();

    // If we are not storing things, clear the results
    if (!Enable_resolve)
    {
      // We no longer need to store the matrix
      Linear_solver_pt->disable_resolve();
      delete Alpha_pt;
      Alpha_pt = 0;
      delete E_pt;
      E_pt = 0;
    }
    // Otherwise also store the pointer to the problem
    else
    {
      Problem_pt = problem_pt;
    }
  }


  //==============================================================
  // Hack the re-solve to use the block factorisation
  //==============================================================
  void AugmentedBlockPitchForkLinearSolver::resolve(const DoubleVector& rhs,
                                                    DoubleVector& result)
  {
    std::cout << "Augmented pitchfork resolve" << std::endl;
    // Check that the factors have been stored
    if (Alpha_pt == 0)
    {
      throw OomphLibError("The required vectors have not been stored",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Get the pointer to the problem
    Problem* const problem_pt = Problem_pt;

    PitchForkHandler* handler_pt =
      static_cast<PitchForkHandler*>(problem_pt->assembly_handler_pt());

    // Switch things to our block solver
    handler_pt->solve_augmented_block_system();
    // We need to find out the number of dofs
    unsigned n_dof = problem_pt->ndof();

    // create the linear algebra distribution for this solver
    // currently only global (non-distributed) distributions are allowed
    LinearAlgebraDistribution dist(problem_pt->communicator_pt(), n_dof, false);
    this->build_distribution(dist);

    // if the result vector is not setup then rebuild with distribution = global
    if (!result.built())
    {
      result.build(this->distribution_pt(), 0.0);
    }


    // Setup storage
    DoubleVector a(this->distribution_pt(), 0.0),
      b(this->distribution_pt(), 0.0);

    // Set the values of the a vector
    for (unsigned n = 0; n < n_dof; n++)
    {
      a[n] = rhs[n];
    }

    Linear_solver_pt->enable_resolve();

    // Copy rhs vector into local storage so it doesn't get overwritten
    // if the linear solver decides to initialise the solution vector, say,
    // which it's quite entitled to do!
    DoubleVector input_a(a);

    Linear_solver_pt->resolve(input_a, a);

    // We can now construct our multipliers
    // Prepare to scale
    double dof_length = 0.0, a_length = 0.0;
    for (unsigned n = 0; n < n_dof; n++)
    {
      if (std::fabs(problem_pt->dof(n)) > dof_length)
      {
        dof_length = std::fabs(problem_pt->dof(n));
      }

      if (std::fabs(a[n]) > a_length)
      {
        a_length = std::fabs(a[n]);
      }
    }
    double a_mult = dof_length / a_length;
    const double FD_step = 1.0e-8;
    a_mult += FD_step;
    a_mult *= FD_step;

    DoubleVector Jprod_a(this->distribution_pt(), 0.0);

    unsigned long n_element = problem_pt->mesh_pt()->nelement();
    for (unsigned long e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the ndofs in each element
      unsigned n_var = handler_pt->ndof(elem_pt);
      // Get some jacobian matrices
      DenseMatrix<double> jac(n_var), jac_a(n_var);
      // the elemental residual
      Vector<double> res_elemental(n_var);
      // Get unperturbed jacobian
      handler_pt->get_jacobian(elem_pt, res_elemental, jac);

      // Backup the dofs
      DoubleVector dof_bac(this->distribution_pt(), 0.0);
      // Perturb the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        dof_bac[n] = problem_pt->dof(eqn_number);
        // Pertub by vector a
        problem_pt->dof(eqn_number) += a_mult * a[eqn_number];
      }

      problem_pt->actions_after_change_in_bifurcation_parameter();

      // Now get the new jacobian
      handler_pt->get_jacobian(elem_pt, res_elemental, jac_a);

      // Reset the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        problem_pt->dof(eqn_number) = dof_bac[n];
      }

      problem_pt->actions_after_change_in_bifurcation_parameter();

      // OK, now work out the products
      for (unsigned n = 0; n < (n_var - 1); n++)
      {
        unsigned eqn_number = handler_pt->eqn_number(elem_pt, n);
        double prod_a = 0.0;
        for (unsigned m = 0; m < (n_var - 1); m++)
        {
          unsigned unknown = handler_pt->eqn_number(elem_pt, m);
          prod_a += (jac_a(n, m) - jac(n, m)) * handler_pt->Y[unknown];
        }
        Jprod_a[eqn_number] += prod_a / a_mult;
      }
    }

    Jprod_a[n_dof - 1] = 0.0;

    // OK, now we can formulate the next vectors
    for (unsigned n = 0; n < n_dof - 1; n++)
    {
      b[n] = rhs[n_dof + n] - Jprod_a[n];
    }
    b[n_dof - 1] = rhs[2 * n_dof - 1];

    DoubleVector f(this->distribution_pt(), 0.0);

    Linear_solver_pt->resolve(b, f);

    // Calculate the final entry in the vector d
    const double d_final = -f[n_dof - 1] / (*E_pt)[n_dof - 1];
    // Assemble the final corrections
    for (unsigned n = 0; n < n_dof - 1; n++)
    {
      result[n] = a[n] - (*Alpha_pt)[n] * d_final;
      result[n_dof + n] = f[n] + (*E_pt)[n] * d_final;
    }

    result[n_dof - 1] = a[n_dof - 1] - (*Alpha_pt)[n_dof - 1] * d_final;
    result[2 * n_dof - 1] = d_final;

    Linear_solver_pt->disable_resolve();

    // Switch things to our block solver
    handler_pt->solve_full_system();
  }

  ///////////////////////////////////////////////////////////////////////
  // Non-inline functions for the PitchForkHandler class
  //////////////////////////////////////////////////////////////////////

  //==================================================================
  /// Constructor: Initialise the PitchForkHandler by setting intial
  /// guesses for Sigma, Y, specifying C and Psi and calculating count.
  //==================================================================
  PitchForkHandler::PitchForkHandler(
    Problem* const& problem_pt,
    AssemblyHandler* const& assembly_handler_pt,
    double* const& parameter_pt,
    const DoubleVector& symmetry_vector)
    : Solve_which_system(Full_augmented), Sigma(0.0), Parameter_pt(parameter_pt)
  {
    // Set the problem pointer
    Problem_pt = problem_pt;
    // Set the assembly handler
    Assembly_handler_pt = assembly_handler_pt;
    // Set the number of degrees of freedom
    Ndof = Problem_pt->ndof();
    // Backup the distribution
    Dof_distribution_pt = Problem_pt->Dof_distribution_pt;
#ifdef OOMPH_HAS_MPI
    // Set the distribution flag
    Distributed = Problem_pt->distributed();
#endif

    // Find the elements in the problem
    unsigned n_element = Problem_pt->mesh_pt()->nelement();

#ifdef OOMPH_HAS_MPI

    // Work out all the global equations to which this processor
    // contributes and store the halo data in the problem
    Problem_pt->setup_dof_halo_scheme();

#endif


    // Now use the dof distribution for all double vectors
    Psi.build(Dof_distribution_pt);
    Y.build(Dof_distribution_pt);
    C.build(Dof_distribution_pt);
    Count.build(Dof_distribution_pt);
#ifdef OOMPH_HAS_MPI
    Count.build_halo_scheme(Problem_pt->Halo_scheme_pt);
#endif

    // Loop over the elements and count the entries
    // and number of (non-halo) elements
    unsigned n_non_halo_element_local = 0;
    for (unsigned e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = Problem_pt->mesh_pt()->element_pt(e);
#ifdef OOMPH_HAS_MPI
      // Ignore halo elements
      if (!elem_pt->is_halo())
      {
#endif
        // Increment the number of non halo elements
        ++n_non_halo_element_local;
        // Now count the number of times the element contributes to a value
        unsigned n_var = assembly_handler_pt->ndof(elem_pt);
        for (unsigned n = 0; n < n_var; n++)
        {
          ++Count.global_value(assembly_handler_pt->eqn_number(elem_pt, n));
        }
#ifdef OOMPH_HAS_MPI
      }
#endif
    }

    // Add together all the counts
#ifdef OOMPH_HAS_MPI
    Count.sum_all_halo_and_haloed_values();

    // If distributed, find the total number of elements in the problem
    if (Distributed)
    {
      // Need to gather the total number of non halo elements
      MPI_Allreduce(&n_non_halo_element_local,
                    &Nelement,
                    1,
                    MPI_UNSIGNED,
                    MPI_SUM,
                    Problem_pt->communicator_pt()->mpi_comm());
    }
    // Otherwise the total number is the same on each processor
    else
#endif
    {
      Nelement = n_non_halo_element_local;
    }

#ifdef OOMPH_HAS_MPI
    // Only add the parameter to the first processor, if distributed
    int my_rank = Problem_pt->communicator_pt()->my_rank();
    if ((!Distributed) || (my_rank == 0))
#endif
    {
      // Add the parameter to the problem
      Problem_pt->Dof_pt.push_back(parameter_pt);
    }

    // Find length of the symmetry vector
    double length = symmetry_vector.norm();

    // Add the unknowns for the null vector to the problem
    // Normalise the symmetry vector and initialise the null vector to the
    // symmetry vector and set the constant vector, c, to the
    // normalised symmetry vector.

    // Need to redistribute the symmetry vector into the (natural)
    // distribution of the Dofs
#ifdef OOMPH_HAS_MPI
    // Check that the symmetry vector has the correct distribution
    if (!(*symmetry_vector.distribution_pt() == *Dof_distribution_pt))
    {
      throw OomphLibError(
        "The symmetry vector must have the same distribution as the dofs\n",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Only loop over the local unknowns
    const unsigned n_dof_local = Dof_distribution_pt->nrow_local();
    for (unsigned n = 0; n < n_dof_local; n++)
    {
      Problem_pt->Dof_pt.push_back(&Y[n]);
      // Psi[n] = symmetry_vector[n];
      Psi[n] = Y[n] = C[n] = symmetry_vector[n] / length;
    }


#ifdef OOMPH_HAS_MPI
    // Set up the required halo schemes (which also synchronises)
    Psi.build_halo_scheme(Problem_pt->Halo_scheme_pt);
    Y.build_halo_scheme(Problem_pt->Halo_scheme_pt);
    C.build_halo_scheme(Problem_pt->Halo_scheme_pt);


    if ((!Distributed) || (my_rank == 0))
#endif
    // Add the slack parameter to the problem on the first processor
    // if distributed
    {
      Problem_pt->Dof_pt.push_back(&Sigma);
    }

#ifdef OOMPH_HAS_MPI
    if (Distributed)
    {
      unsigned augmented_first_row = 0;
      unsigned augmented_n_row_local = 0;

      // Set up the translation scheme on every processor
      Global_eqn_number.resize(2 * Ndof + 2);
      int n_proc = Problem_pt->communicator_pt()->nproc();
      unsigned global_eqn_count = 0;
      for (int d = 0; d < n_proc; d++)
      {
        // Find out the first row of the current processor
        if (my_rank == d)
        {
          augmented_first_row = global_eqn_count;
        }

        const unsigned n_row_local = Dof_distribution_pt->nrow_local(d);
        const unsigned first_row = Dof_distribution_pt->first_row(d);
        // Add the basic equations
        for (unsigned n = 0; n < n_row_local; n++)
        {
          Global_eqn_number[first_row + n] = global_eqn_count;
          ++global_eqn_count;
        }
        // If on the first processor add the pointer to the parameter
        if (d == 0)
        {
          Global_eqn_number[Ndof] = global_eqn_count;
          ++global_eqn_count;
        }
        // Add the eigenfunction
        for (unsigned n = 0; n < n_row_local; n++)
        {
          Global_eqn_number[Ndof + 1 + first_row + n] = global_eqn_count;
          ++global_eqn_count;
        }
        // Finally add the slack parameter
        if (d == 0)
        {
          Global_eqn_number[2 * Ndof + 1] = global_eqn_count;
          ++global_eqn_count;
        }
        // Find out the number of rows of the current processor
        if (my_rank == d)
        {
          augmented_n_row_local = global_eqn_count - augmented_first_row;
        }
      }

      // Make a new linear algebra distribution
      Augmented_dof_distribution_pt =
        new LinearAlgebraDistribution(Problem_pt->communicator_pt(),
                                      augmented_first_row,
                                      augmented_n_row_local);
    }
    else
#endif
    {
      Augmented_dof_distribution_pt = new LinearAlgebraDistribution(
        Problem_pt->communicator_pt(), 2 * Ndof + 2, false);
    }

    // resize
    // Problem_pt->Dof_distribution_pt->build(Problem_pt->communicator_pt(),
    //                                       2*Ndof+2,false);

    Problem_pt->Dof_distribution_pt = Augmented_dof_distribution_pt;

    // Remove all previous sparse storage used during Jacobian assembly
    Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);
  }

  //=====================================================================
  /// Destructor return the problem to its original (non-augmented) state
  //=====================================================================
  PitchForkHandler::~PitchForkHandler()
  {
    // If we are using the block solver reset the problem's linear solver
    // to the original one
    BlockPitchForkLinearSolver* block_pitchfork_solver_pt =
      dynamic_cast<BlockPitchForkLinearSolver*>(Problem_pt->linear_solver_pt());

    if (block_pitchfork_solver_pt)
    {
      // Reset the problem's linear solver
      Problem_pt->linear_solver_pt() =
        block_pitchfork_solver_pt->linear_solver_pt();
      // Delete the block solver
      delete block_pitchfork_solver_pt;
    }

    // If we are using the augmented
    // block solver reset the problem's linear solver
    // to the original one
    AugmentedBlockPitchForkLinearSolver* augmented_block_pitchfork_solver_pt =
      dynamic_cast<AugmentedBlockPitchForkLinearSolver*>(
        Problem_pt->linear_solver_pt());

    if (augmented_block_pitchfork_solver_pt)
    {
      // Reset the problem's linear solver
      Problem_pt->linear_solver_pt() =
        augmented_block_pitchfork_solver_pt->linear_solver_pt();
      // Delete the block solver
      delete augmented_block_pitchfork_solver_pt;
    }

    // Now return the problem to its original size
    Problem_pt->Dof_pt.resize(this->Dof_distribution_pt->nrow_local());

    // Problem_pt->Dof_pt.resize(Ndof);
    // Problem_pt->Dof_distribution_pt->build(Problem_pt->communicator_pt(),
    //                                       Ndof,false);

    Problem_pt->Dof_distribution_pt = this->Dof_distribution_pt;
    // Remove all previous sparse storage used during Jacobian assembly
    Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);
  }

  //================================================================
  /// Get the number of elemental degrees of freedom
  //================================================================
  unsigned PitchForkHandler::ndof(GeneralisedElement* const& elem_pt)
  {
    unsigned raw_ndof = elem_pt->ndof();

    // Return different values depending on the type of block decomposition
    switch (Solve_which_system)
    {
      case Full_augmented:
        return (2 * raw_ndof + 2);
        break;

      case Block_augmented_J:
        return (raw_ndof + 1);
        break;

      case Block_J:
        return raw_ndof;
        break;

      default:
        std::ostringstream error_stream;
        error_stream
          << "The Solve_which_system flag can only take values 0, 1, 2"
          << " not " << Solve_which_system << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  /// Get the global equation number of the local unknown
  unsigned long PitchForkHandler::eqn_number(GeneralisedElement* const& elem_pt,
                                             const unsigned& ieqn_local)
  {
    // Get the raw value
    unsigned raw_ndof = elem_pt->ndof();
    unsigned long global_eqn = 0;

    // Which system are we solving
    switch (Solve_which_system)
    {
        // In the block case, it's just the standard numbering
      case Block_J:
        global_eqn = elem_pt->eqn_number(ieqn_local);
        break;

        // Block augmented not done properly yet
      case Block_augmented_J:
        // Not done the distributed case
#ifdef OOMPH_HAS_MPI
        if (Distributed)
        {
          throw OomphLibError(
            "Block Augmented solver not implemented for distributed case\n",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Full case
      case Full_augmented:
        // The usual equations
        if (ieqn_local < raw_ndof)
        {
          global_eqn = this->global_eqn_number(elem_pt->eqn_number(ieqn_local));
        }
        // The bifurcation parameter equation
        else if (ieqn_local == raw_ndof)
        {
          global_eqn = this->global_eqn_number(Ndof);
        }
        // If we are assembling the full system we also have
        // The components of the null vector
        else if (ieqn_local < (2 * raw_ndof + 1))
        {
          global_eqn = this->global_eqn_number(
            Ndof + 1 + elem_pt->eqn_number(ieqn_local - 1 - raw_ndof));
        }
        // The slack parameter
        else
        {
          global_eqn = this->global_eqn_number(2 * Ndof + 1);
        }
        break;
    } // End of switch
    return global_eqn;
  }

  //==============================================================
  /// Get the residuals
  //==============================================================
  void PitchForkHandler::get_residuals(GeneralisedElement* const& elem_pt,
                                       Vector<double>& residuals)
  {
    // Need to get raw residuals and jacobian
    unsigned raw_ndof = elem_pt->ndof();

    // Find out which system we are solving
    switch (Solve_which_system)
    {
        // If we are solving the original system
      case Block_J:
      {
        // get the basic residuals
        elem_pt->get_residuals(residuals);
        // Now multiply to fill in the residuals for the final term
        for (unsigned i = 0; i < raw_ndof; i++)
        {
          unsigned local_eqn = elem_pt->eqn_number(i);
          // Add the slack parameter to the final residuals
          residuals[i] +=
            Sigma * Psi.global_value(local_eqn) / Count.global_value(local_eqn);
        }
      }
      break;

      // If we are solving the augmented-by-one system
      case Block_augmented_J:
      {
        // Get the basic residuals
        elem_pt->get_residuals(residuals);

        // Zero the final residual
        residuals[raw_ndof] = 0.0;
        // Now multiply to fill in the residuals for the final term
        for (unsigned i = 0; i < raw_ndof; i++)
        {
          unsigned local_eqn = elem_pt->eqn_number(i);
          // Add the slack parameter to the final residuals
          residuals[i] +=
            Sigma * Psi.global_value(local_eqn) / Count.global_value(local_eqn);
          // Final term that specifies the symmetry
          residuals[raw_ndof] += ((*Problem_pt->global_dof_pt(local_eqn)) *
                                  Psi.global_value(local_eqn)) /
                                 Count.global_value(local_eqn);
        }
      }
      break;

      // Otherwise we are solving the fully augemented system
      case Full_augmented:
      {
        DenseMatrix<double> jacobian(raw_ndof);

        // Get the basic residuals and jacobian
        elem_pt->get_jacobian(residuals, jacobian);

        // Initialise the final residuals
        residuals[raw_ndof] = 0.0;
        residuals[2 * raw_ndof + 1] = -1.0 / Nelement;

        // Now multiply to fill in the residuals associated
        // with the null vector condition
        for (unsigned i = 0; i < raw_ndof; i++)
        {
          unsigned local_eqn = elem_pt->eqn_number(i);
          residuals[raw_ndof + 1 + i] = 0.0;
          for (unsigned j = 0; j < raw_ndof; j++)
          {
            unsigned local_unknown = elem_pt->eqn_number(j);
            residuals[raw_ndof + 1 + i] +=
              jacobian(i, j) * Y.global_value(local_unknown);
          }

          // Add the slack parameter to the governing equations
          residuals[i] +=
            Sigma * Psi.global_value(local_eqn) / Count.global_value(local_eqn);

          // Specify the symmetry
          residuals[raw_ndof] += ((*Problem_pt->global_dof_pt(local_eqn)) *
                                  Psi.global_value(local_eqn)) /
                                 Count.global_value(local_eqn);
          // Specify the normalisation
          residuals[2 * raw_ndof + 1] +=
            (Y.global_value(local_eqn) * C.global_value(local_eqn)) /
            Count.global_value(local_eqn);
        }
      }
      break;

      default:
        std::ostringstream error_stream;
        error_stream
          << "The Solve_which_system flag can only take values 0, 1, 2"
          << " not " << Solve_which_system << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  //======================================================================
  /// Calculate the elemental Jacobian matrix "d equation
  /// / d variable".
  //======================================================================
  void PitchForkHandler::get_jacobian(GeneralisedElement* const& elem_pt,
                                      Vector<double>& residuals,
                                      DenseMatrix<double>& jacobian)
  {
    unsigned augmented_ndof = ndof(elem_pt);
    unsigned raw_ndof = elem_pt->ndof();

    // Which system are we solving
    switch (Solve_which_system)
    {
        // If we are solving the original system
      case Block_J:
      {
        // get the raw residuals and jacobian
        elem_pt->get_jacobian(residuals, jacobian);

        // Now multiply to fill in the residuals for the final term
        for (unsigned i = 0; i < raw_ndof; i++)
        {
          unsigned local_eqn = elem_pt->eqn_number(i);
          // Add the slack parameter to the final residuals
          residuals[i] +=
            Sigma * Psi.global_value(local_eqn) / Count.global_value(local_eqn);
        }
      }
      break;

      // If we are solving the augmented-by-one system
      case Block_augmented_J:
      {
        // Get the full residuals, we need them
        get_residuals(elem_pt, residuals);

        // Need to get the raw jacobian (and raw residuals)
        Vector<double> newres(augmented_ndof);
        elem_pt->get_jacobian(newres, jacobian);

        // Now do finite differencing stuff
        const double FD_step = 1.0e-8;
        // Fill in the first lot of finite differences
        {
          // increase the global parameter
          double* unknown_pt = Parameter_pt;
          double init = *unknown_pt;
          *unknown_pt += FD_step;

          // Not do any possible updates
          Problem_pt->actions_after_change_in_bifurcation_parameter();

          // Get the new (modified) residuals
          get_residuals(elem_pt, newres);

          // The final column  is given by the difference
          // between the residuals
          for (unsigned n = 0; n < raw_ndof; n++)
          {
            jacobian(n, augmented_ndof - 1) =
              (newres[n] - residuals[n]) / FD_step;
          }
          // Reset the global parameter
          *unknown_pt = init;

          // Now do any possible updates
          Problem_pt->actions_after_change_in_bifurcation_parameter();
        }

        // Fill in the bottom row
        for (unsigned n = 0; n < raw_ndof; n++)
        {
          unsigned local_eqn = elem_pt->eqn_number(n);
          jacobian(augmented_ndof - 1, n) =
            Psi.global_value(local_eqn) / Count.global_value(local_eqn);
        }
      }
      break;

      // Otherwise solving the full system
      case Full_augmented:
      {
        /// ALICE:
        Problem_pt->actions_before_newton_convergence_check();
        // Get the basic jacobian and residuals
        elem_pt->get_jacobian(residuals, jacobian);
        // get the proper residuals
        get_residuals(elem_pt, residuals);

        // Now fill in the next diagonal jacobian entry
        for (unsigned n = 0; n < raw_ndof; n++)
        {
          for (unsigned m = 0; m < raw_ndof; m++)
          {
            jacobian(raw_ndof + 1 + n, raw_ndof + 1 + m) = jacobian(n, m);
          }
          unsigned local_eqn = elem_pt->eqn_number(n);
          // Add in the sigma contribution
          jacobian(n, 2 * raw_ndof + 1) =
            Psi.global_value(local_eqn) / Count.global_value(local_eqn);
          // Symmetry constraint
          jacobian(raw_ndof, n) =
            Psi.global_value(local_eqn) / Count.global_value(local_eqn);
          // Non-zero constraint
          jacobian(2 * raw_ndof + 1, raw_ndof + 1 + n) =
            C.global_value(local_eqn) / Count.global_value(local_eqn);
        }

        // Finite difference the remaining blocks
        const double FD_step = 1.0e-8;

        Vector<double> newres_p(augmented_ndof);

        // Loop over the ndofs only
        for (unsigned n = 0; n < raw_ndof; ++n)
        {
          // Get the original (unaugmented) global equation number
          unsigned long global_eqn =
            Assembly_handler_pt->eqn_number(elem_pt, n);
          double* unknown_pt = Problem_pt->global_dof_pt(global_eqn);
          double init = *unknown_pt;
          *unknown_pt += FD_step;
          Problem_pt->actions_before_newton_convergence_check();
          // Get the new residuals
          get_residuals(elem_pt, newres_p);
          // Fill in the entries in the block d(Jy)/dx
          for (unsigned m = 0; m < raw_ndof; m++)
          {
            jacobian(raw_ndof + 1 + m, n) =
              (newres_p[raw_ndof + 1 + m] - residuals[raw_ndof + 1 + m]) /
              (FD_step);
          }

          // Reset the unknown
          *unknown_pt = init;
          Problem_pt->actions_before_newton_convergence_check();
        }

        {
          // Now increase the global parameter
          double* unknown_pt = Parameter_pt;
          double init = *unknown_pt;
          *unknown_pt += FD_step;

          Problem_pt->actions_after_change_in_bifurcation_parameter();
          // Get the new residuals
          get_residuals(elem_pt, newres_p);

          // Add in the first block d/dx
          for (unsigned m = 0; m < raw_ndof; m++)
          {
            jacobian(m, raw_ndof) = (newres_p[m] - residuals[m]) / FD_step;
          }
          // Add in the second block d/dy
          for (unsigned m = raw_ndof + 1; m < augmented_ndof - 1; m++)
          {
            jacobian(m, raw_ndof) = (newres_p[m] - residuals[m]) / FD_step;
          }

          // Reset the unknown
          *unknown_pt = init;
          Problem_pt->actions_after_change_in_bifurcation_parameter();
        }
      }
      break;

      default:
        std::ostringstream error_stream;
        error_stream
          << "The Solve_which_system flag can only take values 0, 1, 2"
          << " not " << Solve_which_system << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    /// ALICE
    Problem_pt->actions_after_change_in_bifurcation_parameter();
  }


  //==============================================================
  /// Get the derivatives of the residuals with respect to a parameter
  //==============================================================
  void PitchForkHandler::get_dresiduals_dparameter(
    GeneralisedElement* const& elem_pt,
    double* const& parameter_pt,
    Vector<double>& dres_dparam)
  {
    // Need to get raw residuals and jacobian
    unsigned raw_ndof = elem_pt->ndof();
    Problem_pt->actions_before_newton_convergence_check();
    // Find out which system we are solving
    switch (Solve_which_system)
    {
        // If we are solving the original system
      case Block_J:
      {
        // get the basic residual derivatives
        elem_pt->get_dresiduals_dparameter(parameter_pt, dres_dparam);
        // Slack parameter term does not depened explicitly on the parameter
        // so is not added
      }
      break;

      // If we are solving the augmented-by-one system
      case Block_augmented_J:
      {
        // Get the basic residual derivatives
        elem_pt->get_dresiduals_dparameter(parameter_pt, dres_dparam);

        // Zero the final residuals derivative
        dres_dparam[raw_ndof] = 0.0;

        // Other terms must not depend on the parameter
      }
      break;

      // Otherwise we are solving the fully augemented system
      case Full_augmented:
      {
        DenseMatrix<double> djac_dparam(raw_ndof);

        // Get the basic residuals and jacobian derivatives
        elem_pt->get_djacobian_dparameter(
          parameter_pt, dres_dparam, djac_dparam);

        // The "final" residual derivatives do not depend on the parameter
        dres_dparam[raw_ndof] = 0.0;
        dres_dparam[2 * raw_ndof + 1] = 0.0;

        // Now multiply to fill in the residuals associated
        // with the null vector condition
        for (unsigned i = 0; i < raw_ndof; i++)
        {
          dres_dparam[raw_ndof + 1 + i] = 0.0;
          for (unsigned j = 0; j < raw_ndof; j++)
          {
            unsigned local_unknown = elem_pt->eqn_number(j);
            dres_dparam[raw_ndof + 1 + i] +=
              djac_dparam(i, j) * Y.global_value(local_unknown);
          }
        }
      }
      break;

      default:
        std::ostringstream error_stream;
        error_stream
          << "The Solve_which_system flag can only take values 0, 1, 2"
          << " not " << Solve_which_system << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //========================================================================
  /// Overload the derivative of the residuals and jacobian
  /// with respect to a parameter so that it breaks because it should not
  /// be required
  //========================================================================
  void PitchForkHandler::get_djacobian_dparameter(
    GeneralisedElement* const& elem_pt,
    double* const& parameter_pt,
    Vector<double>& dres_dparam,
    DenseMatrix<double>& djac_dparam)
  {
    std::ostringstream error_stream;
    error_stream
      << "This function has not been implemented because it is not required\n";
    error_stream << "in standard problems.\n";
    error_stream
      << "If you find that you need it, you will have to implement it!\n\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


  //=====================================================================
  /// Overload the hessian vector product function so that
  /// it calls the underlying element's hessian function.
  //========================================================================
  void PitchForkHandler::get_hessian_vector_products(
    GeneralisedElement* const& elem_pt,
    Vector<double> const& Y,
    DenseMatrix<double> const& C,
    DenseMatrix<double>& product)
  {
    elem_pt->get_hessian_vector_products(Y, C, product);
  }


  //==========================================================================
  /// Return the eigenfunction(s) associated with the bifurcation that
  /// has been detected in bifurcation tracking problems
  //==========================================================================
  void PitchForkHandler::get_eigenfunction(Vector<DoubleVector>& eigenfunction)
  {
    // There is only one (real) null vector
    eigenfunction.resize(1);
    // Rebuild the vector
    eigenfunction[0].build(this->Dof_distribution_pt, 0.0);
    // Set the value to be the null vector
    const unsigned n_row_local = eigenfunction[0].nrow_local();
    for (unsigned n = 0; n < n_row_local; n++)
    {
      eigenfunction[0][n] = Y[n];
    }
  }

#ifdef OOMPH_HAS_MPI
  //=====================================================================
  // Synchronise the required data
  //====================================================================
  void PitchForkHandler::synchronise()
  {
    // Only need to bother if the problem is distributed
    if (Distributed)
    {
      // Need only to synchronise the eigenfunction
      Y.synchronise();
      // Also need to synchronise the parameter and the slack parameter
      double broadcast_data[2];
      broadcast_data[0] = *Parameter_pt;
      broadcast_data[1] = Sigma;
      // Broadcast from the root to all processors
      MPI_Bcast(broadcast_data,
                2,
                MPI_DOUBLE,
                0,
                Problem_pt->communicator_pt()->mpi_comm());

      // Now copy received the values back
      *Parameter_pt = broadcast_data[0];
      Sigma = broadcast_data[1];
    }
  }
#endif

  //====================================================================
  /// Set to solve the augmented-by-one block system.
  //===================================================================
  void PitchForkHandler::solve_augmented_block_system()
  {
    // Only bother to do anything if we haven't already set the flag
    if (Solve_which_system != Block_augmented_J)
    {
      // If we were solving the system with the original jacobian add the
      // parameter back
      if (Solve_which_system == Block_J)
      {
        Problem_pt->Dof_pt.push_back(Parameter_pt);
      }

      // Restrict the problem to the standard variables and
      // the bifurcation parameter only
      Problem_pt->Dof_distribution_pt = new LinearAlgebraDistribution(
        Problem_pt->communicator_pt(), Ndof + 1, false);
      Problem_pt->Dof_pt.resize(Ndof + 1);

      // Problem_pt->Dof_pt.resize(Ndof+1);
      //    Problem_pt->Dof_distribution_pt->build(Problem_pt->communicator_pt(),
      //                                         Ndof+1,false);

      // Remove all previous sparse storage used during Jacobian assembly
      Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

      // Set the solve flag
      Solve_which_system = Block_augmented_J;
    }
  }


  //==============================================================
  /// Set to solve the block system. The boolean flag is used
  /// to specify whether the block decomposition should use exactly
  /// the same jacobian as the original system.
  //==============================================================
  void PitchForkHandler::solve_block_system()
  {
    // Only bother to do anything if we haven't already set the flag
    if (Solve_which_system != Block_J)
    {
      // Restrict the problem to the standard variables
      Problem_pt->Dof_pt.resize(this->Dof_distribution_pt->nrow_local());
      Problem_pt->Dof_distribution_pt = this->Dof_distribution_pt;

      // Problem_pt->Dof_distribution_pt->build(Problem_pt->communicator_pt(),
      // Ndof,false);

      // Remove all previous sparse storage used during Jacobian assembly
      Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

      // Set the solve flag
      Solve_which_system = Block_J;
    }
  }

  //==============================================================
  /// Solve non-block system
  //==============================================================
  void PitchForkHandler::solve_full_system()
  {
    // Only do something if we are not solving the full system
    if (Solve_which_system != Full_augmented)
    {
#ifdef OOMPH_HAS_MPI
      int my_rank = Problem_pt->communicator_pt()->my_rank();
#endif
      // If we were solving the problem without any augementation
      // add the parameter again
      if (Solve_which_system == Block_J)
      {
#ifdef OOMPH_HAS_MPI

        if ((!Distributed) || (my_rank == 0))
#endif
        {
          Problem_pt->Dof_pt.push_back(Parameter_pt);
        }
      }

      // Always add the additional unknowns back into the problem
      const unsigned n_dof_local = Dof_distribution_pt->nrow_local();
      for (unsigned n = 0; n < n_dof_local; n++)
      {
        Problem_pt->Dof_pt.push_back(&Y[n]);
      }


#ifdef OOMPH_HAS_MPI
      if ((!Distributed) || (my_rank == 0))
#endif
      // Add the slack parameter to the problem on the first processor
      // if distributed
      {
        Problem_pt->Dof_pt.push_back(&Sigma);
      }


      // Delete the distribtion created in the augmented block solve
      if (Problem_pt->Dof_distribution_pt != this->Dof_distribution_pt)
      {
        delete Problem_pt->Dof_distribution_pt;
      }

      // update the Dof distribution
      Problem_pt->Dof_distribution_pt = this->Augmented_dof_distribution_pt;
      // Problem_pt->Dof_distribution_pt->build(Problem_pt->communicator_pt(),
      //                                       Ndof*2+2,false);
      // Remove all previous sparse storage used during Jacobian assembly
      Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

      // Set the solve flag
      Solve_which_system = Full_augmented;
    }
  }


  ///////////////////////////////////////////////////////////////////////
  // Non-inline functions for the BlockHopfLinearSolver class
  //////////////////////////////////////////////////////////////////////


  //======================================================================
  /// Clean up the memory that may have been allocated by the solver
  //=====================================================================
  BlockHopfLinearSolver::~BlockHopfLinearSolver()
  {
    if (A_pt != 0)
    {
      delete A_pt;
    }
    if (E_pt != 0)
    {
      delete E_pt;
    }
    if (G_pt != 0)
    {
      delete G_pt;
    }
  }

  //===================================================================
  /// Use a block factorisation to solve the augmented system
  /// associated with a Hopf bifurcation.
  //===================================================================
  void BlockHopfLinearSolver::solve(Problem* const& problem_pt,
                                    DoubleVector& result)
  {
    // if the result is setup then it should not be distributed
#ifdef PARANOID
    if (result.built())
    {
      if (result.distributed())
      {
        throw OomphLibError("The result vector must not be distributed",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Locally cache the pointer to the handler.
    HopfHandler* handler_pt =
      static_cast<HopfHandler*>(problem_pt->assembly_handler_pt());

    // Find the number of dofs in the augmented problem
    unsigned n_dof = problem_pt->ndof();

    // create the linear algebra distribution for this solver
    // currently only global (non-distributed) distributions are allowed
    LinearAlgebraDistribution dist(problem_pt->communicator_pt(), n_dof, false);
    this->build_distribution(dist);

    // Firstly, let's calculate the derivative of the residuals wrt
    // the parameter
    DoubleVector dRdparam(this->distribution_pt(), 0.0);

    const double FD_step = 1.0e-8;
    {
      // Cache pointer to the parameter (second last entry in the vector
      double* param_pt = &problem_pt->dof(n_dof - 2);
      // Backup the parameter
      double old_var = *param_pt;
      ;
      // Increment the parameter
      *param_pt += FD_step;
      problem_pt->actions_after_change_in_bifurcation_parameter();

      // Now calculate the new residuals
      problem_pt->get_residuals(dRdparam);

      // Now calculate the difference assume original residuals in resul
      for (unsigned n = 0; n < n_dof; n++)
      {
        dRdparam[n] = (dRdparam[n] - result[n]) / FD_step;
      }

      // Reset the parameter
      *param_pt = old_var;
      problem_pt->actions_after_change_in_bifurcation_parameter();
    }

    // Switch the handler to "standard" mode
    handler_pt->solve_standard_system();

    // Find out the number of dofs in the non-standard problem
    n_dof = problem_pt->ndof();

    // update the distribution pt
    dist.build(problem_pt->communicator_pt(), n_dof, false);
    this->build_distribution(dist);

    // Setup storage for temporary vector
    DoubleVector y1(this->distribution_pt(), 0.0),
      alpha(this->distribution_pt(), 0.0);

    // Allocate storage for A which can be used in the resolve
    if (A_pt != 0)
    {
      delete A_pt;
    }
    A_pt = new DoubleVector(this->distribution_pt(), 0.0);

    // We are going to do resolves using the underlying linear solver
    Linear_solver_pt->enable_resolve();

    // Solve the first system Jy1 = R
    Linear_solver_pt->solve(problem_pt, y1);

    // This should have set the sign of the jacobian, store it
    int sign_of_jacobian = problem_pt->sign_of_jacobian();

    // Now let's get the appropriate bit of alpha
    for (unsigned n = 0; n < n_dof; n++)
    {
      alpha[n] = dRdparam[n];
    }

    // Resolve to find A
    Linear_solver_pt->resolve(alpha, *A_pt);

    // Now set to the complex system
    handler_pt->solve_complex_system();

    // update the distribution
    dist.build(problem_pt->communicator_pt(), n_dof * 2, false);
    this->build_distribution(dist);

    // Resize the stored vector G
    if (G_pt != 0)
    {
      delete G_pt;
    }
    G_pt = new DoubleVector(this->distribution_pt(), 0.0);

    // Solve the first Zg = (Mz -My)
    Linear_solver_pt->solve(problem_pt, *G_pt);

    // This should have set the sign of the complex matrix's determinant,
    // multiply
    sign_of_jacobian *= problem_pt->sign_of_jacobian();

    // We can now construct our multipliers
    // Prepare to scale
    double dof_length = 0.0, a_length = 0.0, y1_length = 0.0;
    // Loop over the standard number of dofs
    for (unsigned n = 0; n < n_dof; n++)
    {
      if (std::fabs(problem_pt->dof(n)) > dof_length)
      {
        dof_length = std::fabs(problem_pt->dof(n));
      }
      if (std::fabs((*A_pt)[n]) > a_length)
      {
        a_length = std::fabs((*A_pt)[n]);
      }
      if (std::fabs(y1[n]) > y1_length)
      {
        y1_length = std::fabs(y1[n]);
      }
    }

    double a_mult = dof_length / a_length;
    double y1_mult = dof_length / y1_length;
    a_mult += FD_step;
    y1_mult += FD_step;
    a_mult *= FD_step;
    y1_mult *= FD_step;

    // Local storage for the product terms
    Vector<double> Jprod_a(2 * n_dof, 0.0), Jprod_y1(2 * n_dof, 0.0);

    // Temporary storage
    Vector<double> rhs(2 * n_dof);

    // find the number of elements
    unsigned long n_element = problem_pt->mesh_pt()->nelement();

    // Calculate the product of the jacobian matrices, etc
    for (unsigned long e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the ndofs in each element
      unsigned n_var = elem_pt->ndof();
      // Get the jacobian matrices
      DenseMatrix<double> jac(n_var), jac_a(n_var), jac_y1(n_var);
      DenseMatrix<double> M(n_var), M_a(n_var), M_y1(n_var);
      // Get unperturbed jacobian
      elem_pt->get_jacobian_and_mass_matrix(rhs, jac, M);

      // Backup the dofs
      Vector<double> dof_bac(n_var);
      // Perturb the original dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = elem_pt->eqn_number(n);
        dof_bac[n] = problem_pt->dof(eqn_number);
        // Pertub by vector A
        problem_pt->dof(eqn_number) += a_mult * (*A_pt)[eqn_number];
      }

      // Now get the new jacobian and mass matrix
      elem_pt->get_jacobian_and_mass_matrix(rhs, jac_a, M_a);

      // Perturb the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = elem_pt->eqn_number(n);
        problem_pt->dof(eqn_number) = dof_bac[n];
        // Pertub by vector y1
        problem_pt->dof(eqn_number) += y1_mult * y1[eqn_number];
      }

      // Now get the new jacobian and mass matrix
      elem_pt->get_jacobian_and_mass_matrix(rhs, jac_y1, M_y1);

      // Reset the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = elem_pt->eqn_number(n);
        problem_pt->dof(eqn_number) = dof_bac[n];
      }

      // OK, now work out the products
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = elem_pt->eqn_number(n);
        double prod_a1 = 0.0, prod_y11 = 0.0;
        double prod_a2 = 0.0, prod_y12 = 0.0;
        for (unsigned m = 0; m < n_var; m++)
        {
          const unsigned unknown = elem_pt->eqn_number(m);
          const double y = handler_pt->Phi[unknown];
          const double z = handler_pt->Psi[unknown];
          const double omega = handler_pt->Omega;
          // Real part (first line)
          prod_a1 +=
            (jac_a(n, m) - jac(n, m)) * y + omega * (M_a(n, m) - M(n, m)) * z;
          prod_y11 +=
            (jac_y1(n, m) - jac(n, m)) * y + omega * (M_y1(n, m) - M(n, m)) * z;
          // Imag par (second line)
          prod_a2 +=
            (jac_a(n, m) - jac(n, m)) * z - omega * (M_a(n, m) - M(n, m)) * y;
          prod_y12 +=
            (jac_y1(n, m) - jac(n, m)) * z - omega * (M_y1(n, m) - M(n, m)) * y;
        }
        Jprod_a[eqn_number] += prod_a1 / a_mult;
        Jprod_y1[eqn_number] += prod_y11 / y1_mult;
        Jprod_a[n_dof + eqn_number] += prod_a2 / a_mult;
        Jprod_y1[n_dof + eqn_number] += prod_y12 / y1_mult;
      }
    }

    // The assumption here is still that the result has been set to the
    // residuals.
    for (unsigned n = 0; n < 2 * n_dof; n++)
    {
      rhs[n] = result[n_dof + n] - Jprod_y1[n];
    }

    // Temporary storage
    DoubleVector y2(this->distribution_pt(), 0.0);

    // DoubleVector for rhs
    DoubleVector rhs2(this->distribution_pt(), 0.0);
    for (unsigned i = 0; i < 2 * n_dof; i++)
    {
      rhs2[i] = rhs[i];
    }

    // Solve it
    Linear_solver_pt->resolve(rhs2, y2);

    // Assemble the next RHS
    for (unsigned n = 0; n < 2 * n_dof; n++)
    {
      rhs[n] = dRdparam[n_dof + n] - Jprod_a[n];
    }

    // Resive the storage
    if (E_pt != 0)
    {
      delete E_pt;
    }
    E_pt = new DoubleVector(this->distribution_pt(), 0.0);

    // Solve for the next RHS
    for (unsigned i = 0; i < 2 * n_dof; i++)
    {
      rhs2[i] = rhs[i];
    }
    Linear_solver_pt->resolve(rhs2, *E_pt);

    // We can now calculate the final corrections
    // We need to work out a large number of dot products
    double dot_c = 0.0, dot_d = 0.0, dot_e = 0.0, dot_f = 0.0, dot_g = 0.0;
    double dot_h = 0.0;

    for (unsigned n = 0; n < n_dof; n++)
    {
      // Get the appopriate entry
      const double Cn = handler_pt->C[n];
      dot_c += Cn * y2[n];
      dot_d += Cn * y2[n_dof + n];
      dot_e += Cn * (*E_pt)[n];
      dot_f += Cn * (*E_pt)[n_dof + n];
      dot_g += Cn * (*G_pt)[n];
      dot_h += Cn * (*G_pt)[n_dof + n];
    }

    // Now we should be able to work out the corrections
    double denom = dot_e * dot_h - dot_g * dot_f;

    // Copy the previous residuals
    double R31 = result[3 * n_dof], R32 = result[3 * n_dof + 1];
    // Delta parameter
    const double delta_param =
      ((R32 - dot_d) * dot_g - (R31 - dot_c) * dot_h) / denom;
    // Delta frequency
    const double delta_w = -((R32 - dot_d) + dot_f * delta_param) / (dot_h);

    // Load into the result vector
    result[3 * n_dof] = delta_param;
    result[3 * n_dof + 1] = delta_w;

    // The corrections to the null vector
    for (unsigned n = 0; n < 2 * n_dof; n++)
    {
      result[n_dof + n] =
        y2[n] - (*E_pt)[n] * delta_param - (*G_pt)[n] * delta_w;
    }

    // Finally add the corrections to the unknowns
    for (unsigned n = 0; n < n_dof; n++)
    {
      result[n] = y1[n] - (*A_pt)[n] * delta_param;
    }

    // The sign of the jacobian is the previous signs multiplied by the
    // sign of the denominator
    problem_pt->sign_of_jacobian() =
      sign_of_jacobian * static_cast<int>(std::fabs(denom) / denom);

    // Switch things to our full solver
    handler_pt->solve_full_system();

    // If we are not storing things, clear the results
    if (!Enable_resolve)
    {
      // We no longer need to store the matrix
      Linear_solver_pt->disable_resolve();
      delete A_pt;
      A_pt = 0;
      delete E_pt;
      E_pt = 0;
      delete G_pt;
      G_pt = 0;
    }
    // Otherwise, also store the problem pointer
    else
    {
      Problem_pt = problem_pt;
    }
  }

  //=====================================================================
  /// Solve for two right hand sides, required for (efficient) continuation
  /// because otherwise we have to store the inverse of the jacobian
  /// and the complex equivalent ... nasty.
  //=====================================================================
  void BlockHopfLinearSolver::solve_for_two_rhs(Problem* const& problem_pt,
                                                DoubleVector& result,
                                                const DoubleVector& rhs2,
                                                DoubleVector& result2)
  {
    // if the result is setup then it should not be distributed
#ifdef PARANOID
    if (result.built())
    {
      if (result.distributed())
      {
        throw OomphLibError("The result vector must not be distributed",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
    if (result2.built())
    {
      if (result2.distributed())
      {
        throw OomphLibError("The result2 vector must not be distributed",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Locally cache the pointer to the handler.
    HopfHandler* handler_pt =
      static_cast<HopfHandler*>(problem_pt->assembly_handler_pt());

    // Find the number of dofs in the augmented problem
    unsigned n_dof = problem_pt->ndof();

    // create the linear algebra distribution for this solver
    // currently only global (non-distributed) distributions are allowed
    LinearAlgebraDistribution dist(problem_pt->communicator_pt(), n_dof, false);
    this->build_distribution(dist);

    // if the result vector is not setup then rebuild with distribution = global
    if (!result.built())
    {
      result.build(this->distribution_pt(), 0.0);
    }
    if (!result2.built())
    {
      result2.build(this->distribution_pt(), 0.0);
    }

    // Firstly, let's calculate the derivative of the residuals wrt
    // the parameter
    DoubleVector dRdparam(this->distribution_pt(), 0.0);

    const double FD_step = 1.0e-8;
    {
      // Cache pointer to the parameter (second last entry in the vector
      double* param_pt = &problem_pt->dof(n_dof - 2);
      // Backup the parameter
      double old_var = *param_pt;
      ;
      // Increment the parameter
      *param_pt += FD_step;
      problem_pt->actions_after_change_in_bifurcation_parameter();

      // Now calculate the new residuals
      problem_pt->get_residuals(dRdparam);

      // Now calculate the difference assume original residuals in resul
      for (unsigned n = 0; n < n_dof; n++)
      {
        dRdparam[n] = (dRdparam[n] - result[n]) / FD_step;
      }

      // Reset the parameter
      *param_pt = old_var;
      problem_pt->actions_after_change_in_bifurcation_parameter();
    }

    // Switch the handler to "standard" mode
    handler_pt->solve_standard_system();

    // Find out the number of dofs in the non-standard problem
    n_dof = problem_pt->ndof();

    //
    dist.build(problem_pt->communicator_pt(), n_dof, false);
    this->build_distribution(dist);

    // Setup storage for temporary vector
    DoubleVector y1(this->distribution_pt(), 0.0),
      alpha(this->distribution_pt(), 0.0),
      y1_resolve(this->distribution_pt(), 0.0);

    // Allocate storage for A which can be used in the resolve
    if (A_pt != 0)
    {
      delete A_pt;
    }
    A_pt = new DoubleVector(this->distribution_pt(), 0.0);

    // We are going to do resolves using the underlying linear solver
    Linear_solver_pt->enable_resolve();

    // Solve the first system Jy1 =
    Linear_solver_pt->solve(problem_pt, y1);

    // This should have set the sign of the jacobian, store it
    int sign_of_jacobian = problem_pt->sign_of_jacobian();

    // Now let's get the appropriate bit of alpha
    for (unsigned n = 0; n < n_dof; n++)
    {
      alpha[n] = dRdparam[n];
    }

    // Resolve to find A
    Linear_solver_pt->resolve(alpha, *A_pt);

    // Get the solution for the second rhs
    for (unsigned n = 0; n < n_dof; n++)
    {
      alpha[n] = rhs2[n];
    }

    // Resolve to find y1_resolve
    Linear_solver_pt->resolve(alpha, y1_resolve);

    // Now set to the complex system
    handler_pt->solve_complex_system();

    // rebuild the Distribution
    dist.build(problem_pt->communicator_pt(), n_dof * 2, false);
    this->build_distribution(dist);

    // Resize the stored vector G
    if (G_pt != 0)
    {
      delete G_pt;
    }
    G_pt = new DoubleVector(this->distribution_pt(), 0.0);

    // Solve the first Zg = (Mz -My)
    Linear_solver_pt->solve(problem_pt, *G_pt);

    // This should have set the sign of the complex matrix's determinant,
    // multiply
    sign_of_jacobian *= problem_pt->sign_of_jacobian();

    // We can now construct our multipliers
    // Prepare to scale
    double dof_length = 0.0, a_length = 0.0, y1_length = 0.0,
           y1_resolve_length = 0.0;
    // Loop over the standard number of dofs
    for (unsigned n = 0; n < n_dof; n++)
    {
      if (std::fabs(problem_pt->dof(n)) > dof_length)
      {
        dof_length = std::fabs(problem_pt->dof(n));
      }
      if (std::fabs((*A_pt)[n]) > a_length)
      {
        a_length = std::fabs((*A_pt)[n]);
      }
      if (std::fabs(y1[n]) > y1_length)
      {
        y1_length = std::fabs(y1[n]);
      }
      if (std::fabs(y1_resolve[n]) > y1_resolve_length)
      {
        y1_resolve_length = std::fabs(y1[n]);
      }
    }


    double a_mult = dof_length / a_length;
    double y1_mult = dof_length / y1_length;
    double y1_resolve_mult = dof_length / y1_resolve_length;
    a_mult += FD_step;
    y1_mult += FD_step;
    y1_resolve_mult += FD_step;
    a_mult *= FD_step;
    y1_mult *= FD_step;
    y1_resolve_mult *= FD_step;

    // Local storage for the product terms
    Vector<double> Jprod_a(2 * n_dof, 0.0), Jprod_y1(2 * n_dof, 0.0);
    Vector<double> Jprod_y1_resolve(2 * n_dof, 0.0);

    // Temporary storage
    Vector<double> rhs(2 * n_dof);

    // find the number of elements
    unsigned long n_element = problem_pt->mesh_pt()->nelement();

    // Calculate the product of the jacobian matrices, etc
    for (unsigned long e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the ndofs in each element
      unsigned n_var = elem_pt->ndof();
      // Get the jacobian matrices
      DenseMatrix<double> jac(n_var), jac_a(n_var), jac_y1(n_var),
        jac_y1_resolve(n_var);
      DenseMatrix<double> M(n_var), M_a(n_var), M_y1(n_var),
        M_y1_resolve(n_var);
      // Get unperturbed jacobian
      elem_pt->get_jacobian_and_mass_matrix(rhs, jac, M);

      // Backup the dofs
      Vector<double> dof_bac(n_var);
      // Perturb the original dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = elem_pt->eqn_number(n);
        dof_bac[n] = problem_pt->dof(eqn_number);
        // Pertub by vector A
        problem_pt->dof(eqn_number) += a_mult * (*A_pt)[eqn_number];
      }

      // Now get the new jacobian and mass matrix
      elem_pt->get_jacobian_and_mass_matrix(rhs, jac_a, M_a);

      // Perturb the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = elem_pt->eqn_number(n);
        problem_pt->dof(eqn_number) = dof_bac[n];
        // Pertub by vector y1
        problem_pt->dof(eqn_number) += y1_mult * y1[eqn_number];
      }

      // Now get the new jacobian and mass matrix
      elem_pt->get_jacobian_and_mass_matrix(rhs, jac_y1, M_y1);

      // Perturb the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = elem_pt->eqn_number(n);
        problem_pt->dof(eqn_number) = dof_bac[n];
        // Pertub by vector y1
        problem_pt->dof(eqn_number) += y1_resolve_mult * y1_resolve[eqn_number];
      }

      // Now get the new jacobian and mass matrix
      elem_pt->get_jacobian_and_mass_matrix(rhs, jac_y1_resolve, M_y1_resolve);

      // Reset the dofs
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = elem_pt->eqn_number(n);
        problem_pt->dof(eqn_number) = dof_bac[n];
      }

      // OK, now work out the products
      for (unsigned n = 0; n < n_var; n++)
      {
        unsigned eqn_number = elem_pt->eqn_number(n);
        double prod_a1 = 0.0, prod_y11 = 0.0, prod_y1_resolve1 = 0.0;
        double prod_a2 = 0.0, prod_y12 = 0.0, prod_y1_resolve2 = 0.0;
        for (unsigned m = 0; m < n_var; m++)
        {
          const unsigned unknown = elem_pt->eqn_number(m);
          const double y = handler_pt->Phi[unknown];
          const double z = handler_pt->Psi[unknown];
          const double omega = handler_pt->Omega;
          // Real part (first line)
          prod_a1 +=
            (jac_a(n, m) - jac(n, m)) * y + omega * (M_a(n, m) - M(n, m)) * z;
          prod_y11 +=
            (jac_y1(n, m) - jac(n, m)) * y + omega * (M_y1(n, m) - M(n, m)) * z;
          prod_y1_resolve1 += (jac_y1_resolve(n, m) - jac(n, m)) * y +
                              omega * (M_y1_resolve(n, m) - M(n, m)) * z;
          // Imag par (second line)
          prod_a2 +=
            (jac_a(n, m) - jac(n, m)) * z - omega * (M_a(n, m) - M(n, m)) * y;
          prod_y12 +=
            (jac_y1(n, m) - jac(n, m)) * z - omega * (M_y1(n, m) - M(n, m)) * y;
          prod_y1_resolve2 += (jac_y1_resolve(n, m) - jac(n, m)) * z -
                              omega * (M_y1_resolve(n, m) - M(n, m)) * y;
        }
        Jprod_a[eqn_number] += prod_a1 / a_mult;
        Jprod_y1[eqn_number] += prod_y11 / y1_mult;
        Jprod_y1_resolve[eqn_number] += prod_y1_resolve1 / y1_resolve_mult;
        Jprod_a[n_dof + eqn_number] += prod_a2 / a_mult;
        Jprod_y1[n_dof + eqn_number] += prod_y12 / y1_mult;
        Jprod_y1_resolve[n_dof + eqn_number] +=
          prod_y1_resolve2 / y1_resolve_mult;
      }
    }

    // The assumption here is still that the result has been set to the
    // residuals.
    for (unsigned n = 0; n < 2 * n_dof; n++)
    {
      rhs[n] = result[n_dof + n] - Jprod_y1[n];
    }

    // Temporary storage
    DoubleVector y2(this->distribution_pt(), 0.0);

    // Solve it
    DoubleVector temp_rhs(this->distribution_pt(), 0.0);
    for (unsigned i = 0; i < 2 * n_dof; i++)
    {
      temp_rhs[i] = rhs[i];
    }
    Linear_solver_pt->resolve(temp_rhs, y2);

    // Assemble the next RHS
    for (unsigned n = 0; n < 2 * n_dof; n++)
    {
      rhs[n] = dRdparam[n_dof + n] - Jprod_a[n];
    }

    // Resive the storage
    if (E_pt != 0)
    {
      delete E_pt;
    }
    E_pt = new DoubleVector(this->distribution_pt(), 0.0);

    // Solve for the next RHS
    for (unsigned i = 0; i < 2 * n_dof; i++)
    {
      temp_rhs[i] = rhs[i];
    }
    Linear_solver_pt->resolve(temp_rhs, *E_pt);

    // Assemble the next RHS
    for (unsigned n = 0; n < 2 * n_dof; n++)
    {
      rhs[n] = rhs2[n_dof + n] - Jprod_y1_resolve[n];
    }

    DoubleVector y2_resolve(this->distribution_pt(), 0.0);
    for (unsigned i = 0; i < 2 * n_dof; i++)
    {
      temp_rhs[i] = rhs[i];
    }
    Linear_solver_pt->resolve(temp_rhs, y2_resolve);


    // We can now calculate the final corrections
    // We need to work out a large number of dot products
    double dot_c = 0.0, dot_d = 0.0, dot_e = 0.0, dot_f = 0.0, dot_g = 0.0;
    double dot_h = 0.0;

    double dot_c_resolve = 0.0, dot_d_resolve = 0.0;

    for (unsigned n = 0; n < n_dof; n++)
    {
      // Get the appopriate entry
      const double Cn = handler_pt->C[n];
      dot_c += Cn * y2[n];
      dot_d += Cn * y2[n_dof + n];
      dot_c_resolve += Cn * y2_resolve[n];
      dot_d_resolve += Cn * y2_resolve[n_dof + n];
      dot_e += Cn * (*E_pt)[n];
      dot_f += Cn * (*E_pt)[n_dof + n];
      dot_g += Cn * (*G_pt)[n];
      dot_h += Cn * (*G_pt)[n_dof + n];
    }

    // Now we should be able to work out the corrections
    double denom = dot_e * dot_h - dot_g * dot_f;

    // Copy the previous residuals
    double R31 = result[3 * n_dof], R32 = result[3 * n_dof + 1];
    // Delta parameter
    const double delta_param =
      ((R32 - dot_d) * dot_g - (R31 - dot_c) * dot_h) / denom;
    // Delta frequency
    const double delta_w = -((R32 - dot_d) + dot_f * delta_param) / (dot_h);

    // Corrections
    double R31_resolve = rhs2[3 * n_dof], R32_resolve = rhs2[3 * n_dof + 1];
    // Delta parameter
    const double delta_param_resolve = ((R32_resolve - dot_d_resolve) * dot_g -
                                        (R31_resolve - dot_c_resolve) * dot_h) /
                                       denom;
    // Delta frequency
    const double delta_w_resolve =
      -((R32_resolve - dot_d_resolve) + dot_f * delta_param_resolve) / (dot_h);


    // Load into the result vector
    result[3 * n_dof] = delta_param;
    result[3 * n_dof + 1] = delta_w;

    // The corrections to the null vector
    for (unsigned n = 0; n < 2 * n_dof; n++)
    {
      result[n_dof + n] =
        y2[n] - (*E_pt)[n] * delta_param - (*G_pt)[n] * delta_w;
    }

    // Finally add the corrections to the unknowns
    for (unsigned n = 0; n < n_dof; n++)
    {
      result[n] = y1[n] - (*A_pt)[n] * delta_param;
    }

    // Load into the result vector
    result2[3 * n_dof] = delta_param_resolve;
    result2[3 * n_dof + 1] = delta_w_resolve;

    // The corrections to the null vector
    for (unsigned n = 0; n < 2 * n_dof; n++)
    {
      result2[n_dof + n] = y2_resolve[n] - (*E_pt)[n] * delta_param_resolve -
                           (*G_pt)[n] * delta_w_resolve;
    }

    // Finally add the corrections to the unknowns
    for (unsigned n = 0; n < n_dof; n++)
    {
      result2[n] = y1_resolve[n] - (*A_pt)[n] * delta_param_resolve;
    }


    // The sign of the jacobian is the previous signs multiplied by the
    // sign of the denominator
    problem_pt->sign_of_jacobian() =
      sign_of_jacobian * static_cast<int>(std::fabs(denom) / denom);

    // Switch things to our full solver
    handler_pt->solve_full_system();

    // If we are not storing things, clear the results
    if (!Enable_resolve)
    {
      // We no longer need to store the matrix
      Linear_solver_pt->disable_resolve();
      delete A_pt;
      A_pt = 0;
      delete E_pt;
      E_pt = 0;
      delete G_pt;
      G_pt = 0;
    }
    // Otherwise, also store the problem pointer
    else
    {
      Problem_pt = problem_pt;
    }
  }


  //======================================================================
  // Hack the re-solve to use the block factorisation
  //======================================================================
  void BlockHopfLinearSolver::resolve(const DoubleVector& rhs,
                                      DoubleVector& result)
  {
    throw OomphLibError("resolve() is not implemented for this solver",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }


  ///////////////////////////////////////////////////////////////////////
  // Non-inline functions for the HopfHandler class
  //////////////////////////////////////////////////////////////////////

  //====================================================================
  /// Constructor: Initialise the hopf handler,
  /// by setting initial guesses for Phi, Psi  and calculating Count.
  /// If the system changes, a new  handler must be constructed.
  //===================================================================
  HopfHandler::HopfHandler(Problem* const& problem_pt,
                           double* const& parameter_pt)
    : Solve_which_system(0), Parameter_pt(parameter_pt), Omega(0.0)
  {
    // Set the problem pointer
    Problem_pt = problem_pt;
    // Set the number of non-augmented degrees of freedom
    Ndof = problem_pt->ndof();

    // create the linear algebra distribution for this solver
    // currently only global (non-distributed) distributions are allowed
    LinearAlgebraDistribution* dist_pt =
      new LinearAlgebraDistribution(problem_pt->communicator_pt(), Ndof, false);

    // Resize the vectors of additional dofs
    Phi.resize(Ndof);
    Psi.resize(Ndof);
    C.resize(Ndof);
    Count.resize(Ndof, 0);

    // Loop over all the elements in the problem
    unsigned n_element = problem_pt->mesh_pt()->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the local freedoms in an element
      unsigned n_var = elem_pt->ndof();
      for (unsigned n = 0; n < n_var; n++)
      {
        // Increase the associated global equation number counter
        ++Count[elem_pt->eqn_number(n)];
      }
    }

    // Calculate the value Phi by
    // solving the system JPhi = dF/dlambda

    // Locally cache the linear solver
    LinearSolver* const linear_solver_pt = problem_pt->linear_solver_pt();

    // Save the status before entry to this routine
    bool enable_resolve = linear_solver_pt->is_resolve_enabled();

    // We need to do a resolve
    linear_solver_pt->enable_resolve();

    // Storage for the solution
    DoubleVector x(dist_pt, 0.0);

    // Solve the standard problem, we only want to make sure that
    // we factorise the matrix, if it has not been factorised. We shall
    // ignore the return value of x.
    linear_solver_pt->solve(problem_pt, x);

    // Get the vector dresiduals/dparameter
    problem_pt->get_derivative_wrt_global_parameter(parameter_pt, x);

    // Copy rhs vector into local storage so it doesn't get overwritten
    // if the linear solver decides to initialise the solution vector, say,
    // which it's quite entitled to do!
    DoubleVector input_x(x);

    // Now resolve the system with the new RHS and overwrite the solution
    linear_solver_pt->resolve(input_x, x);

    // Restore the storage status of the linear solver
    if (enable_resolve)
    {
      linear_solver_pt->enable_resolve();
    }
    else
    {
      linear_solver_pt->disable_resolve();
    }

    // Normalise the solution x
    double length = 0.0;
    for (unsigned n = 0; n < Ndof; n++)
    {
      length += x[n] * x[n];
    }
    length = sqrt(length);

    // Now add the real part of the null space components to the problem
    // unknowns and initialise it all
    // This is dumb at the moment ... fix with eigensolver?
    for (unsigned n = 0; n < Ndof; n++)
    {
      problem_pt->Dof_pt.push_back(&Phi[n]);
      C[n] = Phi[n] = -x[n] / length;
    }

    // Set the imaginary part so that the appropriate residual is
    // zero initially (eigensolvers)
    for (unsigned n = 0; n < Ndof; n += 2)
    {
      // Make sure that we are not at the end of an array of odd length
      if (n != Ndof - 1)
      {
        Psi[n] = C[n + 1];
        Psi[n + 1] = -C[n];
      }
      // If it's odd set the final entry to zero
      else
      {
        Psi[n] = 0.0;
      }
    }

    // Next add the imaginary parts of the null space components to the problem
    for (unsigned n = 0; n < Ndof; n++)
    {
      problem_pt->Dof_pt.push_back(&Psi[n]);
    }
    // Now add the parameter
    problem_pt->Dof_pt.push_back(parameter_pt);
    // Finally add the frequency
    problem_pt->Dof_pt.push_back(&Omega);

    // rebuild the dof dist
    Problem_pt->Dof_distribution_pt->build(
      Problem_pt->communicator_pt(), Ndof * 3 + 2, false);
    // Remove all previous sparse storage used during Jacobian assembly
    Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

    // delete the dist_pt
    delete dist_pt;
  }

  //====================================================================
  /// Constructor: Initialise the hopf handler,
  /// by setting initial guesses for Phi, Psi, Omega  and calculating Count.
  /// If the system changes, a new  handler must be constructed.
  //===================================================================
  HopfHandler::HopfHandler(Problem* const& problem_pt,
                           double* const& parameter_pt,
                           const double& omega,
                           const DoubleVector& phi,
                           const DoubleVector& psi)
    : Solve_which_system(0), Parameter_pt(parameter_pt), Omega(omega)
  {
    // Set the problem pointer
    Problem_pt = problem_pt;
    // Set the number of non-augmented degrees of freedom
    Ndof = problem_pt->ndof();

    // Resize the vectors of additional dofs
    Phi.resize(Ndof);
    Psi.resize(Ndof);
    C.resize(Ndof);
    Count.resize(Ndof, 0);

    // Loop over all the elements in the problem
    unsigned n_element = problem_pt->mesh_pt()->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the local freedoms in an element
      unsigned n_var = elem_pt->ndof();
      for (unsigned n = 0; n < n_var; n++)
      {
        // Increase the associated global equation number counter
        ++Count[elem_pt->eqn_number(n)];
      }
    }

    // Normalise the guess for phi
    double length = 0.0;
    for (unsigned n = 0; n < Ndof; n++)
    {
      length += phi[n] * phi[n];
    }
    length = sqrt(length);

    // Now add the real part of the null space components to the problem
    // unknowns and initialise it all
    for (unsigned n = 0; n < Ndof; n++)
    {
      problem_pt->Dof_pt.push_back(&Phi[n]);
      C[n] = Phi[n] = phi[n] / length;
      Psi[n] = psi[n] / length;
    }

    // Next add the imaginary parts of the null space components to the problem
    for (unsigned n = 0; n < Ndof; n++)
    {
      problem_pt->Dof_pt.push_back(&Psi[n]);
    }

    // Now add the parameter
    problem_pt->Dof_pt.push_back(parameter_pt);
    // Finally add the frequency
    problem_pt->Dof_pt.push_back(&Omega);

    // rebuild the Dof_distribution_pt
    Problem_pt->Dof_distribution_pt->build(
      Problem_pt->communicator_pt(), Ndof * 3 + 2, false);
    // Remove all previous sparse storage used during Jacobian assembly
    Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);
  }


  //=======================================================================
  /// Destructor return the problem to its original (non-augmented) state
  //=======================================================================
  HopfHandler::~HopfHandler()
  {
    // If we are using the block solver reset the problem's linear solver
    // to the original one
    BlockHopfLinearSolver* block_hopf_solver_pt =
      dynamic_cast<BlockHopfLinearSolver*>(Problem_pt->linear_solver_pt());
    if (block_hopf_solver_pt)
    {
      // Reset the problem's linear solver
      Problem_pt->linear_solver_pt() = block_hopf_solver_pt->linear_solver_pt();
      // Delete the block solver
      delete block_hopf_solver_pt;
    }
    // Now return the problem to its original size
    Problem_pt->Dof_pt.resize(Ndof);
    Problem_pt->Dof_distribution_pt->build(
      Problem_pt->communicator_pt(), Ndof, false);
    // Remove all previous sparse storage used during Jacobian assembly
    Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);
  }


  //=============================================================
  /// Get the number of elemental degrees of freedom
  //=============================================================
  unsigned HopfHandler::ndof(GeneralisedElement* const& elem_pt)
  {
    unsigned raw_ndof = elem_pt->ndof();
    switch (Solve_which_system)
    {
        // Full augmented system
      case 0:
        return (3 * raw_ndof + 2);
        break;
        // Standard non-augmented system
      case 1:
        return raw_ndof;
        break;
        // Complex system
      case 2:
        return (2 * raw_ndof);
        break;

      default:
        throw OomphLibError("Solve_which_system can only be 0,1 or 2",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
  }

  //=============================================================
  /// Get the global equation number of the local unknown
  //============================================================
  unsigned long HopfHandler::eqn_number(GeneralisedElement* const& elem_pt,
                                        const unsigned& ieqn_local)
  {
    // Get the raw value
    unsigned raw_ndof = elem_pt->ndof();
    unsigned long global_eqn;
    if (ieqn_local < raw_ndof)
    {
      global_eqn = elem_pt->eqn_number(ieqn_local);
    }
    else if (ieqn_local < 2 * raw_ndof)
    {
      global_eqn = Ndof + elem_pt->eqn_number(ieqn_local - raw_ndof);
    }
    else if (ieqn_local < 3 * raw_ndof)
    {
      global_eqn = 2 * Ndof + elem_pt->eqn_number(ieqn_local - 2 * raw_ndof);
    }
    else if (ieqn_local == 3 * raw_ndof)
    {
      global_eqn = 3 * Ndof;
    }
    else
    {
      global_eqn = 3 * Ndof + 1;
    }
    return global_eqn;
  }

  //==================================================================
  /// Get the residuals
  //=================================================================
  void HopfHandler::get_residuals(GeneralisedElement* const& elem_pt,
                                  Vector<double>& residuals)
  {
    // Should only call get residuals for the full system
    if (Solve_which_system == 0)
    {
      // Need to get raw residuals and jacobian
      unsigned raw_ndof = elem_pt->ndof();

      DenseMatrix<double> jacobian(raw_ndof), M(raw_ndof);
      // Get the basic residuals, jacobian and mass matrix
      elem_pt->get_jacobian_and_mass_matrix(residuals, jacobian, M);

      // Initialise the pen-ultimate residual
      residuals[3 * raw_ndof] =
        -1.0 / (double)(Problem_pt->mesh_pt()->nelement());
      residuals[3 * raw_ndof + 1] = 0.0;

      // Now multiply to fill in the residuals
      for (unsigned i = 0; i < raw_ndof; i++)
      {
        residuals[raw_ndof + i] = 0.0;
        residuals[2 * raw_ndof + i] = 0.0;
        for (unsigned j = 0; j < raw_ndof; j++)
        {
          unsigned global_unknown = elem_pt->eqn_number(j);
          // Real part
          residuals[raw_ndof + i] += jacobian(i, j) * Phi[global_unknown] +
                                     Omega * M(i, j) * Psi[global_unknown];
          // Imaginary part
          residuals[2 * raw_ndof + i] += jacobian(i, j) * Psi[global_unknown] -
                                         Omega * M(i, j) * Phi[global_unknown];
        }
        // Get the global equation number
        unsigned global_eqn = elem_pt->eqn_number(i);

        // Real part
        residuals[3 * raw_ndof] +=
          (Phi[global_eqn] * C[global_eqn]) / Count[global_eqn];
        // Imaginary part
        residuals[3 * raw_ndof + 1] +=
          (Psi[global_eqn] * C[global_eqn]) / Count[global_eqn];
      }
    }
    else
    {
      throw OomphLibError("Solve_which_system can only be 0",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  }


  //===============================================================
  /// Calculate the elemental Jacobian matrix "d equation
  /// / d variable".
  //==================================================================
  void HopfHandler::get_jacobian(GeneralisedElement* const& elem_pt,
                                 Vector<double>& residuals,
                                 DenseMatrix<double>& jacobian)
  {
    // The standard case
    if (Solve_which_system == 0)
    {
      unsigned augmented_ndof = ndof(elem_pt);
      unsigned raw_ndof = elem_pt->ndof();

      // Get the basic residuals and jacobian
      DenseMatrix<double> M(raw_ndof);
      elem_pt->get_jacobian_and_mass_matrix(residuals, jacobian, M);
      // Now fill in the actual residuals
      get_residuals(elem_pt, residuals);

      // Now the jacobian appears in other entries
      for (unsigned n = 0; n < raw_ndof; ++n)
      {
        for (unsigned m = 0; m < raw_ndof; ++m)
        {
          jacobian(raw_ndof + n, raw_ndof + m) = jacobian(n, m);
          jacobian(raw_ndof + n, 2 * raw_ndof + m) = Omega * M(n, m);
          jacobian(2 * raw_ndof + n, 2 * raw_ndof + m) = jacobian(n, m);
          jacobian(2 * raw_ndof + n, raw_ndof + m) = -Omega * M(n, m);
          unsigned global_eqn = elem_pt->eqn_number(m);
          jacobian(raw_ndof + n, 3 * raw_ndof + 1) += M(n, m) * Psi[global_eqn];
          jacobian(2 * raw_ndof + n, 3 * raw_ndof + 1) -=
            M(n, m) * Phi[global_eqn];
        }

        unsigned local_eqn = elem_pt->eqn_number(n);
        jacobian(3 * raw_ndof, raw_ndof + n) = C[local_eqn] / Count[local_eqn];
        jacobian(3 * raw_ndof + 1, 2 * raw_ndof + n) =
          C[local_eqn] / Count[local_eqn];
      }

      const double FD_step = 1.0e-8;

      Vector<double> newres_p(augmented_ndof), newres_m(augmented_ndof);

      // Loop over the dofs
      for (unsigned n = 0; n < raw_ndof; n++)
      {
        // Just do the x's
        unsigned long global_eqn = eqn_number(elem_pt, n);
        double* unknown_pt = Problem_pt->Dof_pt[global_eqn];
        double init = *unknown_pt;
        *unknown_pt += FD_step;

        // Get the new residuals
        get_residuals(elem_pt, newres_p);

        // Reset
        *unknown_pt = init;

        // Subtract
        *unknown_pt -= FD_step;
        get_residuals(elem_pt, newres_m);

        for (unsigned m = 0; m < raw_ndof; m++)
        {
          jacobian(raw_ndof + m, n) =
            (newres_p[raw_ndof + m] - residuals[raw_ndof + m]) / (FD_step);
          jacobian(2 * raw_ndof + m, n) =
            (newres_p[2 * raw_ndof + m] - residuals[2 * raw_ndof + m]) /
            (FD_step);
        }
        // Reset the unknown
        *unknown_pt = init;
      }

      {
        // Now do the global parameter
        double* unknown_pt = Problem_pt->Dof_pt[3 * Ndof];
        double init = *unknown_pt;
        *unknown_pt += FD_step;

        Problem_pt->actions_after_change_in_bifurcation_parameter();
        // Get the new residuals
        get_residuals(elem_pt, newres_p);

        // Reset
        *unknown_pt = init;

        // Subtract
        *unknown_pt -= FD_step;
        get_residuals(elem_pt, newres_m);

        for (unsigned m = 0; m < augmented_ndof - 2; m++)
        {
          jacobian(m, 3 * raw_ndof) = (newres_p[m] - residuals[m]) / FD_step;
        }
        // Reset the unknown
        *unknown_pt = init;
        Problem_pt->actions_after_change_in_bifurcation_parameter();
      }
    } // End of standard case
    // Normal case
    else if (Solve_which_system == 1)
    {
      // Just get the normal jacobian and residuals
      elem_pt->get_jacobian(residuals, jacobian);
    }
    // Otherwise the augmented complex case
    else if (Solve_which_system == 2)
    {
      unsigned raw_ndof = elem_pt->ndof();

      // Get the basic residuals and jacobian
      DenseMatrix<double> M(raw_ndof);
      elem_pt->get_jacobian_and_mass_matrix(residuals, jacobian, M);

      // We now need to fill in the other blocks
      for (unsigned n = 0; n < raw_ndof; n++)
      {
        for (unsigned m = 0; m < raw_ndof; m++)
        {
          jacobian(n, raw_ndof + m) = Omega * M(n, m);
          jacobian(raw_ndof + n, m) = -Omega * M(n, m);
          jacobian(raw_ndof + n, raw_ndof + m) = jacobian(n, m);
        }
      }

      // Now overwrite to fill in the residuals
      // The decision take is to solve for the mass matrix multiplied
      // terms in the residuals because they require no additional
      // information to assemble.
      for (unsigned n = 0; n < raw_ndof; n++)
      {
        residuals[n] = 0.0;
        residuals[raw_ndof + n] = 0.0;
        for (unsigned m = 0; m < raw_ndof; m++)
        {
          unsigned global_unknown = elem_pt->eqn_number(m);
          // Real part
          residuals[n] += M(n, m) * Psi[global_unknown];
          // Imaginary part
          residuals[raw_ndof + n] -= M(n, m) * Phi[global_unknown];
        }
      }
    } // End of complex augmented case
    else
    {
      throw OomphLibError("Solve_which_system can only be 0,1 or 2",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  }


  //==================================================================
  /// Get the derivatives of the augmented residuals with respect to
  /// a parameter
  //=================================================================
  void HopfHandler::get_dresiduals_dparameter(
    GeneralisedElement* const& elem_pt,
    double* const& parameter_pt,
    Vector<double>& dres_dparam)
  {
    // Should only call get residuals for the full system
    if (Solve_which_system == 0)
    {
      // Need to get raw residuals and jacobian
      unsigned raw_ndof = elem_pt->ndof();

      DenseMatrix<double> djac_dparam(raw_ndof), dM_dparam(raw_ndof);
      // Get the basic residuals, jacobian and mass matrix
      elem_pt->get_djacobian_and_dmass_matrix_dparameter(
        parameter_pt, dres_dparam, djac_dparam, dM_dparam);

      // Initialise the pen-ultimate residual, which does not
      // depend on the parameter
      dres_dparam[3 * raw_ndof] = 0.0;
      dres_dparam[3 * raw_ndof + 1] = 0.0;

      // Now multiply to fill in the residuals
      for (unsigned i = 0; i < raw_ndof; i++)
      {
        dres_dparam[raw_ndof + i] = 0.0;
        dres_dparam[2 * raw_ndof + i] = 0.0;
        for (unsigned j = 0; j < raw_ndof; j++)
        {
          unsigned global_unknown = elem_pt->eqn_number(j);
          // Real part
          dres_dparam[raw_ndof + i] +=
            djac_dparam(i, j) * Phi[global_unknown] +
            Omega * dM_dparam(i, j) * Psi[global_unknown];
          // Imaginary part
          dres_dparam[2 * raw_ndof + i] +=
            djac_dparam(i, j) * Psi[global_unknown] -
            Omega * dM_dparam(i, j) * Phi[global_unknown];
        }
      }
    }
    else
    {
      throw OomphLibError("Solve_which_system can only be 0",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  }


  //========================================================================
  /// Overload the derivative of the residuals and jacobian
  /// with respect to a parameter so that it breaks because it should not
  /// be required
  //========================================================================
  void HopfHandler::get_djacobian_dparameter(GeneralisedElement* const& elem_pt,
                                             double* const& parameter_pt,
                                             Vector<double>& dres_dparam,
                                             DenseMatrix<double>& djac_dparam)
  {
    std::ostringstream error_stream;
    error_stream
      << "This function has not been implemented because it is not required\n";
    error_stream << "in standard problems.\n";
    error_stream
      << "If you find that you need it, you will have to implement it!\n\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


  //=====================================================================
  /// Overload the hessian vector product function so that
  /// it breaks because it should not be required
  //========================================================================
  void HopfHandler::get_hessian_vector_products(
    GeneralisedElement* const& elem_pt,
    Vector<double> const& Y,
    DenseMatrix<double> const& C,
    DenseMatrix<double>& product)
  {
    std::ostringstream error_stream;
    error_stream
      << "This function has not been implemented because it is not required\n";
    error_stream << "in standard problems.\n";
    error_stream
      << "If you find that you need it, you will have to implement it!\n\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


  //==========================================================================
  /// Return the eigenfunction(s) associated with the bifurcation that
  /// has been detected in bifurcation tracking problems
  //==========================================================================
  void HopfHandler::get_eigenfunction(Vector<DoubleVector>& eigenfunction)
  {
    // There is a real and imaginary part of the null vector
    eigenfunction.resize(2);
    LinearAlgebraDistribution dist(Problem_pt->communicator_pt(), Ndof, false);
    // Rebuild the vector
    eigenfunction[0].build(&dist, 0.0);
    eigenfunction[1].build(&dist, 0.0);
    // Set the value to be the null vector
    for (unsigned n = 0; n < Ndof; n++)
    {
      eigenfunction[0][n] = Phi[n];
      eigenfunction[1][n] = Psi[n];
    }
  }


  //====================================================================
  /// Set to solve the standard (underlying jacobian)  system
  //===================================================================
  void HopfHandler::solve_standard_system()
  {
    if (Solve_which_system != 1)
    {
      Solve_which_system = 1;
      // Restrict the problem to the standard variables only
      Problem_pt->Dof_pt.resize(Ndof);
      Problem_pt->Dof_distribution_pt->build(
        Problem_pt->communicator_pt(), Ndof, false);
      // Remove all previous sparse storage used during Jacobian assembly
      Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);
    }
  }

  //====================================================================
  /// Set to solve the complex (jacobian and mass matrix)  system
  //===================================================================
  void HopfHandler::solve_complex_system()
  {
    // If we were not solving the complex system resize the unknowns
    // accordingly
    if (Solve_which_system != 2)
    {
      Solve_which_system = 2;
      // Resize to the first Ndofs (will work whichever system we were
      // solving before)
      Problem_pt->Dof_pt.resize(Ndof);
      // Add the first (real) part of the eigenfunction back into the problem
      for (unsigned n = 0; n < Ndof; n++)
      {
        Problem_pt->Dof_pt.push_back(&Phi[n]);
      }
      Problem_pt->Dof_distribution_pt->build(
        Problem_pt->communicator_pt(), Ndof * 2, false);
      // Remove all previous sparse storage used during Jacobian assembly
      Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);
    }
  }


  //=================================================================
  /// Set to Solve full system system
  //=================================================================
  void HopfHandler::solve_full_system()
  {
    // If we are starting from another system
    if (Solve_which_system)
    {
      Solve_which_system = 0;

      // Resize to the first Ndofs (will work whichever system we were
      // solving before)
      Problem_pt->Dof_pt.resize(Ndof);
      // Add the additional unknowns back into the problem
      for (unsigned n = 0; n < Ndof; n++)
      {
        Problem_pt->Dof_pt.push_back(&Phi[n]);
      }
      for (unsigned n = 0; n < Ndof; n++)
      {
        Problem_pt->Dof_pt.push_back(&Psi[n]);
      }
      // Now add the parameter
      Problem_pt->Dof_pt.push_back(Parameter_pt);
      // Finally add the frequency
      Problem_pt->Dof_pt.push_back(&Omega);

      //
      Problem_pt->Dof_distribution_pt->build(
        Problem_pt->communicator_pt(), 3 * Ndof + 2, false);
      // Remove all previous sparse storage used during Jacobian assembly
      Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);
    }
  }
} // namespace oomph
