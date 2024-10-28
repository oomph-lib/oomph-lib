#define BOOST_TEST_MODULE axisym_cap_test_module
#include <boost/test/included/unit_test.hpp>

#include "run_tests.h"

using namespace std;
using namespace oomph;

// ***
// Test matrix comparison

BOOST_AUTO_TEST_CASE(compare_matrix_same)
{
  const unsigned N = 4;
  DenseDoubleMatrix A(N, N, 0.0);
  DenseDoubleMatrix B(N, N, 0.0);
  A(1, 1) = 1.0;
  B(1, 1) = 1.0;
  B(2, 1) = 1e-10;
  BOOST_TEST(compare_matrices(A, B) == 1);
}

BOOST_AUTO_TEST_CASE(compare_matrix_different)
{
  const unsigned N = 4;
  DenseDoubleMatrix A(N, N, 0.0);
  DenseDoubleMatrix B(N, N, 0.0);
  A(2, 1) = 1.0;
  B(2, 1) = 1e-10;
  BOOST_TEST(compare_matrices(A, B) == 0);
}

// ***
// Test base problem

BOOST_AUTO_TEST_CASE(nonlinear_problem_creation_with_default_parameters)
{
  createBaseProblem();
}

BOOST_AUTO_TEST_CASE(default_problem_solves)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  problem_pt->steady_newton_solve();
}

BOOST_AUTO_TEST_CASE(default_problem_solves_with_adapt)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  problem_pt->steady_newton_solve(1);
}

BOOST_AUTO_TEST_CASE(nonlinear_analytic_J_equals_fd_J)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  problem_pt->steady_newton_solve();
  DoubleVector residuals;
  DenseMatrix<double> expected_jacobian;
  DenseDoubleMatrix actual_jacobian;
  problem_pt->get_fd_jacobian(residuals, expected_jacobian);
  problem_pt->get_jacobian(residuals, actual_jacobian);

  const unsigned N = expected_jacobian.nrow();
  const unsigned M = expected_jacobian.ncol();

  BOOST_TEST(actual_jacobian.nrow() == N);
  BOOST_TEST(actual_jacobian.ncol() == M);

  bool jacobians_are_equal =
    compare_matrices(expected_jacobian, actual_jacobian);
  BOOST_TEST(jacobians_are_equal);
}

BOOST_AUTO_TEST_CASE(nonlinear_analytic_M_equals_fd_M)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  problem_pt->steady_newton_solve();
  DoubleVector residuals;
  CRDoubleMatrix actual_jacobian;
  CRDoubleMatrix actual_mass_matrix;
  problem_pt->get_eigenproblem_matrices(actual_mass_matrix, actual_jacobian);

  problem_pt->time_stepper_pt()->make_steady();
  DenseMatrix<double> expected_jacobian;
  problem_pt->get_fd_jacobian(residuals, expected_jacobian);
  problem_pt->time_stepper_pt()->undo_make_steady();
  DenseMatrix<double> unsteady_jacobian;
  problem_pt->get_fd_jacobian(residuals, unsteady_jacobian);
  const double dt = problem_pt->time_stepper_pt()->time_pt()->dt();
  const unsigned n_dof = problem_pt->ndof();
  DenseMatrix<double> expected_mass_matrix(n_dof, n_dof, 0.0);
  for (unsigned i = 0; i < n_dof; i++)
  {
    for (unsigned j = 0; j < n_dof; j++)
    {
      /// The 2/3 is due to the BDF<2> scheme.
      expected_mass_matrix(j, i) +=
        (2.0 / 3.0) * dt * (-unsteady_jacobian(j, i) + expected_jacobian(j, i));
    }
  }

  const unsigned N = expected_jacobian.nrow();
  const unsigned M = expected_jacobian.ncol();

  BOOST_TEST(actual_jacobian.nrow() == N);
  BOOST_TEST(actual_jacobian.ncol() == M);

  bool jacobians_are_equal =
    compare_matrices(expected_jacobian, actual_jacobian);
  BOOST_TEST(jacobians_are_equal);
  bool mass_matrix_are_equal =
    compare_matrices(expected_mass_matrix, actual_mass_matrix);
  BOOST_TEST(mass_matrix_are_equal);
}

BOOST_AUTO_TEST_CASE(nonlinear_problem_solve_eigenproblem)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  problem_pt->steady_newton_solve();
  problem_pt->reset_lagrange();
  problem_pt->assign_initial_values_impulsive();
  Vector<std::complex<double>> eigenvalue =
    problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  printf("%16.11g \n", eigenvalue[0].real());
  BOOST_TEST(abs(eigenvalue[0].real() - (-0.24783965941010913)) < 1e-6);
}

BOOST_AUTO_TEST_CASE(nonlinear_problem_solve_eigenproblem_with_fluid_flow)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  problem_pt->steady_newton_solve();
  problem_pt->reset_lagrange();
  problem_pt->assign_initial_values_impulsive();
  Vector<std::complex<double>> eigenvalue =
    problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  printf("%16.11g \n", eigenvalue[0].real());
  BOOST_TEST(abs(eigenvalue[0].real() - (-0.27630985141964381)) < 1e-6);
}

// Unpin the horizontal dofs
// We don't expect the these to match. The inner boundary conditions imply
// dR/ds = 0 but this is not imposed.
// BOOST_AUTO_TEST_CASE(
//  nonlinear_problem_solve_eigenproblem_with_and_without_x_deformation)
//{
//  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//  problem_pt->steady_newton_solve();
//  problem_pt->reset_lagrange();
//  problem_pt->assign_initial_values_impulsive();
//  Vector<std::complex<double>> expected_eigenvalue =
//    problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
//  printf("%16.11g \n", expected_eigenvalue[0].real());
//
//  problem_pt->pin_horizontal_displacement();
//
//  Vector<std::complex<double>> actual_eigenvalue =
//    problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
//  printf("%16.11g \n", actual_eigenvalue[0].real());
//
//  BOOST_TEST(abs(expected_eigenvalue[0].real() - actual_eigenvalue[0].real())
//  <
//             1e-3);
//}

BOOST_AUTO_TEST_CASE(nonlinear_problem_unsteady_run)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Flux_Parameters::Withdraw_speed = 0.1;
  const double dt = 0.005;
  const double ft = 3 * dt;
  problem_pt->timestep(dt, ft);

  const double actual_centre_height = problem_pt->get_centre_point_z();
  BOOST_TEST(abs(actual_centre_height - 0.0028304837510845986) < 1e-5);
}

// ***
// Test linear problem

BOOST_AUTO_TEST_CASE(linear_problem_creation_from_base_problem)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  createLinearProblem(problem_pt);
}

BOOST_AUTO_TEST_CASE(linear_analytic_J_equals_fd_J)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  problem_pt->steady_newton_solve();
  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt);
  DoubleVector residuals;
  DenseMatrix<double> expected_jacobian;
  DenseDoubleMatrix actual_jacobian;
  linear_problem_pt->get_jacobian(residuals, actual_jacobian);
  linear_problem_pt->get_fd_jacobian(residuals, expected_jacobian);

  const unsigned N = expected_jacobian.nrow();
  const unsigned M = expected_jacobian.ncol();

  BOOST_TEST(actual_jacobian.nrow() == N);
  BOOST_TEST(actual_jacobian.ncol() == M);

  bool jacobians_are_equal =
    compare_matrices(expected_jacobian, actual_jacobian);
  BOOST_TEST(jacobians_are_equal);
}

BOOST_AUTO_TEST_CASE(linear_analytic_M_equals_fd_M)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  problem_pt->steady_newton_solve();
  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt);
  DoubleVector residuals;
  CRDoubleMatrix actual_jacobian;
  CRDoubleMatrix actual_mass_matrix;
  linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
                                               actual_jacobian);

  linear_problem_pt->time_stepper_pt()->make_steady();
  DenseMatrix<double> expected_jacobian;
  linear_problem_pt->get_fd_jacobian(residuals, expected_jacobian);
  linear_problem_pt->time_stepper_pt()->undo_make_steady();
  DenseMatrix<double> unsteady_jacobian;
  linear_problem_pt->get_fd_jacobian(residuals, unsteady_jacobian);
  const double dt = linear_problem_pt->time_stepper_pt()->time_pt()->dt();
  const unsigned n_dof = linear_problem_pt->ndof();
  DenseMatrix<double> expected_mass_matrix(n_dof, n_dof, 0.0);
  for (unsigned i = 0; i < n_dof; i++)
  {
    for (unsigned j = 0; j < n_dof; j++)
    {
      /// The 2/3 is due to the BDF<2> scheme.
      expected_mass_matrix(j, i) +=
        (2.0 / 3.0) * dt * (-unsteady_jacobian(j, i) + expected_jacobian(j, i));
    }
  }

  const unsigned N = expected_jacobian.nrow();
  const unsigned M = expected_jacobian.ncol();

  BOOST_TEST(actual_jacobian.nrow() == N);
  BOOST_TEST(actual_jacobian.ncol() == M);

  bool jacobians_are_equal =
    compare_matrices(expected_jacobian, actual_jacobian);
  BOOST_TEST(jacobians_are_equal);
  bool mass_matrix_are_equal =
    compare_matrices(expected_mass_matrix, actual_mass_matrix);
  BOOST_TEST(mass_matrix_are_equal);
}

BOOST_AUTO_TEST_CASE(linear_analytic_M_equals_fd_M_mode_1)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  problem_pt->steady_newton_solve();
  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt, 1);
  DoubleVector residuals;
  CRDoubleMatrix actual_jacobian;
  CRDoubleMatrix actual_mass_matrix;
  linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
                                               actual_jacobian);

  linear_problem_pt->time_stepper_pt()->make_steady();
  DenseMatrix<double> expected_jacobian;
  linear_problem_pt->get_fd_jacobian(residuals, expected_jacobian);
  linear_problem_pt->time_stepper_pt()->undo_make_steady();
  DenseMatrix<double> unsteady_jacobian;
  linear_problem_pt->get_fd_jacobian(residuals, unsteady_jacobian);
  const double dt = linear_problem_pt->time_stepper_pt()->time_pt()->dt();
  const unsigned n_dof = linear_problem_pt->ndof();
  DenseMatrix<double> expected_mass_matrix(n_dof, n_dof, 0.0);
  for (unsigned i = 0; i < n_dof; i++)
  {
    for (unsigned j = 0; j < n_dof; j++)
    {
      /// The 2/3 is due to the BDF<2> scheme.
      expected_mass_matrix(j, i) +=
        (2.0 / 3.0) * dt * (-unsteady_jacobian(j, i) + expected_jacobian(j, i));
    }
  }

  const unsigned N = expected_jacobian.nrow();
  const unsigned M = expected_jacobian.ncol();

  BOOST_TEST(actual_jacobian.nrow() == N);
  BOOST_TEST(actual_jacobian.ncol() == M);
  save_dofs_types(linear_problem_pt, "dofs.txt");

  bool jacobians_are_equal =
    compare_matrices(expected_jacobian, actual_jacobian);
  BOOST_TEST(jacobians_are_equal);
  bool mass_matrix_are_equal =
    compare_matrices(expected_mass_matrix, actual_mass_matrix);
  BOOST_TEST(mass_matrix_are_equal);
}


BOOST_AUTO_TEST_CASE(linear_problem_solve_eigenproblem)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  problem_pt->steady_newton_solve();
  problem_pt->reset_lagrange();
  problem_pt->assign_initial_values_impulsive();
  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt);

  linear_problem_pt->assign_initial_values_impulsive();
  Vector<std::complex<double>> eigenvalue =
    linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  printf("%16.11g \n", eigenvalue[0].real());

  BOOST_TEST(abs(eigenvalue[0].real() - (-0.24783965941010913)) < 1e-6);
}

BOOST_AUTO_TEST_CASE(
  linear_problem_solve_eigenproblem_with_and_without_x_deformation)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  problem_pt->steady_newton_solve();
  problem_pt->reset_lagrange();
  problem_pt->assign_initial_values_impulsive();
  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt);
  linear_problem_pt->assign_initial_values_impulsive();
  Vector<std::complex<double>> expected_eigenvalue =
    linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  printf("%16.11g \n", expected_eigenvalue[0].real());

  linear_problem_pt->pin_horizontal_mesh_deformation();

  Vector<std::complex<double>> actual_eigenvalue =
    linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  printf("%16.11g \n", actual_eigenvalue[0].real());

  BOOST_TEST(abs(expected_eigenvalue[0].real() - actual_eigenvalue[0].real()) <
             1e-3);
}

BOOST_AUTO_TEST_CASE(linear_problem_solve_eigenproblem_with_fluid_flow)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  problem_pt->steady_newton_solve();
  problem_pt->reset_lagrange();
  problem_pt->assign_initial_values_impulsive();

  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt);
  linear_problem_pt->assign_initial_values_impulsive();
  Vector<std::complex<double>> eigenvalue =
    linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  printf("%16.11g \n", eigenvalue[0].real());

  BOOST_TEST(abs(eigenvalue[0].real() - (-0.27630985141964381)) < 1e-6);
}

BOOST_AUTO_TEST_CASE(mode_1_problem)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  problem_pt->steady_newton_solve();
  problem_pt->reset_lagrange();
  problem_pt->assign_initial_values_impulsive();

  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt, 1);
  linear_problem_pt->assign_initial_values_impulsive();

  Vector<std::complex<double>> eigenvalue =
    linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  printf("%16.11g \n", eigenvalue[0].real());

  linear_problem_pt->doc_solution();

  BOOST_TEST(abs(eigenvalue[0].real() - (-0.054927275769394331)) < 1e-6);
}

// We don't expect this test to pass. Missing a constraint at the centre node.
// BOOST_AUTO_TEST_CASE(
//   linear_problem_solve_eigenproblem_with_and_without_x_deformation_mode_1)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->steady_newton_solve();
//   problem_pt->reset_lagrange();
//   problem_pt->assign_initial_values_impulsive();
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt, 1);
//   linear_problem_pt->assign_initial_values_impulsive();
//   Vector<std::complex<double>> expected_eigenvalue =
//     linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
//   printf("%16.11g \n", expected_eigenvalue[0].real());
//
//   linear_problem_pt->pin_horizontal_mesh_deformation();
//
//   Vector<std::complex<double>> actual_eigenvalue =
//     linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
//   printf("%16.11g \n", actual_eigenvalue[0].real());
//
//   BOOST_TEST(abs(expected_eigenvalue[0].real() - actual_eigenvalue[0].real())
//   <
//              1e-3);
// }

BOOST_AUTO_TEST_CASE(mode_1_problem_velocity_gradient_at_the_centre)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  problem_pt->steady_newton_solve(1);
  problem_pt->reset_lagrange();
  problem_pt->assign_initial_values_impulsive();

  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt, 1);
  linear_problem_pt->assign_initial_values_impulsive();

  linear_problem_pt->make_unsteady();
  linear_problem_pt->pin_horizontal_mesh_deformation();

  linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(8);

  linear_problem_pt->perturb_with_eigensolution(0.01);

  linear_problem_pt->doc_solution();

  Vector<double> gradu = linear_problem_pt->estimate_gradient_at_origin();

  BOOST_TEST(abs(gradu[0]) < 1e-3);
  BOOST_TEST(abs(gradu[1]) < 1e-3);
  BOOST_TEST(abs(gradu[4]) < 1e-3);
  BOOST_TEST(abs(gradu[5]) < 1e-3);
}

BOOST_AUTO_TEST_CASE(mode_1_problem_velocity_gradient_at_the_centre_flow)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  problem_pt->steady_newton_solve(1);
  problem_pt->reset_lagrange();
  problem_pt->assign_initial_values_impulsive();

  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt, 1);
  linear_problem_pt->assign_initial_values_impulsive();

  linear_problem_pt->make_unsteady();
  linear_problem_pt->pin_horizontal_mesh_deformation();

  linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(8);

  linear_problem_pt->perturb_with_eigensolution(0.01);

  linear_problem_pt->doc_solution();

  Vector<double> gradu = linear_problem_pt->estimate_gradient_at_origin();

  BOOST_TEST(abs(gradu[0]) < 1e-3);
  BOOST_TEST(abs(gradu[1]) < 1e-3);
  BOOST_TEST(abs(gradu[4]) < 1e-3);
  BOOST_TEST(abs(gradu[5]) < 1e-3);
}

BOOST_AUTO_TEST_CASE(mode_2_problem)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  problem_pt->steady_newton_solve();
  problem_pt->reset_lagrange();
  problem_pt->assign_initial_values_impulsive();

  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt, 2);
  linear_problem_pt->assign_initial_values_impulsive();
  Vector<std::complex<double>> eigenvalue =
    linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  printf("%16.11g \n", eigenvalue[0].real());

  BOOST_TEST(abs(eigenvalue[0].real() - (-0.13394533183580556)) < 1e-6);
}

// This test doesn't pass. Needs more investigation. 10% error in eigenvalue
// BOOST_AUTO_TEST_CASE(
//   linear_problem_solve_eigenproblem_with_and_without_x_deformation_mode_2)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->steady_newton_solve();
//   problem_pt->reset_lagrange();
//   problem_pt->assign_initial_values_impulsive();
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt, 2);
//   linear_problem_pt->assign_initial_values_impulsive();
//   Vector<std::complex<double>> expected_eigenvalue =
//     linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
//   printf("%16.11g \n", expected_eigenvalue[0].real());
//
//   linear_problem_pt->pin_horizontal_mesh_deformation();
//
//   Vector<std::complex<double>> actual_eigenvalue =
//     linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
//   printf("%16.11g \n", actual_eigenvalue[0].real());
//
//   BOOST_TEST(abs(expected_eigenvalue[0].real() - actual_eigenvalue[0].real())
//   <
//              1e-4);
// }

BOOST_AUTO_TEST_CASE(linear_problem_unsteady_run)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  problem_pt->steady_newton_solve();
  problem_pt->reset_lagrange();
  problem_pt->assign_initial_values_impulsive();

  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt);
  Flux_Parameters::Withdraw_speed = 0.1;
  const double dt = 0.005;
  const double ft = 3 * dt;
  linear_problem_pt->timestep(dt, ft);

  // linear_problem_pt->doc_solution();
  const double actual_centre_height = linear_problem_pt->get_centre_point_ZC();
  BOOST_TEST(abs(actual_centre_height - 0.00045495963962505068) < 1e-6);
}


// // ***
// // Test axisym vs linear
//
// // --- TEST SOLID
//
// BOOST_AUTO_TEST_CASE(compare_jacobians_interior_solid_nodes)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->pin_fluid();
//   problem_pt->pin_solid_boundaries();
//   problem_pt->steady_newton_solve();
//   problem_pt->reset_lagrange();
//   DoubleVector residuals;
//   DenseDoubleMatrix expected_jacobian;
//   problem_pt->get_jacobian(residuals, expected_jacobian);
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   DenseDoubleMatrix actual_jacobian;
//   linear_problem_pt->pin_fluid();
//   linear_problem_pt->pin_solid_boundaries();
//   linear_problem_pt->get_jacobian(residuals, actual_jacobian);
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//
//   save_dofs_types(problem_pt, "expected_dofs.txt");
//   save_dofs_types(linear_problem_pt, "actual_dofs.txt");
//
//   BOOST_TEST(jacobians_are_equal);
// }
//
// BOOST_AUTO_TEST_CASE(compare_jacobians_interior_solid_nodes_and_inner_boundary)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->pin_fluid();
//   problem_pt->pin_solid_boundaries(lower);
//   problem_pt->pin_solid_boundaries(upper);
//   problem_pt->pin_solid_boundaries(outer);
//   problem_pt->steady_newton_solve();
//   problem_pt->reset_lagrange();
//   DoubleVector residuals;
//   DenseDoubleMatrix expected_jacobian;
//   problem_pt->get_jacobian(residuals, expected_jacobian);
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->pin_fluid();
//   linear_problem_pt->pin_solid_boundaries(lower);
//   linear_problem_pt->pin_solid_boundaries(upper);
//   linear_problem_pt->pin_solid_boundaries(outer);
//   DenseDoubleMatrix actual_jacobian;
//   linear_problem_pt->get_jacobian(residuals, actual_jacobian);
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//
//   save_dofs_types(problem_pt, "expected_dofs.txt");
//   save_dofs_types(linear_problem_pt, "actual_dofs.txt");
//
//   BOOST_TEST(jacobians_are_equal);
// }
//
// BOOST_AUTO_TEST_CASE(compare_jacobians_interior_solid_nodes_and_outer_boundary)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->pin_fluid();
//   problem_pt->pin_solid_boundaries(lower);
//   problem_pt->pin_solid_boundaries(inner);
//   problem_pt->pin_solid_boundaries(upper);
//   problem_pt->steady_newton_solve();
//   problem_pt->reset_lagrange();
//   DoubleVector residuals;
//   DenseDoubleMatrix expected_jacobian;
//   problem_pt->get_jacobian(residuals, expected_jacobian);
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->pin_fluid();
//   linear_problem_pt->pin_solid_boundaries(lower);
//   linear_problem_pt->pin_solid_boundaries(inner);
//   linear_problem_pt->pin_solid_boundaries(upper);
//   DenseDoubleMatrix actual_jacobian;
//   linear_problem_pt->get_jacobian(residuals, actual_jacobian);
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//
//   save_dofs_types(problem_pt, "expected_dofs.txt");
//   save_dofs_types(linear_problem_pt, "actual_dofs.txt");
//
//   BOOST_TEST(jacobians_are_equal);
// }
//
// BOOST_AUTO_TEST_CASE(compare_jacobians_interior_solid_nodes_and_upper_boundary)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->pin_fluid();
//   problem_pt->pin_solid_boundaries(lower);
//   problem_pt->pin_solid_boundaries(inner);
//   problem_pt->pin_solid_boundaries(outer);
//   problem_pt->steady_newton_solve();
//   problem_pt->reset_lagrange();
//   DoubleVector residuals;
//   DenseDoubleMatrix expected_jacobian;
//   problem_pt->get_jacobian(residuals, expected_jacobian);
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->pin_fluid();
//   linear_problem_pt->pin_solid_boundaries(lower);
//   linear_problem_pt->pin_solid_boundaries(inner);
//   linear_problem_pt->pin_solid_boundaries(outer);
//   DenseDoubleMatrix actual_jacobian;
//   linear_problem_pt->get_jacobian(residuals, actual_jacobian);
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//
//   save_dofs_types(problem_pt, "expected_dofs.txt");
//   save_dofs_types(linear_problem_pt, "actual_dofs.txt");
//
//   BOOST_TEST(jacobians_are_equal);
// }
//
// BOOST_AUTO_TEST_CASE(compare_jacobians_interior_solid_nodes_and_lower_boundary)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->pin_fluid();
//   problem_pt->pin_solid_boundaries(upper);
//   problem_pt->pin_solid_boundaries(inner);
//   problem_pt->pin_solid_boundaries(outer);
//   problem_pt->steady_newton_solve();
//   problem_pt->reset_lagrange();
//   DoubleVector residuals;
//   DenseDoubleMatrix expected_jacobian;
//   problem_pt->get_jacobian(residuals, expected_jacobian);
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->pin_fluid();
//   linear_problem_pt->pin_solid_boundaries(upper);
//   linear_problem_pt->pin_solid_boundaries(inner);
//   linear_problem_pt->pin_solid_boundaries(outer);
//   DenseDoubleMatrix actual_jacobian;
//   linear_problem_pt->get_jacobian(residuals, actual_jacobian);
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//
//   save_dofs_types(problem_pt, "expected_dofs.txt");
//   save_dofs_types(linear_problem_pt, "actual_dofs.txt");
//
//   BOOST_TEST(jacobians_are_equal);
// }
//
// // TEST FLUID JACOBIANS
// ------------------------------------------------------
//
// BOOST_AUTO_TEST_CASE(
//   compare_jacobians_and_mass_matrix_interior_fluid_nodes_no_flow)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->steady_newton_solve();
//   problem_pt->pin_solid();
//   problem_pt->pin_fluid_boundary(upper);
//   problem_pt->pin_fluid_boundary(lower);
//   problem_pt->pin_fluid_boundary(inner);
//   problem_pt->pin_fluid_boundary(outer);
//   CRDoubleMatrix expected_jacobian;
//   CRDoubleMatrix expected_mass_matrix;
//   problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
//                                         expected_jacobian);
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->pin_solid();
//   linear_problem_pt->pin_fluid_boundary(upper);
//   linear_problem_pt->pin_fluid_boundary(lower);
//   linear_problem_pt->pin_fluid_boundary(inner);
//   linear_problem_pt->pin_fluid_boundary(outer);
//   CRDoubleMatrix actual_jacobian;
//   CRDoubleMatrix actual_mass_matrix;
//   linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
//                                                actual_jacobian);
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//   BOOST_TEST(jacobians_are_equal);
//   bool mass_matrix_are_equal =
//     compare_matrices(expected_mass_matrix, actual_mass_matrix);
//   BOOST_TEST(mass_matrix_are_equal);
// }
//
// BOOST_AUTO_TEST_CASE(check_base_flow_passed_to_perturbed_problem_correctly)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->steady_newton_solve();
//   problem_pt->perturb_fluid();
//   problem_pt->pin_solid();
//   problem_pt->pin_fluid_boundary(upper);
//   problem_pt->pin_fluid_boundary(lower);
//   problem_pt->pin_fluid_boundary(inner);
//   problem_pt->pin_fluid_boundary(outer);
//
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->pin_solid();
//   linear_problem_pt->pin_fluid_boundary(upper);
//   linear_problem_pt->pin_fluid_boundary(lower);
//   linear_problem_pt->pin_fluid_boundary(inner);
//   linear_problem_pt->pin_fluid_boundary(outer);
// }
//
// BOOST_AUTO_TEST_CASE(
//   compare_jacobians_and_mass_matrix_interior_fluid_nodes_with_flow)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->steady_newton_solve();
//   problem_pt->perturb_fluid();
//   problem_pt->pin_solid();
//   problem_pt->pin_fluid_boundary(upper);
//   problem_pt->pin_fluid_boundary(lower);
//   problem_pt->pin_fluid_boundary(inner);
//   problem_pt->pin_fluid_boundary(outer);
//   CRDoubleMatrix expected_jacobian;
//   CRDoubleMatrix expected_mass_matrix;
//   problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
//                                         expected_jacobian);
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->pin_solid();
//   linear_problem_pt->pin_fluid_boundary(upper);
//   linear_problem_pt->pin_fluid_boundary(lower);
//   linear_problem_pt->pin_fluid_boundary(inner);
//   linear_problem_pt->pin_fluid_boundary(outer);
//   CRDoubleMatrix actual_jacobian;
//   CRDoubleMatrix actual_mass_matrix;
//   linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
//                                                actual_jacobian);
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//   BOOST_TEST(jacobians_are_equal);
//   bool mass_matrix_are_equal =
//     compare_matrices(expected_mass_matrix, actual_mass_matrix);
//   BOOST_TEST(mass_matrix_are_equal);
// }
//
//
// BOOST_AUTO_TEST_CASE(compare_jacobians_and_mass_matrix_interior_and_upper_fluid)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->steady_newton_solve();
//   problem_pt->pin_solid();
//   // problem_pt->pin_fluid_boundary(upper);
//   problem_pt->pin_fluid_boundary(lower);
//   problem_pt->pin_fluid_boundary(inner);
//   problem_pt->pin_fluid_boundary(outer);
//   CRDoubleMatrix expected_jacobian;
//   CRDoubleMatrix expected_mass_matrix;
//   problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
//                                         expected_jacobian);
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->pin_solid();
//   // linear_problem_pt->pin_fluid_boundary(upper);
//   linear_problem_pt->pin_fluid_boundary(lower);
//   linear_problem_pt->pin_fluid_boundary(inner);
//   linear_problem_pt->pin_fluid_boundary(outer);
//   CRDoubleMatrix actual_jacobian;
//   CRDoubleMatrix actual_mass_matrix;
//   linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
//                                                actual_jacobian);
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//   BOOST_TEST(jacobians_are_equal);
//   bool mass_matrix_are_equal =
//     compare_matrices(expected_mass_matrix, actual_mass_matrix);
//   BOOST_TEST(mass_matrix_are_equal);
// }
//
// BOOST_AUTO_TEST_CASE(compare_jacobians_and_mass_matrix_interior_and_outer_fluid)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->steady_newton_solve();
//   problem_pt->pin_solid();
//   problem_pt->pin_fluid_boundary(upper);
//   problem_pt->pin_fluid_boundary(lower);
//   problem_pt->pin_fluid_boundary(inner);
//   // problem_pt->pin_fluid_boundary(outer);
//   CRDoubleMatrix expected_jacobian;
//   CRDoubleMatrix expected_mass_matrix;
//   problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
//                                         expected_jacobian);
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->pin_solid();
//   linear_problem_pt->pin_fluid_boundary(upper);
//   linear_problem_pt->pin_fluid_boundary(lower);
//   linear_problem_pt->pin_fluid_boundary(inner);
//   // linear_problem_pt->pin_fluid_boundary(outer);
//   CRDoubleMatrix actual_jacobian;
//   CRDoubleMatrix actual_mass_matrix;
//   linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
//                                                actual_jacobian);
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//   BOOST_TEST(jacobians_are_equal);
//   bool mass_matrix_are_equal =
//     compare_matrices(expected_mass_matrix, actual_mass_matrix);
//   BOOST_TEST(mass_matrix_are_equal);
// }
//
// BOOST_AUTO_TEST_CASE(compare_jacobians_and_mass_matrix_interior_and_inner_fluid)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->steady_newton_solve();
//   problem_pt->pin_solid();
//   problem_pt->pin_fluid_boundary(upper);
//   problem_pt->pin_fluid_boundary(lower);
//   // problem_pt->pin_fluid_boundary(inner);
//   problem_pt->pin_fluid_boundary(outer);
//   CRDoubleMatrix expected_jacobian;
//   CRDoubleMatrix expected_mass_matrix;
//   problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
//                                         expected_jacobian);
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->pin_solid();
//   linear_problem_pt->pin_fluid_boundary(upper);
//   linear_problem_pt->pin_fluid_boundary(lower);
//   // linear_problem_pt->pin_fluid_boundary(inner);
//   linear_problem_pt->pin_fluid_boundary(outer);
//   CRDoubleMatrix actual_jacobian;
//   CRDoubleMatrix actual_mass_matrix;
//   linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
//                                                actual_jacobian);
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//   BOOST_TEST(jacobians_are_equal);
//   bool mass_matrix_are_equal =
//     compare_matrices(expected_mass_matrix, actual_mass_matrix);
//   BOOST_TEST(mass_matrix_are_equal);
// }
//
// BOOST_AUTO_TEST_CASE(compare_jacobians_and_mass_matrix_interior_and_lower_fluid)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->steady_newton_solve();
//   problem_pt->pin_solid();
//   problem_pt->pin_fluid_boundary(upper);
//   // problem_pt->pin_fluid_boundary(lower);
//   problem_pt->pin_fluid_boundary(inner);
//   problem_pt->pin_fluid_boundary(outer);
//   CRDoubleMatrix expected_jacobian;
//   CRDoubleMatrix expected_mass_matrix;
//   problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
//                                         expected_jacobian);
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->pin_solid();
//   linear_problem_pt->pin_fluid_boundary(upper);
//   // linear_problem_pt->pin_fluid_boundary(lower);
//   linear_problem_pt->pin_fluid_boundary(inner);
//   linear_problem_pt->pin_fluid_boundary(outer);
//   CRDoubleMatrix actual_jacobian;
//   CRDoubleMatrix actual_mass_matrix;
//   linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
//                                                actual_jacobian);
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//   BOOST_TEST(jacobians_are_equal);
//   bool mass_matrix_are_equal =
//     compare_matrices(expected_mass_matrix, actual_mass_matrix);
//   BOOST_TEST(mass_matrix_are_equal);
// }
//
// // TEST FULL JACOBIANS
// -------------------------------------------------------
//
// BOOST_AUTO_TEST_CASE(compare_matrices_pinned_boundaries)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   // Slip_Parameters::wall_velocity = 0.01;
//   problem_pt->steady_newton_solve();
//   problem_pt->reset_lagrange();
//   problem_pt->assign_initial_values_impulsive();
//   problem_pt->pin_solid_boundaries(inner);
//   problem_pt->pin_solid_boundaries(upper);
//   problem_pt->pin_solid_boundaries(outer);
//   problem_pt->pin_solid_boundaries(lower);
//   problem_pt->pin_fluid_boundary(upper);
//   problem_pt->pin_fluid_boundary(inner);
//   problem_pt->pin_fluid_boundary(outer);
//   problem_pt->pin_fluid_boundary(lower);
//   problem_pt->pin_kinematic_lagrange_multiplier(0.0);
//   problem_pt->pin_interior_pressure();
//
//   CRDoubleMatrix expected_jacobian;
//   CRDoubleMatrix expected_mass_matrix;
//   problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
//                                         expected_jacobian);
//
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->assign_initial_values_impulsive();
//   linear_problem_pt->pin_solid_boundaries(inner);
//   linear_problem_pt->pin_solid_boundaries(upper);
//   linear_problem_pt->pin_solid_boundaries(outer);
//   linear_problem_pt->pin_solid_boundaries(lower);
//   linear_problem_pt->pin_fluid_boundary(lower);
//   linear_problem_pt->pin_fluid_boundary(upper);
//   linear_problem_pt->pin_fluid_boundary(inner);
//   linear_problem_pt->pin_fluid_boundary(outer);
//   linear_problem_pt->set_constant_lagrange_free_surface_boundary_condition(0.0,
//                                                                            0.0);
//   linear_problem_pt->pin_interior_pressure();
//   CRDoubleMatrix actual_jacobian;
//   CRDoubleMatrix actual_mass_matrix;
//   linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
//                                                actual_jacobian);
//
//   save_dofs_types(problem_pt, "expected_dofs.txt");
//   save_dofs_types(linear_problem_pt, "actual_dofs.txt");
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//   BOOST_TEST(jacobians_are_equal);
//   bool mass_matrix_are_equal =
//     compare_matrices(expected_mass_matrix, actual_mass_matrix);
//   BOOST_TEST(mass_matrix_are_equal);
// }
//
// BOOST_AUTO_TEST_CASE(compare_matrices_test_kinematic_condition_pinned_lagrange)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   // Slip_Parameters::wall_velocity = 0.01;
//   problem_pt->steady_newton_solve();
//   problem_pt->reset_lagrange();
//   problem_pt->assign_initial_values_impulsive();
//   problem_pt->pin_solid_boundaries(inner);
//   problem_pt->pin_solid_boundaries(upper);
//   problem_pt->pin_solid_boundaries(outer);
//   problem_pt->pin_fluid_boundary(upper);
//   problem_pt->pin_fluid_boundary(inner);
//   problem_pt->pin_fluid_boundary(outer);
//   problem_pt->pin_kinematic_lagrange_multiplier(0.0);
//
//   CRDoubleMatrix expected_jacobian;
//   CRDoubleMatrix expected_mass_matrix;
//   problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
//                                         expected_jacobian);
//
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->assign_initial_values_impulsive();
//   linear_problem_pt->pin_solid_boundaries(inner);
//   linear_problem_pt->pin_solid_boundaries(upper);
//   linear_problem_pt->pin_solid_boundaries(outer);
//   linear_problem_pt->pin_fluid_boundary(upper);
//   linear_problem_pt->pin_fluid_boundary(inner);
//   linear_problem_pt->pin_fluid_boundary(outer);
//   linear_problem_pt->set_constant_lagrange_free_surface_boundary_condition(0.0,
//                                                                            0.0);
//   CRDoubleMatrix actual_jacobian;
//   CRDoubleMatrix actual_mass_matrix;
//   linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
//                                                actual_jacobian);
//
//   save_dofs_types(problem_pt, "expected_dofs.txt");
//   save_dofs_types(linear_problem_pt, "actual_dofs.txt");
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//   BOOST_TEST(jacobians_are_equal);
//   bool mass_matrix_are_equal =
//     compare_matrices(expected_mass_matrix, actual_mass_matrix);
//   BOOST_TEST(mass_matrix_are_equal);
// }
//
// BOOST_AUTO_TEST_CASE(compare_matrices_test_kinematic_condition)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   Slip_Parameters::wall_velocity = 0.01;
//   problem_pt->steady_newton_solve();
//   problem_pt->reset_lagrange();
//   problem_pt->assign_initial_values_impulsive();
//   problem_pt->pin_solid_boundaries(inner);
//   problem_pt->pin_solid_boundaries(upper);
//   problem_pt->pin_solid_boundaries(outer);
//   problem_pt->pin_fluid_boundary(upper);
//   problem_pt->pin_fluid_boundary(inner);
//   problem_pt->pin_fluid_boundary(outer);
//
//   CRDoubleMatrix expected_jacobian;
//   CRDoubleMatrix expected_mass_matrix;
//   problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
//                                         expected_jacobian);
//
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->assign_initial_values_impulsive();
//   linear_problem_pt->pin_solid_boundaries(inner);
//   linear_problem_pt->pin_solid_boundaries(upper);
//   linear_problem_pt->pin_solid_boundaries(outer);
//   linear_problem_pt->pin_fluid_boundary(upper);
//   linear_problem_pt->pin_fluid_boundary(inner);
//   linear_problem_pt->pin_fluid_boundary(outer);
//   CRDoubleMatrix actual_jacobian;
//   CRDoubleMatrix actual_mass_matrix;
//   linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
//                                                actual_jacobian);
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//   BOOST_TEST(jacobians_are_equal);
//   bool mass_matrix_are_equal =
//     compare_matrices(expected_mass_matrix, actual_mass_matrix);
//   BOOST_TEST(mass_matrix_are_equal);
// }
//
// BOOST_AUTO_TEST_CASE(compare_matrices_test_all_no_slip)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   Slip_Parameters::wall_velocity = 0.01;
//   Slip_Parameters::slip_length = -1;
//   problem_pt->steady_newton_solve();
//   problem_pt->reset_lagrange();
//   problem_pt->assign_initial_values_impulsive();
//
//   CRDoubleMatrix expected_jacobian;
//   CRDoubleMatrix expected_mass_matrix;
//   problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
//                                         expected_jacobian);
//
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->assign_initial_values_impulsive();
//
//   CRDoubleMatrix actual_jacobian;
//   CRDoubleMatrix actual_mass_matrix;
//   linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
//                                                actual_jacobian);
//
//   save_dofs_types(problem_pt, "expected_dofs.txt");
//   save_dofs_types(linear_problem_pt, "actual_dofs.txt");
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//   BOOST_TEST(jacobians_are_equal);
//   bool mass_matrix_are_equal =
//     compare_matrices(expected_mass_matrix, actual_mass_matrix);
//   BOOST_TEST(mass_matrix_are_equal);
//
//   output_matrices(expected_jacobian, "expected_jacobian.dat");
//   output_matrices(actual_jacobian, "actual_jacobian.dat");
// }

BOOST_AUTO_TEST_CASE(compare_eigenvalues)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  Slip_Parameters::slip_length = 1;
  problem_pt->steady_newton_solve();
  problem_pt->reset_lagrange();
  problem_pt->assign_initial_values_impulsive();
  // problem_pt->doc_solution();

  // CRDoubleMatrix expected_jacobian;
  // CRDoubleMatrix expected_mass_matrix;
  // problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
  //                                       expected_jacobian);


  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt);
  linear_problem_pt->assign_initial_values_impulsive();
  // linear_problem_pt->doc_solution();

  // CRDoubleMatrix actual_jacobian;
  // CRDoubleMatrix actual_mass_matrix;
  // linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
  //                                              actual_jacobian);

  // const unsigned N = expected_jacobian.nrow();
  // const unsigned M = expected_jacobian.ncol();
  // BOOST_TEST(actual_jacobian.nrow() == N);
  // BOOST_TEST(actual_jacobian.ncol() == M);
  //  save_dofs_types(problem_pt, "expected_dofs.txt");
  //  save_dofs_types(linear_problem_pt, "actual_dofs.txt");

  // bool jacobians_are_equal =
  //   compare_matrices(expected_jacobian, actual_jacobian);
  // bool mass_matrix_are_equal =
  //   compare_matrices(expected_mass_matrix, actual_mass_matrix);
  //  Note: Matrices are not expected to be the same. I have divided the
  //  forcing by r to be able to pin the lagrange multiplier to zero.
  //  BOOST_TEST(jacobians_are_equal);
  //  BOOST_TEST(mass_matrix_are_equal);
  //  output_matrices(expected_jacobian, "expected_jacobian.dat");
  //  output_matrices(actual_jacobian, "actual_jacobian.dat");

  Vector<std::complex<double>> expected_eigenvalue =
    problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  Vector<std::complex<double>> actual_eigenvalue =
    linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  BOOST_TEST(abs(actual_eigenvalue[0].real() - expected_eigenvalue[0].real()) <
             1e-6);
}


BOOST_AUTO_TEST_CASE(compare_eigenvalues_60degrees)
{
  Global_Physical_Parameters::Equilibrium_contact_angle =
    60.0 * oomph::MathematicalConstants::Pi / 180.0;
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  Slip_Parameters::slip_length = 1;
  problem_pt->steady_newton_solve();
  problem_pt->reset_lagrange();
  problem_pt->assign_initial_values_impulsive();
  // problem_pt->doc_solution();

  // CRDoubleMatrix expected_jacobian;
  // CRDoubleMatrix expected_mass_matrix;
  // problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
  //                                       expected_jacobian);


  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt);
  linear_problem_pt->assign_initial_values_impulsive();
  // linear_problem_pt->doc_solution();

  // CRDoubleMatrix actual_jacobian;
  // CRDoubleMatrix actual_mass_matrix;
  // linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
  //                                              actual_jacobian);

  // const unsigned N = expected_jacobian.nrow();
  // const unsigned M = expected_jacobian.ncol();
  // BOOST_TEST(actual_jacobian.nrow() == N);
  // BOOST_TEST(actual_jacobian.ncol() == M);
  //  save_dofs_types(problem_pt, "expected_dofs.txt");
  //  save_dofs_types(linear_problem_pt, "actual_dofs.txt");

  // bool jacobians_are_equal =
  //   compare_matrices(expected_jacobian, actual_jacobian, 5e-5);
  // bool mass_matrix_are_equal =
  //   compare_matrices(expected_mass_matrix, actual_mass_matrix);
  //  Note: Matrices are not expected to be the same. I have divided the
  //  forcing by r to be able to pin the lagrange multiplier to zero.
  //  BOOST_TEST(jacobians_are_equal);
  //  BOOST_TEST(mass_matrix_are_equal);
  //  output_matrices(expected_jacobian, "expected_jacobian.dat");
  //  output_matrices(actual_jacobian, "actual_jacobian.dat");

  Vector<std::complex<double>> expected_eigenvalue =
    problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  Vector<std::complex<double>> actual_eigenvalue =
    linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  BOOST_TEST(abs(actual_eigenvalue[0].real() - expected_eigenvalue[0].real()) <
             1e-5);
  Global_Physical_Parameters::Equilibrium_contact_angle =
    90.0 * oomph::MathematicalConstants::Pi / 180.0;
}

BOOST_AUTO_TEST_CASE(compare_eigenvalues_120degrees)
{
  Global_Physical_Parameters::Equilibrium_contact_angle =
    120.0 * oomph::MathematicalConstants::Pi / 180.0;
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  Slip_Parameters::slip_length = 1;
  problem_pt->steady_newton_solve();
  problem_pt->reset_lagrange();
  problem_pt->assign_initial_values_impulsive();
  // problem_pt->doc_solution();

  // CRDoubleMatrix expected_jacobian;
  // CRDoubleMatrix expected_mass_matrix;
  // problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
  //                                       expected_jacobian);


  shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt);
  linear_problem_pt->assign_initial_values_impulsive();
  // linear_problem_pt->doc_solution();

  //  CRDoubleMatrix actual_jacobian;
  //  CRDoubleMatrix actual_mass_matrix;
  //  linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
  //                                               actual_jacobian);
  // const unsigned N = expected_jacobian.nrow();
  // const unsigned M = expected_jacobian.ncol();
  // BOOST_TEST(actual_jacobian.nrow() == N);
  // BOOST_TEST(actual_jacobian.ncol() == M);
  // save_dofs_types(problem_pt, "expected_dofs.txt");
  // save_dofs_types(linear_problem_pt, "actual_dofs.txt");

  // Note: Matrices are not expected to be the same. I have divided the
  // forcing by r to be able to pin the lagrange multiplier to zero.
  // bool jacobians_are_equal =
  //  compare_matrices(expected_jacobian, actual_jacobian, 5e-5);
  // BOOST_TEST(jacobians_are_equal);
  // bool mass_matrix_are_equal =
  //  compare_matrices(expected_mass_matrix, actual_mass_matrix);
  // BOOST_TEST(mass_matrix_are_equal);
  //  output_matrices(expected_jacobian, "expected_jacobian.dat");
  //  output_matrices(actual_jacobian, "actual_jacobian.dat");

  Vector<std::complex<double>> expected_eigenvalue =
    problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  Vector<std::complex<double>> actual_eigenvalue =
    linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  BOOST_TEST(abs(actual_eigenvalue[0].real() - expected_eigenvalue[0].real()) <
             1e-2);
  Global_Physical_Parameters::Equilibrium_contact_angle =
    90.0 * oomph::MathematicalConstants::Pi / 180.0;
}

// Note: Matrices are not expected to be the same. I have divided the
// forcing by r to be able to pin the lagrange multiplier to zero.
// BOOST_AUTO_TEST_CASE(unsteady_compare_matrices)
// {
//   Flux_Parameters::Withdraw_speed = 0.1;
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   // problem_pt->steady_newton_solve();
//   problem_pt->make_unsteady();
//   problem_pt->reset_lagrange();
//   problem_pt->assign_initial_values_impulsive();
//
//
//   CRDoubleMatrix expected_jacobian;
//   CRDoubleMatrix expected_mass_matrix;
//   problem_pt->get_eigenproblem_matrices(expected_mass_matrix,
//                                         expected_jacobian);
//
//   shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
//     createLinearProblem(problem_pt);
//   linear_problem_pt->make_unsteady();
//   linear_problem_pt->assign_initial_values_impulsive();
//
//   CRDoubleMatrix actual_jacobian;
//   CRDoubleMatrix actual_mass_matrix;
//   linear_problem_pt->get_eigenproblem_matrices(actual_mass_matrix,
//                                                actual_jacobian);
//
//   save_dofs_types(problem_pt, "expected_dofs.txt");
//   save_dofs_types(linear_problem_pt, "actual_dofs.txt");
//
//   const unsigned N = expected_jacobian.nrow();
//   const unsigned M = expected_jacobian.ncol();
//
//   BOOST_TEST(actual_jacobian.nrow() == N);
//   BOOST_TEST(actual_jacobian.ncol() == M);
//
//   bool jacobians_are_equal =
//     compare_matrices(expected_jacobian, actual_jacobian);
//   BOOST_TEST(jacobians_are_equal);
//   bool mass_matrix_are_equal =
//     compare_matrices(expected_mass_matrix, actual_mass_matrix);
//   BOOST_TEST(mass_matrix_are_equal);
//
//   output_matrices(expected_jacobian, "expected_jacobian.dat");
//   output_matrices(actual_jacobian, "actual_jacobian.dat");
// }

// ***
// Function definitions

namespace oomph
{
  std::shared_ptr<AXISYM_PROBLEM> createBaseProblem()
  {
    Slip_Parameters::wall_velocity = 0.0;
    Slip_Parameters::slip_length = 5e-3;

    oomph_info.stream_pt() = &oomph_nullstream;
    const bool has_restart = false;
    return std::shared_ptr<AXISYM_PROBLEM>(new AXISYM_PROBLEM(
      Global_Physical_Parameters::Equilibrium_contact_angle, has_restart));
  }

  std::shared_ptr<PERTURBED_PROBLEM> createLinearProblem(
    std::shared_ptr<AXISYM_PROBLEM> base_problem_pt,
    const unsigned& azimuthal_mode_number)
  {
    return std::shared_ptr<PERTURBED_PROBLEM>(
      new PERTURBED_PROBLEM(base_problem_pt->bulk_mesh_pt(),
                            base_problem_pt->free_surface_mesh_pt(),
                            base_problem_pt->slip_surface_mesh_pt(),
                            azimuthal_mode_number));
  }
} // namespace oomph
