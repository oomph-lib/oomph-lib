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

BOOST_AUTO_TEST_CASE(nonlinear_analytic_J_equals_fd_J_no_solve)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  Global_Physical_Parameters::Bo = 1.0;
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

BOOST_AUTO_TEST_CASE(nonlinear_analytic_J_equals_fd_J_gravity)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  // Slip_Parameters::wall_velocity = 0.01;
  Global_Physical_Parameters::Bo = 1.0;
  problem_pt->steady_newton_solve();
  DoubleVector residuals;
  DenseMatrix<double> expected_jacobian;
  DenseDoubleMatrix actual_jacobian;
  problem_pt->get_fd_jacobian(residuals, expected_jacobian);
  problem_pt->get_jacobian(residuals, actual_jacobian);

  std::ofstream output_stream;
  output_stream.open("RESLT/dofs.txt");
  problem_pt->describe_dofs(output_stream);
  output_stream.close();

  const unsigned N = expected_jacobian.nrow();
  const unsigned M = expected_jacobian.ncol();

  BOOST_TEST(actual_jacobian.nrow() == N);
  BOOST_TEST(actual_jacobian.ncol() == M);

  bool jacobians_are_equal =
    compare_matrices(expected_jacobian, actual_jacobian);
  BOOST_TEST(jacobians_are_equal);
}

BOOST_AUTO_TEST_CASE(nonlinear_analytic_J_equals_fd_J_gravity_60)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt =
    createBaseProblem(60.0 * MathematicalConstants::Pi / 180.0);
  // Slip_Parameters::wall_velocity = 0.01;
  Global_Physical_Parameters::Bo = 1.0;
  problem_pt->steady_newton_solve();
  DoubleVector residuals;
  DenseMatrix<double> expected_jacobian;
  DenseDoubleMatrix actual_jacobian;
  problem_pt->get_fd_jacobian(residuals, expected_jacobian);
  problem_pt->get_jacobian(residuals, actual_jacobian);

  std::ofstream output_stream;
  output_stream.open("RESLT/dofs.txt");
  problem_pt->describe_dofs(output_stream);
  output_stream.close();

  const unsigned N = expected_jacobian.nrow();
  const unsigned M = expected_jacobian.ncol();

  BOOST_TEST(actual_jacobian.nrow() == N);
  BOOST_TEST(actual_jacobian.ncol() == M);

  bool jacobians_are_equal =
    compare_matrices(expected_jacobian, actual_jacobian);
  BOOST_TEST(jacobians_are_equal);
}

BOOST_AUTO_TEST_CASE(debug_elemental_jacobian)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  problem_pt->steady_newton_solve();
  problem_pt->debug_elemental_jacobian();
}

BOOST_AUTO_TEST_CASE(nonlinear_analytic_J_equals_fd_J_ca)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  problem_pt->steady_newton_solve();
  problem_pt->debug_elemental_jacobian();
  DoubleVector residuals;
  DenseMatrix<double> expected_jacobian;
  DenseDoubleMatrix actual_jacobian;
  problem_pt->get_fd_jacobian(residuals, expected_jacobian);
  problem_pt->get_jacobian(residuals, actual_jacobian);

  std::ofstream output_stream;
  output_stream.open("RESLT/dofs.txt");
  problem_pt->describe_dofs(output_stream);
  output_stream.close();

  problem_pt->doc_solution();

  const unsigned N = expected_jacobian.nrow();
  const unsigned M = expected_jacobian.ncol();

  BOOST_TEST(actual_jacobian.nrow() == N);
  BOOST_TEST(actual_jacobian.ncol() == M);

  bool jacobians_are_equal =
    compare_matrices(expected_jacobian, actual_jacobian);
  BOOST_TEST(jacobians_are_equal);
}

BOOST_AUTO_TEST_CASE(nonlinear_analytic_J_equals_fd_J_ca_60)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt =
    createBaseProblem(60.0 * MathematicalConstants::Pi / 180.0);
  Slip_Parameters::wall_velocity = 0.01;
  problem_pt->steady_newton_solve();
  problem_pt->debug_elemental_jacobian();
  DoubleVector residuals;
  DenseMatrix<double> expected_jacobian;
  DenseDoubleMatrix actual_jacobian;
  problem_pt->get_fd_jacobian(residuals, expected_jacobian);
  problem_pt->get_jacobian(residuals, actual_jacobian);

  std::ofstream output_stream;
  output_stream.open("RESLT/dofs.txt");
  problem_pt->describe_dofs(output_stream);
  output_stream.close();

  problem_pt->doc_solution();

  const unsigned N = expected_jacobian.nrow();
  const unsigned M = expected_jacobian.ncol();

  BOOST_TEST(actual_jacobian.nrow() == N);
  BOOST_TEST(actual_jacobian.ncol() == M);

  bool jacobians_are_equal =
    compare_matrices(expected_jacobian, actual_jacobian);
  BOOST_TEST(jacobians_are_equal);
}

BOOST_AUTO_TEST_CASE(fd_solve)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.1;
  problem_pt->linear_solver_pt() = new FD_LU;
  problem_pt->steady_newton_solve();
  problem_pt->debug_elemental_jacobian();
  DoubleVector residuals;
  DenseMatrix<double> expected_jacobian;
  DenseDoubleMatrix actual_jacobian;
  problem_pt->get_fd_jacobian(residuals, expected_jacobian);
  problem_pt->get_jacobian(residuals, actual_jacobian);

  std::ofstream output_stream;
  output_stream.open("RESLT/dofs.txt");
  problem_pt->describe_dofs(output_stream);
  output_stream.close();

  const unsigned N = expected_jacobian.nrow();
  const unsigned M = expected_jacobian.ncol();

  BOOST_TEST(actual_jacobian.nrow() == N);
  BOOST_TEST(actual_jacobian.ncol() == M);

  bool jacobians_are_equal =
    compare_matrices(expected_jacobian, actual_jacobian);
  BOOST_TEST(jacobians_are_equal);
}


BOOST_AUTO_TEST_CASE(nonlinear_analytic_J_equals_fd_J)
{
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.01;
  Global_Physical_Parameters::Bo = 1.0;
  problem_pt->steady_newton_solve();
  DoubleVector residuals;
  DenseMatrix<double> expected_jacobian;
  DenseDoubleMatrix actual_jacobian;
  problem_pt->get_fd_jacobian(residuals, expected_jacobian);
  problem_pt->get_jacobian(residuals, actual_jacobian);

  std::ofstream output_stream;
  output_stream.open("RESLT/dofs.txt");
  problem_pt->describe_dofs(output_stream);
  output_stream.close();

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

// BOOST_AUTO_TEST_CASE(nonlinear_problem_solve_eigenproblem)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   problem_pt->steady_newton_solve();
//   problem_pt->reset_lagrange();
//   problem_pt->assign_initial_values_impulsive();
//   Vector<std::complex<double>> eigenvalue =
//     problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
//   printf("%16.11g \n", eigenvalue[0].real());
//   BOOST_TEST(abs(eigenvalue[0].real() - (-0.24783965941010913)) < 1e-6);
// }

// BOOST_AUTO_TEST_CASE(nonlinear_problem_solve_eigenproblem_with_fluid_flow)
// {
//   shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
//   Slip_Parameters::wall_velocity = 0.01;
//   problem_pt->steady_newton_solve();
//   problem_pt->reset_lagrange();
//   problem_pt->assign_initial_values_impulsive();
//   Vector<std::complex<double>> eigenvalue =
//     problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
//   printf("%16.11g \n", eigenvalue[0].real());
//   BOOST_TEST(abs(eigenvalue[0].real() - (-0.27630985141964381)) < 1e-6);
// }

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
  BOOST_TEST(abs(actual_centre_height - 0.0028304837510845986) < 1.1e-3);
}

// ***
// Function definitions

namespace oomph
{
  std::shared_ptr<AXISYM_PROBLEM> createBaseProblem(const double& contact_angle)
  {
    Slip_Parameters::wall_velocity = 0.0;
    Slip_Parameters::slip_length = 1.0;
    Global_Physical_Parameters::Equilibrium_contact_angle = contact_angle;

    // oomph_info.stream_pt() = &oomph_nullstream;
    const bool has_restart = false;
    return std::shared_ptr<AXISYM_PROBLEM>(new AXISYM_PROBLEM(
      Global_Physical_Parameters::Equilibrium_contact_angle, has_restart));
  }
} // namespace oomph
