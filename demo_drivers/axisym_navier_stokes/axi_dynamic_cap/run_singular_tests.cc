#define BOOST_TEST_MODULE axisym_cap_test_module
#include <boost/test/included/unit_test.hpp>

#include "run_singular_tests.h"

using namespace std;
using namespace oomph;

// ***
// Test fix on a refined mesh
BOOST_AUTO_TEST_SUITE(fix)
BOOST_AUTO_TEST_CASE(singularity_without_fix)
{
  Global_Physical_Parameters::Equilibrium_contact_angle =
    120.0 * MathematicalConstants::Pi / 180.0;
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  problem_pt->set_bond_number(0.0);
  problem_pt->steady_newton_solve_adapt_if_needed(5);
  Vector<double> pressure = problem_pt->get_pressure_around_corner();
  for (unsigned i = 0; i < 8; i++)
  {
    cout << pressure[i] << " ";
  }
  cout << endl;
  problem_pt->doc_solution();

  // BOOST_TEST(jacobians_are_equal);
}

BOOST_AUTO_TEST_CASE(singularity_with_fix)
{
  Global_Physical_Parameters::Equilibrium_contact_angle =
    120.0 * MathematicalConstants::Pi / 180.0;
  shared_ptr<SINGULAR_PROBLEM> problem_pt = createSingularProblem();
  problem_pt->set_bond_number(0.0);
  problem_pt->steady_newton_solve_adapt_if_needed(5);
  Vector<double> pressure = problem_pt->get_pressure_around_corner();
  for (unsigned i = 0; i < 8; i++)
  {
    cout << pressure[i] << " ";
  }
  cout << endl;
  problem_pt->set_directory("RESLT2");
  problem_pt->doc_solution();

  // BOOST_TEST(jacobians_are_equal);
}
BOOST_AUTO_TEST_SUITE_END()

// ***
// Test the system close to 90 degrees
BOOST_AUTO_TEST_SUITE(close_to_90)
BOOST_AUTO_TEST_CASE(slightly_acute_normal)
{
  Global_Physical_Parameters::Equilibrium_contact_angle =
    89.0 * MathematicalConstants::Pi / 180.0;
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  problem_pt->set_bond_number(0.0);
  problem_pt->steady_newton_solve();
  Vector<double> pressure = problem_pt->get_pressure_around_corner();
  for (unsigned i = 0; i < 8; i++)
  {
    cout << pressure[i] << " ";
  }
  cout << endl;

  BOOST_TEST(abs(problem_pt->get_height_drop() + 0.0018670362255485331) < 1e-6);
}

BOOST_AUTO_TEST_CASE(slightly_acute_fix_fixed_c)
{
  Global_Physical_Parameters::Equilibrium_contact_angle =
    89.0 * MathematicalConstants::Pi / 180.0;
  shared_ptr<SINGULAR_PROBLEM> problem_pt = createSingularProblem();
  problem_pt->set_bond_number(0.0);
  problem_pt->fix_c(0.0);
  problem_pt->steady_newton_solve();
  Vector<double> pressure = problem_pt->get_pressure_around_corner();
  for (unsigned i = 0; i < 8; i++)
  {
    cout << pressure[i] << " ";
  }
  cout << endl;

  BOOST_TEST(abs(problem_pt->get_height_drop() + 0.0018670362255485331) < 1e-6);
}

BOOST_AUTO_TEST_CASE(slightly_acute_fix)
{
  Global_Physical_Parameters::Equilibrium_contact_angle =
    89.0 * MathematicalConstants::Pi / 180.0;
  shared_ptr<SINGULAR_PROBLEM> problem_pt = createSingularProblem();
  problem_pt->set_bond_number(0.0);
  problem_pt->steady_newton_solve();
  Vector<double> pressure = problem_pt->get_pressure_around_corner();
  for (unsigned i = 0; i < 8; i++)
  {
    cout << pressure[i] << " ";
  }
  cout << endl;

  BOOST_TEST(abs(problem_pt->get_height_drop() + 0.0018670362255485331) < 2e-5);
}

BOOST_AUTO_TEST_CASE(slightly_obtuse_normal)
{
  Global_Physical_Parameters::Equilibrium_contact_angle =
    91.0 * MathematicalConstants::Pi / 180.0;
  shared_ptr<AXISYM_PROBLEM> problem_pt = createBaseProblem();
  problem_pt->set_bond_number(0.0);
  problem_pt->steady_newton_solve();
  Vector<double> pressure = problem_pt->get_pressure_around_corner();
  for (unsigned i = 0; i < 8; i++)
  {
    cout << pressure[i] << " ";
  }
  cout << endl;

  BOOST_TEST(abs(problem_pt->get_height_drop() + 0.019074557268186201) < 1e-6);
}

BOOST_AUTO_TEST_CASE(slightly_obtuse_fix)
{
  Global_Physical_Parameters::Equilibrium_contact_angle =
    91.0 * MathematicalConstants::Pi / 180.0;
  shared_ptr<SINGULAR_PROBLEM> problem_pt = createSingularProblem();
  problem_pt->set_bond_number(0.0);
  problem_pt->steady_newton_solve();
  Vector<double> pressure = problem_pt->get_pressure_around_corner();
  for (unsigned i = 0; i < 8; i++)
  {
    cout << pressure[i] << " ";
  }
  cout << endl;

  problem_pt->doc_solution();

  BOOST_TEST(abs(problem_pt->get_height_drop() + 0.019074557268186201) < 4e-5);
}
BOOST_AUTO_TEST_SUITE_END()

// ***
// Test singular problem
BOOST_AUTO_TEST_CASE(singular_problem_creation_with_default_parameters)
{
  shared_ptr<SINGULAR_PROBLEM> problem_pt = createSingularProblem();
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
    compare_matrices(expected_jacobian, actual_jacobian, 1e-5);
  BOOST_TEST(jacobians_are_equal);
  output_matrices(actual_jacobian, "actual_jacobian.dat");
  save_dofs_types(problem_pt, "dofs.txt");
  problem_pt->steady_newton_solve();
  problem_pt->doc_solution();
  BOOST_TEST(abs(problem_pt->get_height_drop() + 0.29795664727577137) < 1e-6);
}

BOOST_AUTO_TEST_CASE(singular_problem_creation_obtuse)
{
  Global_Physical_Parameters::Equilibrium_contact_angle =
    120.0 * MathematicalConstants::Pi / 180.0;
  shared_ptr<SINGULAR_PROBLEM> problem_pt = createSingularProblem();
  // DoubleVector residuals;
  // DenseMatrix<double> expected_jacobian;
  // DenseDoubleMatrix actual_jacobian;
  // problem_pt->get_fd_jacobian(residuals, expected_jacobian);
  // problem_pt->get_jacobian(residuals, actual_jacobian);
  // const unsigned N = expected_jacobian.nrow();
  // const unsigned M = expected_jacobian.ncol();
  // BOOST_TEST(actual_jacobian.nrow() == N);
  // BOOST_TEST(actual_jacobian.ncol() == M);
  // bool jacobians_are_equal =
  //   compare_matrices(expected_jacobian, actual_jacobian, 2e-5);
  // BOOST_TEST(jacobians_are_equal);
  // output_matrices(actual_jacobian, "actual_jacobian.dat");
  // save_dofs_types(problem_pt, "dofs.txt");
  problem_pt->steady_newton_solve();
  problem_pt->doc_solution();
  BOOST_TEST(abs(problem_pt->get_height_drop() - 0.0) < 1e-6);
  Global_Physical_Parameters::Equilibrium_contact_angle =
    90.0 * MathematicalConstants::Pi / 180.0;
}

// ***
// Test linear problem from a singular problem
BOOST_AUTO_TEST_SUITE(linear_stability)
BOOST_AUTO_TEST_CASE(linear_singular_problem)
{
  shared_ptr<SINGULAR_PROBLEM> problem_pt = createSingularProblem();
  problem_pt->steady_newton_solve();
  std::shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt);
  Vector<std::complex<double>> eigenvalue =
    linear_problem_pt->solve_and_document_n_most_unstable_eigensolutions(1);
  printf("%16.11g \n", eigenvalue[0].real());
  BOOST_TEST(abs(eigenvalue[0].real() - (-0.24783965941010913)) < 1e-6);
}

BOOST_AUTO_TEST_CASE(
  fixed_c_90_degree_linear_stability_zero_wall_velocity_zero_c)
{
  // Expected solution
  shared_ptr<AXISYM_PROBLEM> regular_problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.0;
  regular_problem_pt->steady_newton_solve();
  std::cout << "Height drop: " << regular_problem_pt->get_height_drop()
            << std::endl;
  std::shared_ptr<PERTURBED_BASE_PROBLEM> linear_regular_problem_pt =
    createLinearBaseProblem(regular_problem_pt);
  Vector<std::complex<double>> expected_eigenvalue =
    linear_regular_problem_pt->solve_n_most_unstable_eigensolutions(1);

  // Actual solution
  shared_ptr<SINGULAR_PROBLEM> problem_pt = createSingularProblem();
  Slip_Parameters::wall_velocity = 0.0;
  problem_pt->fix_c(0.0);
  problem_pt->steady_newton_solve();
  std::cout << "Height drop: " << problem_pt->get_height_drop() << std::endl;
  std::shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt);
  Vector<std::complex<double>> actual_eigenvalue =
    linear_problem_pt->solve_n_most_unstable_eigensolutions(1);

  // Compare
  BOOST_TEST(abs(actual_eigenvalue[0].real() - expected_eigenvalue[0].real()) <
             1e-6);
}

BOOST_AUTO_TEST_CASE(fixed_c_90_degree_linear_stability_zero_wall_velocity)
{
  // Expected solution
  shared_ptr<AXISYM_PROBLEM> regular_problem_pt = createBaseProblem();
  Slip_Parameters::wall_velocity = 0.0;
  regular_problem_pt->steady_newton_solve();
  std::cout << "Height drop: " << regular_problem_pt->get_height_drop()
            << std::endl;
  std::shared_ptr<PERTURBED_BASE_PROBLEM> linear_regular_problem_pt =
    createLinearBaseProblem(regular_problem_pt);
  Vector<std::complex<double>> expected_eigenvalue =
    linear_regular_problem_pt->solve_n_most_unstable_eigensolutions(1);

  // Actual solution
  shared_ptr<SINGULAR_PROBLEM> problem_pt = createSingularProblem();
  Slip_Parameters::wall_velocity = 0.0;
  problem_pt->fix_c(1.0);
  problem_pt->steady_newton_solve();
  std::cout << "Height drop: " << problem_pt->get_height_drop() << std::endl;
  std::shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt);
  Vector<std::complex<double>> actual_eigenvalue =
    linear_problem_pt->solve_n_most_unstable_eigensolutions(1);

  // Compare
  BOOST_TEST(abs(actual_eigenvalue[0].real() - expected_eigenvalue[0].real()) <
             1e-6);
}

BOOST_AUTO_TEST_CASE(fixed_c_90_degree_linear_stability)
{
  // Expected solution
  shared_ptr<AXISYM_PROBLEM> regular_problem_pt = createBaseProblem();
  regular_problem_pt->steady_newton_solve();
  std::cout << "Height drop: " << regular_problem_pt->get_height_drop()
            << std::endl;
  std::shared_ptr<PERTURBED_BASE_PROBLEM> linear_regular_problem_pt =
    createLinearBaseProblem(regular_problem_pt);
  Vector<std::complex<double>> expected_eigenvalue =
    linear_regular_problem_pt->solve_n_most_unstable_eigensolutions(1);

  // Actual solution
  shared_ptr<SINGULAR_PROBLEM> problem_pt = createSingularProblem();
  problem_pt->fix_c(1.0);
  problem_pt->steady_newton_solve();
  std::cout << "Height drop: " << problem_pt->get_height_drop() << std::endl;
  std::shared_ptr<PERTURBED_PROBLEM> linear_problem_pt =
    createLinearProblem(problem_pt);
  Vector<std::complex<double>> actual_eigenvalue =
    linear_problem_pt->solve_n_most_unstable_eigensolutions(1);

  // Compare
  BOOST_TEST(abs(actual_eigenvalue[0].real() - expected_eigenvalue[0].real()) <
             1e-6);
}
BOOST_AUTO_TEST_SUITE_END()

// ***
// Height continuation
BOOST_AUTO_TEST_SUITE(height_continuation)
BOOST_AUTO_TEST_CASE(bo_height_control_singular_problem_creation_obtuse)
{
  Global_Physical_Parameters::Equilibrium_contact_angle =
    120.0 * MathematicalConstants::Pi / 180.0;
  BoHeightControlSingularAxisymDynamicCapProblem<
    SingularAxisymNavierStokesElement<
      HijackedProjectableAxisymmetricTTaylorHoodPVDElement>,
    BDF<2>>* problem_pt =
    new BoHeightControlSingularAxisymDynamicCapProblem<
      SingularAxisymNavierStokesElement<
        HijackedProjectableAxisymmetricTTaylorHoodPVDElement>,
      BDF<2>>(Global_Physical_Parameters::Equilibrium_contact_angle, 0);

  problem_pt->set_height_step_parameter_to_bond_number();

  problem_pt->pre_height_solve(0.1);

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
    compare_matrices(expected_jacobian, actual_jacobian, 2e-5);
  BOOST_TEST(jacobians_are_equal);
  output_matrices(actual_jacobian, "actual_jacobian.dat");
  save_dofs_types(problem_pt, "dofs.txt");

  problem_pt->post_height_solve();

  problem_pt->height_step_solve(1e-3);
  problem_pt->doc_solution();
  Global_Physical_Parameters::Equilibrium_contact_angle =
    90.0 * MathematicalConstants::Pi / 180.0;
}

BOOST_AUTO_TEST_CASE(test_ca_height_step_jacobian)
{
  Global_Physical_Parameters::Equilibrium_contact_angle =
    120.0 * MathematicalConstants::Pi / 180.0;
  CaHeightControlSingularAxisymDynamicCapProblem<
    SingularAxisymNavierStokesElement<
      HijackedProjectableAxisymmetricTTaylorHoodPVDElement>,
    BDF<2>>* problem_pt =
    new CaHeightControlSingularAxisymDynamicCapProblem<
      SingularAxisymNavierStokesElement<
        HijackedProjectableAxisymmetricTTaylorHoodPVDElement>,
      BDF<2>>(Global_Physical_Parameters::Equilibrium_contact_angle, 0);
  problem_pt->set_height_step_parameter_to_wall_velocity();

  problem_pt->pre_height_solve(0.1);

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
    compare_matrices(expected_jacobian, actual_jacobian, 2e-5);
  BOOST_TEST(jacobians_are_equal);
  output_matrices(actual_jacobian, "actual_jacobian.dat");
  save_dofs_types(problem_pt, "dofs.txt");

  problem_pt->post_height_solve();

  problem_pt->height_step_solve(1e-3);
  problem_pt->doc_solution();
  Global_Physical_Parameters::Equilibrium_contact_angle =
    90.0 * MathematicalConstants::Pi / 180.0;
}

BOOST_AUTO_TEST_CASE(ca_height_control_singular_problem_creation_obtuse)
{
  Global_Physical_Parameters::Equilibrium_contact_angle =
    120.0 * MathematicalConstants::Pi / 180.0;
  CaHeightControlSingularAxisymDynamicCapProblem<
    SingularAxisymNavierStokesElement<
      HijackedProjectableAxisymmetricTTaylorHoodPVDElement>,
    BDF<2>>* problem_pt =
    new CaHeightControlSingularAxisymDynamicCapProblem<
      SingularAxisymNavierStokesElement<
        HijackedProjectableAxisymmetricTTaylorHoodPVDElement>,
      BDF<2>>(Global_Physical_Parameters::Equilibrium_contact_angle, 0);
  // DoubleVector residuals;
  // DenseMatrix<double> expected_jacobian;
  // DenseDoubleMatrix actual_jacobian;
  // problem_pt->get_fd_jacobian(residuals, expected_jacobian);
  // problem_pt->get_jacobian(residuals, actual_jacobian);
  // const unsigned N = expected_jacobian.nrow();
  // const unsigned M = expected_jacobian.ncol();
  // BOOST_TEST(actual_jacobian.nrow() == N);
  // BOOST_TEST(actual_jacobian.ncol() == M);
  // bool jacobians_are_equal =
  //   compare_matrices(expected_jacobian, actual_jacobian, 2e-5);
  // BOOST_TEST(jacobians_are_equal);
  // output_matrices(actual_jacobian, "actual_jacobian.dat");
  // save_dofs_types(problem_pt, "dofs.txt");
  problem_pt->steady_newton_solve();
  problem_pt->set_height_step_parameter_to_wall_velocity();
  problem_pt->steady_newton_solve();
  problem_pt->height_step_solve(1e-3);
  problem_pt->doc_solution();
  Global_Physical_Parameters::Equilibrium_contact_angle =
    90.0 * MathematicalConstants::Pi / 180.0;
}
BOOST_AUTO_TEST_SUITE_END()

// ***
// Function definitions

namespace oomph
{
  std::shared_ptr<AXISYM_PROBLEM> createBaseProblem()
  {
    Slip_Parameters::wall_velocity = 0.01;
    Slip_Parameters::slip_length = 0.1;

    // oomph_info.stream_pt() = &oomph_nullstream;
    const bool has_restart = false;
    return std::shared_ptr<AXISYM_PROBLEM>(new AXISYM_PROBLEM(
      Global_Physical_Parameters::Equilibrium_contact_angle, has_restart));
  }

  std::shared_ptr<SINGULAR_PROBLEM> createSingularProblem()
  {
    Slip_Parameters::wall_velocity = 0.01;
    Slip_Parameters::slip_length = 0.1;
    Global_Physical_Parameters::Bo = 0.0;
    Global_Physical_Parameters::Equilibrium_contact_angle =
      120.0 * MathematicalConstants::Pi / 180.0;

    const bool has_restart = false;
    return std::shared_ptr<SINGULAR_PROBLEM>(new SINGULAR_PROBLEM(
      Global_Physical_Parameters::Equilibrium_contact_angle, has_restart));
  }

  std::shared_ptr<PERTURBED_PROBLEM> createLinearProblem(
    std::shared_ptr<SINGULAR_PROBLEM> base_problem_pt,
    const unsigned& azimuthal_mode_number)
  {
    return std::shared_ptr<PERTURBED_PROBLEM>(
      new PERTURBED_PROBLEM(base_problem_pt->bulk_mesh_pt(),
                            base_problem_pt->free_surface_mesh_pt(),
                            base_problem_pt->slip_surface_mesh_pt(),
                            azimuthal_mode_number));
  }

  std::shared_ptr<PERTURBED_BASE_PROBLEM> createLinearBaseProblem(
    std::shared_ptr<AXISYM_PROBLEM> base_problem_pt,
    const unsigned& azimuthal_mode_number)
  {
    return std::shared_ptr<PERTURBED_BASE_PROBLEM>(
      new PERTURBED_BASE_PROBLEM(base_problem_pt->bulk_mesh_pt(),
                                 base_problem_pt->free_surface_mesh_pt(),
                                 base_problem_pt->slip_surface_mesh_pt(),
                                 azimuthal_mode_number));
  }


} // namespace oomph
