#define BOOST_TEST_MODULE obtuse_contact_angle_test_module
#include <boost/test/included/unit_test.hpp>

#include "run_tests.h"

using namespace std;
using namespace oomph;
BOOST_AUTO_TEST_SUITE(UtilityFunctions)
BOOST_AUTO_TEST_CASE(test_true)
{
  BOOST_TEST(true);
}

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

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(AxisymSectorProblemTest)

BOOST_AUTO_TEST_CASE(axisym_sector_problem)
{
  Parameters parameters;

  // Create problem
  AxisymSectorProblem<
    ProjectableAxisymmetricTaylorHoodElement<AxisymmetricTTaylorHoodElement>>
    problem;
  problem.pin_far_field_elements();

  // Unsteady problem
  const double dt = 2e-2;
  const double ft = 5 * dt;
  const unsigned nt = std::ceil(ft / dt);

  for (unsigned it = 0; it < nt; it++)
  {
    std::cout << "unsteady_newton_solve" << std::endl;
    problem.unsteady_newton_solve(dt);
    problem.doc_solution();
  }

  // Steady problem
  // problem.make_steady();
  problem.steady_newton_solve();
  problem.doc_solution();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SingularAxisymSectorProblemTest)

// BOOST_AUTO_TEST_CASE(singular_axisym_sector_problem)
//{
//   // Set parameter values
//   Parameters parameters;
//
//   // Create problem
//   SingularAxisymSectorProblem<AxisymNavierStokesElementWithSingularity<
//     Hijacked<ProjectableAxisymmetricTaylorHoodElement<
//       PseudoSolidNodeUpdateElement<AxisymmetricTTaylorHoodElement,
//                                    TPVDElement<2, 3>>>>>>
//     problem;
//
//   // Unsteady problem
//   const double dt = 2e-2;
//   const double ft = 5 * dt;
//   const unsigned nt = std::ceil(ft / dt);
//
//   for (unsigned it = 0; it < nt; it++)
//   {
//     std::cout << "unsteady_newton_solve" << std::endl;
//     problem.unsteady_newton_solve(dt);
//     problem.doc_solution();
//   }
//
//   // Steady problem
//   problem.make_steady();
//   problem.steady_newton_solve();
//   problem.doc_solution();
// }
//
// BOOST_AUTO_TEST_CASE(debug)
//{
//   // Set parameter values
//   Parameters parameters;
//
//   // Create problem
//   SingularAxisymSectorProblem<AxisymNavierStokesElementWithSingularity<
//     Hijacked<ProjectableAxisymmetricTaylorHoodElement<
//       PseudoSolidNodeUpdateElement<AxisymmetricTTaylorHoodElement,
//                                    TPVDElement<2, 3>>>>>>
//     problem;
//
//   //problem.debug_residuals();
//   problem.debug_jacobian();
// }

BOOST_AUTO_TEST_SUITE_END()
