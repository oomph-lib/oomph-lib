#include "generic.h"

#include "constitutive/constitutive_laws.h"
#include "constitutive/plastic_constitutive_laws.h"

#include "solid/solid_plastic_elements.h"

#include "constitutive.h"

#include "meshes/simple_cubic_mesh.h"

using namespace std;
using namespace oomph;

// Define the solid mesh
template<class ELEMENT>
class ElasticSimpleCubicMesh : public virtual SimpleCubicMesh<ELEMENT>,
                               public virtual SolidMesh,
                               public virtual RefineableBrickMesh<ELEMENT>
{
public:
  ElasticSimpleCubicMesh(
    const unsigned nx,
    const unsigned ny,
    const unsigned nz,
    TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
    : SimpleCubicMesh<ELEMENT>(nx, ny, nz, 1.0, 1.0, 1.0, time_stepper_pt),
      SolidMesh(),
      RefineableBrickMesh<ELEMENT>()
  {
    set_lagrangian_nodal_coordinates();
    this->setup_octree_forest();
  }
};

namespace Parameters
{
  const double max_loading_height = 2.0;
  const unsigned n_step_loading = 30;
  const unsigned n_step_unloading = 120;
  const unsigned print_freq = 2;

  double Pinned_height = 1.0;
  bool controlling_height = true;
  const unsigned dirichlet_boundary_top = 0;
  const unsigned dirichlet_boundary_bottom = 5;
  const double dt = 0.1;

  double E = 1.0;
  double nu = 0.3;
  double mu = 0.384615;
  double lambda = 0.576923;

  StrainEnergyFunction* elastic_strain_energy_function_pt =
    new ModifiedNeoHookean(&lambda, &mu);
  ConstitutiveLaw* constitutive_law_pt =
    new IsotropicStrainEnergyFunctionConstitutiveLaw(
      elastic_strain_energy_function_pt);

  double isotropic_hardening_factor = std::sqrt(3.0 / 2.0);
  double isotropic_hardening_f0 = 500e-3;
  double isotropic_hardening_h1 = 0.8;
  double isotropic_hardening_h2 = 50;

  double zero = 0.0;
  double kinematic_hardening_stress_c = 0.1;

  double elastic_core_stress_c = 0.1;

  double normal_yield_ratio_elastic = 0.2;
  double normal_yield_ratio_constant_u = 200;
  double normal_yield_ratio_elasic_core_u = 1.0;

  double eta_p = 0.5;

  double kinematic_hardening_b = 0.1;
  double kinematic_hardening_eta = 0.5;

  double elastic_core_x = 0.1;
  double elastic_core_eta = 0.5;

  PlasticConstitutiveLaw* plastic_constitutive_law_pt =
    new PlasticConstitutiveLaw();

  const unsigned max_adapt = 1;
  const double min_error = 1e-3;
  const double max_error = 1e-1;
  Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
} // namespace Parameters

template<class ELEMENT>
class CubicAnisotropicSolidProblem : public Problem
{
public:
  CubicAnisotropicSolidProblem(const unsigned& nx,
                               const unsigned& ny,
                               const unsigned& nz);
  ~CubicAnisotropicSolidProblem() {}
  void actions_before_newton_solve();
  void actions_after_newton_solve();
  void unpin_bottom_boundary();
  void doc_solution(const unsigned& nplot);
};

template<class ELEMENT>
CubicAnisotropicSolidProblem<ELEMENT>::CubicAnisotropicSolidProblem(
  const unsigned& nx, const unsigned& ny, const unsigned& nz)
{
  // Problem::max_residuals() = 100.0;
  // Problem::max_newton_iterations() = 50;

  Problem::add_time_stepper_pt(new Newmark<2>);

  Problem::mesh_pt() =
    new ElasticSimpleCubicMesh<ELEMENT>(nx, ny, nz, Problem::time_stepper_pt());

  dynamic_cast<ElasticSimpleCubicMesh<ELEMENT>*>(Problem::mesh_pt())
    ->spatial_error_estimator_pt() = Parameters::error_estimator_pt;
  dynamic_cast<ElasticSimpleCubicMesh<ELEMENT>*>(Problem::mesh_pt())
    ->min_permitted_error() = Parameters::min_error;
  dynamic_cast<ElasticSimpleCubicMesh<ELEMENT>*>(Problem::mesh_pt())
    ->max_permitted_error() = Parameters::max_error;

  // Finish setting up the plastic constitutive law
  Parameters::plastic_constitutive_law_pt->isotropic_hardening_law_pt =
    new ExponentialIsotropicHardeningLaw(
      &Parameters::isotropic_hardening_factor,
      &Parameters::isotropic_hardening_f0,
      &Parameters::isotropic_hardening_h1,
      &Parameters::isotropic_hardening_h2);
  Parameters::plastic_constitutive_law_pt->yield_criterion_pt =
    new VonMisesYieldCriterion();
  Parameters::plastic_constitutive_law_pt->kinematic_hardening_law_pt =
    new IsotropicStrainEnergyFunctionConstitutiveLaw(new ModifiedNeoHookean(
      &Parameters::zero, &Parameters::kinematic_hardening_stress_c));
  Parameters::plastic_constitutive_law_pt->elastic_core_law_pt =
    new IsotropicStrainEnergyFunctionConstitutiveLaw(new ModifiedNeoHookean(
      &Parameters::zero, &Parameters::elastic_core_stress_c));
  Parameters::plastic_constitutive_law_pt->normal_yield_ratio_law_pt =
    new NormalYieldRatioLaw(&Parameters::normal_yield_ratio_elastic,
                            &Parameters::normal_yield_ratio_constant_u,
                            &Parameters::normal_yield_ratio_elasic_core_u);
  Parameters::plastic_constitutive_law_pt->eta_p_pt = &Parameters::eta_p;
  Parameters::plastic_constitutive_law_pt->kinematic_hardening_b_pt =
    &Parameters::kinematic_hardening_b;
  Parameters::plastic_constitutive_law_pt->kinematic_hardening_eta_pt =
    &Parameters::kinematic_hardening_eta;
  Parameters::plastic_constitutive_law_pt->elastic_core_x_pt =
    &Parameters::elastic_core_x;
  Parameters::plastic_constitutive_law_pt->elastic_core_eta_pt =
    &Parameters::elastic_core_eta;
  // END Finish setting up the plastic constitutive law


  const unsigned n_element = Problem::mesh_pt()->nelement();
  for (unsigned i = 0; i < n_element; i++)
  {
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Problem::mesh_pt()->element_pt(i));

    el_pt->plastic_model_type() =
      PlasticEquationsBase<3>::PlasticModel::Conventional;

    // Construct the internal data based on the plastic model type
    el_pt->construct_plastic_data();

    el_pt->constitutive_law_pt() = Parameters::constitutive_law_pt;

    el_pt->plastic_constitutive_law_pt() =
      Parameters::plastic_constitutive_law_pt;

    el_pt->assign_default_values_based_on_constitutive_law();

    el_pt->assign_plastic_timestepper_pt(Problem::time_stepper_pt(), true);
  }

  // Apply boundary conditions - pin all of them on one side
  unsigned n_bound_node =
    Problem::mesh_pt()->nboundary_node(Parameters::dirichlet_boundary_top);
  for (unsigned i = 0; i < n_bound_node; i++)
  {
    SolidNode* node_pt =
      dynamic_cast<SolidNode*>(Problem::mesh_pt()->boundary_node_pt(
        Parameters::dirichlet_boundary_top, i));

    for (unsigned i = 0; i < 3; i++)
    {
      node_pt->pin_position(i);
    }
  }

  n_bound_node =
    Problem::mesh_pt()->nboundary_node(Parameters::dirichlet_boundary_bottom);
  for (unsigned i = 0; i < n_bound_node; i++)
  {
    SolidNode* node_pt =
      dynamic_cast<SolidNode*>(Problem::mesh_pt()->boundary_node_pt(
        Parameters::dirichlet_boundary_bottom, i));

    for (unsigned i = 0; i < 3; i++)
    {
      node_pt->pin_position(i);
    }
  }

  initialise_dt(Parameters::dt);
  assign_initial_values_impulsive(Parameters::dt);

  cout << "There are " << assign_eqn_numbers() << " degrees of freedom"
       << std::endl;
}

template<class ELEMENT>
void CubicAnisotropicSolidProblem<ELEMENT>::actions_before_newton_solve()
{
  if (!Parameters::controlling_height) return;
  const unsigned n_bound_node =
    Problem::mesh_pt()->nboundary_node(Parameters::dirichlet_boundary_bottom);
  for (unsigned i = 0; i < n_bound_node; i++)
  {
    SolidNode* node_pt =
      dynamic_cast<SolidNode*>(Problem::mesh_pt()->boundary_node_pt(
        Parameters::dirichlet_boundary_bottom, i));

    node_pt->x(2) = Parameters::Pinned_height;
  }
}

template<class ELEMENT>
void CubicAnisotropicSolidProblem<ELEMENT>::actions_after_newton_solve()
{
  const unsigned n_node = Problem::mesh_pt()->nnode();
  for (unsigned l = 0; l < n_node; l++)
  {
    Node* node_pt = Problem::mesh_pt()->node_pt(l);
    node_pt->position_time_stepper_pt()->assign_initial_positions_impulsive(
      node_pt);
  }
}

template<class ELEMENT>
void CubicAnisotropicSolidProblem<ELEMENT>::unpin_bottom_boundary()
{
  Parameters::controlling_height = false;

  const unsigned n_bound_node =
    Problem::mesh_pt()->nboundary_node(Parameters::dirichlet_boundary_bottom);
  for (unsigned i = 0; i < n_bound_node; i++)
  {
    SolidNode* node_pt =
      dynamic_cast<SolidNode*>(Problem::mesh_pt()->boundary_node_pt(
        Parameters::dirichlet_boundary_bottom, i));

    for (unsigned i = 0; i < 3; i++)
    {
      node_pt->unpin_position(i);
    }
  }
  assign_eqn_numbers();
}

template<class ELEMENT>
void CubicAnisotropicSolidProblem<ELEMENT>::doc_solution(const unsigned& label)
{
  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts;
  npts = 5;

  // Output solution with specified number of plot points per element
  sprintf(filename, "RESLT_REFINEABLE_CUBE_HANG/soln%i.dat", label);
  some_file.open(filename);
  Problem::mesh_pt()->output(some_file, npts);
  some_file.close();
}

int main()
{
  unsigned nx = 2;
  unsigned ny = 2;
  unsigned nz = 2;

  CubicAnisotropicSolidProblem<RefineableQPlasticPVDElement<3, 2>> problem(
    nx, ny, nz);

  // Check whether the problem can be solved
  cout << "\n\n\nProblem self-test ";
  if (problem.self_test() == 0)
  {
    cout << "passed: Problem can be solved." << std::endl;
  }
  else
  {
    throw OomphLibError(
      "failed!", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  // Output solution
  unsigned out_num = 0;
  problem.doc_solution(out_num++);


  for (unsigned i = 0; i < Parameters::n_step_loading; i++)
  {
    // Manually move the bottom face of the cube
    const double s = (double)(i + 1) / (double)Parameters::n_step_loading;
    Parameters::Pinned_height =
      1.0 * (1.0 - s) + s * Parameters::max_loading_height;

    problem.unsteady_newton_solve(
      Parameters::dt, Parameters::max_adapt, false, true);
    if (i % Parameters::print_freq == 0) problem.doc_solution(out_num++);
  }

  problem.unpin_bottom_boundary();
  for (unsigned i = 0; i < Parameters::n_step_unloading; i++)
  {
    problem.unsteady_newton_solve(
      Parameters::dt, Parameters::max_adapt, false, true);
    if (i % Parameters::print_freq == 0) problem.doc_solution(out_num++);
  }
}