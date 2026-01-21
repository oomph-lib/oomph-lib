// #include "generic/generic.h"

#include "constitutive/constitutive_laws.h"
#include "constitutive/plastic_constitutive_laws.h"

#include "solid/solid_elements.h"
#include "solid/solid_plastic_elements.h"


#include "solid_cubic_mesh.h"

#include "meshes/rectangular_quadmesh.h"

using namespace oomph;

#define PLASTIC

// Define Dimension here (2 or 3)
#define DIM 3

template<class ELEMENT>
class PlasticUniaxialStrainTestProblem : public Problem
{
private:
  SolidMesh* plasticMesh;

  // Boundary ID definitions change between Quad and Cubic meshes
#if DIM == 2
  // RectangularQuadMesh: 0:bottom, 1:right, 2:top, 3:left
  const unsigned pin_boundary = 3; // Left edge
  const unsigned moving_boundary = 1; // Right edge
#else
  const unsigned pin_boundary = 4; // Left face
  const unsigned moving_boundary = 2; // Right face
#endif

  double E = 100e9;
  double nu = 0.3;
  double density = 1000;

  // Properties of the cube / strain
  double l = 1e-1; // m
  double v = 0.001; // m / s

  double depsilon = 25e-4;
  // double depsilon = 55e-4;
  // double dt = 1e-2;
  double dt = depsilon * l / v;
  double tmax = dt * 0.1 / depsilon;

  // Choose things to normalize with
  double l0 = l;
  double rho0 = density;
  double E0 = E;
  double v0 = std::sqrt(E0 / rho0);
  double t0 = l0 / v0;

  // double l0 = 1;
  // double rho0 = 1;
  // double E0 = 1;
  // double v0 = 1;
  // double t0 = 1;

  double EOne = E / E0;
  double dtOne = dt / t0;
  double tmaxOne = tmax / t0;
  double vOne = v / v0;
  double lOne = l / l0;
  double rhoOne = density / rho0;

  ConstitutiveLaw* constitutive_law_pt;
  PlasticConstitutiveLaw* plastic_constitutive_law_pt;

  TimeStepper* my_time_stepper_pt;

  std::ofstream stress_strain_file;

  // Doc info.
  DocInfo doc_info;

  // Cyclic control parameters
  bool isLoading = false;
  double totalStrain = 0;
  double lastStartingStrain = 0;
  double lastStrainStartTime = 0;
  double initialStrain = 0.05;
  double consecutiveStrain = 0.25 * initialStrain;

  bool steady_solve = false;

  void setup()
  {
    if (steady_solve)
    {
      my_time_stepper_pt = new Steady<1>;
    }
    else
    {
      my_time_stepper_pt = new Newmark<2>;
    }
    add_time_stepper_pt(my_time_stepper_pt);

    // Create Mesh based on DIM
#if DIM == 2
    plasticMesh = new ElasticRectangularQuadMesh<ELEMENT>(
      1, 1, lOne, lOne, time_stepper_pt());
#else
    plasticMesh = new SolidCubicMesh<ELEMENT>(2,
                                              2,
                                              2, // Elements in X, Y, Z
                                              0.0,
                                              lOne, // X min, X max
                                              0.0,
                                              lOne, // Y min, Y max
                                              0.0,
                                              lOne, // Z min, Z max
                                              time_stepper_pt());
#endif

    oomph_info << "Created a mesh with " << plasticMesh->nelement()
               << " elements in " << DIM << "D." << std::endl;


    // Create elastic constitutive law
    constitutive_law_pt = new GeneralisedHookean(&nu, &EOne);

    plastic_constitutive_law_pt = new PlasticConstitutiveLaw();
    plastic_constitutive_law_pt->isotropic_hardening_f0 = 500e6 / E0;
    plastic_constitutive_law_pt->isotropic_hardening_h1 = 0.8;
    plastic_constitutive_law_pt->isotropic_hardening_h2 = 50;
    plastic_constitutive_law_pt->normal_yield_ratio_constant_u = 200;

    plastic_constitutive_law_pt->normal_yield_ratio_elastic = 0.2;

    plastic_constitutive_law_pt->eta_p = 0.5;

    plastic_constitutive_law_pt->kinematic_hardening_stress_c = 1e8 / E0;
    plastic_constitutive_law_pt->kinematic_hardening_b = 0.1;
    plastic_constitutive_law_pt->kinematic_hardening_eta = 0.5;

    plastic_constitutive_law_pt->elastic_core_stress_c = 1e8 / E0;
    plastic_constitutive_law_pt->elastic_core_x = 0.1;
    plastic_constitutive_law_pt->elastic_core_eta = 0.5;

    Max_newton_iterations = 50;
    Newton_solver_tolerance = 1e-8;
    Max_residuals = 1e10;
    Relaxation_factor = 1;
    Shut_up_in_newton_solve = 1;
    disable_info_in_newton_solve();

    // Apply the constitutive laws to the elements
    for (unsigned e = 0; e < plasticMesh->nelement(); e++)
    {
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(plasticMesh->element_pt(e));

      // Assign constitutive law
      el_pt->constitutive_law_pt() = constitutive_law_pt;

#ifdef PLASTIC
      el_pt->plastic_constitutive_law_pt() = plastic_constitutive_law_pt;

      el_pt->assign_plastic_timestepper(my_time_stepper_pt);

      el_pt->assign_default_values_based_on_constitutive_law();

      el_pt->plastic_newton_solver_tolerance() = Newton_solver_tolerance;

      // el_pt->enable_plastic_solve_by_fd();
#endif
      el_pt->lambda_sq_pt() = &rhoOne;
    }

    add_sub_mesh(plasticMesh);
    build_global_mesh();

    initialise_dt(dtOne);
    assign_initial_values_impulsive(dtOne);


    doc_info.set_directory("PlasticTestProblemOut");

    setup_initial_condition();

    // Do equation numbering
    oomph_info << "# of dofs " << assign_eqn_numbers() << std::endl;

    stress_strain_file.open("stress_strain_result.dat");

#if DIM == 2
    stress_strain_file << "Time Strain Stress_XX Stress_YY Nnewton_iter_taken"
                       << std::endl;
#else
    stress_strain_file
      << "Time Strain Stress_XX Stress_YY Stress_ZZ Nnewton_iter_taken"
      << std::endl;
#endif
  }

  void record_stress_strain()
  {
    // Get Stress
    DenseMatrix<double> sigma;
    getStress(sigma);

    // Get strain
    double strain = getStrain();

    // Write to file
    stress_strain_file << time_pt()->time() << " " << strain << " "
                       << sigma(0, 0) * E0 << " " // Sigma_xx
                       << sigma(1, 1) * E0; // Sigma_yy

#if DIM == 3
    stress_strain_file << " " << sigma(2, 2) * E0; // Sigma_zz
#endif

    stress_strain_file << " " << Nnewton_iter_taken;

    stress_strain_file << std::endl;
  }

  void output()
  {
    int currentOutputNumber = doc_info.number()++;
    bool writeRestart = true;

    // Output the tube shape in the most strongly collapsed configuration
    if (plasticMesh)
    {
      std::string filename;
      filename = doc_info.directory() + "/" + "FEM_" +
                 std::to_string(currentOutputNumber) + ".dat";
      std::ofstream file(filename);
      plasticMesh->output(file, 2);
      file.close();
      convert_dat_to_vtu(filename);
    }

    if (writeRestart)
    {
      // Now actually dump the whole problem to be able to restart
      std::string filenameRestart = doc_info.directory() + "/" + "FEM_" +
                                    std::to_string(currentOutputNumber) +
                                    ".restart";
      std::ofstream fileRestart(filenameRestart);

      // Write step number
      fileRestart << doc_info.number() << " # doc info number" << std::endl;

      const std::streamsize oldPrecision = fileRestart.precision();
      fileRestart.precision(std::numeric_limits<double>::digits10);
      Problem::dump(fileRestart);
      fileRestart.precision(oldPrecision);
    }
  }

  void convert_dat_to_vtu(const std::string& filename)
  {
    // Call external script to convert .dat file to .vtu file. With option -o to
    // overwrite existing file.
    const std::string cmd = "./oomph-convert.py -o " + filename;
    std::ignore = std::system(cmd.c_str());
  }

  void setup_initial_condition()
  {
    unsigned nnode = plasticMesh->nboundary_node(pin_boundary);
    for (unsigned l = 0; l < nnode; l++)
    {
      SolidNode* node_pt = dynamic_cast<SolidNode*>(
        plasticMesh->boundary_node_pt(pin_boundary, l));
      node_pt->pin_position(0);

#if DIM == 2
      if (l == 0)
      {
        node_pt->pin_position(1);
      }
#else
      if (std::abs(node_pt->x(1)) < 1e-9 && std::abs(node_pt->x(2)) < 1e-9)
      {
        node_pt->pin_position(1); // Pin Y
        node_pt->pin_position(2); // Pin Z
      }

      if (std::abs(node_pt->x(1) - 1.0) < 1e-9 &&
          std::abs(node_pt->x(2)) < 1e-9)
      {
        node_pt->pin_position(2);
      }

      oomph_info << "Node " << l << " on left boundary has coords "
                 << node_pt->x(0) << ", " << node_pt->x(1) << ", "
                 << node_pt->x(2) << std::endl;
#endif
    }

    nnode = plasticMesh->nboundary_node(moving_boundary);
    for (unsigned l = 0; l < nnode; l++)
    {
      SolidNode* node_pt = dynamic_cast<SolidNode*>(
        plasticMesh->boundary_node_pt(moving_boundary, l));

      node_pt->pin_position(0);
    }
  }

  void actions_before_newton_solve()
  {
    {
      const unsigned nnode = plasticMesh->nboundary_node(moving_boundary);
      for (unsigned l = 0; l < nnode; l++)
      {
        SolidNode* node_pt = dynamic_cast<SolidNode*>(
          plasticMesh->boundary_node_pt(moving_boundary, l));

        // Move X coordinate: Original Length (1.0) + time
        double totalStrain = 1 + lastStartingStrain;
        if (isLoading)
        {
          totalStrain +=
            depsilon * (time_pt()->time() - lastStrainStartTime) * t0;
        }
        else
        {
          totalStrain -=
            depsilon * (time_pt()->time() - lastStrainStartTime) * t0;
        }

        node_pt->x(0) = lOne * totalStrain;
      }
    }
  }

  double getStrain()
  {
    const unsigned nnode = plasticMesh->nboundary_node(moving_boundary);
    for (unsigned l = 0; l < nnode; l++)
    {
      SolidNode* node_pt = dynamic_cast<SolidNode*>(
        plasticMesh->boundary_node_pt(moving_boundary, l));

      return node_pt->x(0) / lOne - 1;
    }

    return 0.0;
  }

  void getStress(DenseMatrix<double>& sigma)
  {
    sigma.resize(DIM, DIM, 0.0);
    sigma.initialise(0.0);

    // Get the first element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(plasticMesh->element_pt(0));

    // Use the central integration point
    unsigned ipt = 0;

#ifdef PLASTIC
    // Get Stress
    el_pt->get_cauchy_stress(ipt, sigma);
#endif
  }


  void solve()
  {
    output();
    unsigned nsolve = 0;

    record_stress_strain();

    for (unsigned int e = 0; e < plasticMesh->nelement(); e++)
    {
      ELEMENT* elPt = dynamic_cast<ELEMENT*>(plasticMesh->element_pt(e));
      if (!elPt)
      {
        oomph_info << "OOps something went wrong " << elPt << std::endl;
        exit(0);
      }
      elPt->check_initial_condition();
      // exit(0);
    }

    unsigned int num_cycles = 3;

    for (unsigned int cycle = 0; cycle < num_cycles; cycle++)
    {
      oomph_info << "Cycle " << cycle + 1 << ": LOADING" << std::endl;

      // 1. Ensure nodes are PINNED
      isLoading = true;
      lastStrainStartTime = time_pt()->time();
      lastStartingStrain = getStrain();


      double targetStrain = initialStrain;
      if (cycle > 0)
      {
        targetStrain = consecutiveStrain;
      }
      while (getStrain() - lastStartingStrain < targetStrain)
      {
        solve_step();

        oomph_info << " Last solve took " << Nnewton_iter_taken << " iterations"
                   << std::endl;

#if defined(PLASTIC) && (DIM == 2)
        for (unsigned int e = 0; e < plasticMesh->nelement(); e++)
          dynamic_cast<ELEMENT*>(plasticMesh->element_pt(e))
            ->plastic_newton_solve();
#endif

        record_stress_strain();
        nsolve++;
        if (nsolve % 10 == 0) output();
      }

      oomph_info << "Cycle " << cycle + 1 << ": Unloading" << std::endl;

      isLoading = false;
      lastStartingStrain = getStrain();
      lastStrainStartTime = time_pt()->time();

      DenseMatrix<double> sigma;
      do
      {
        solve_step();

#if defined(PLASTIC) && (DIM == 2)
        for (unsigned int e = 0; e < plasticMesh->nelement(); e++)
          dynamic_cast<ELEMENT*>(plasticMesh->element_pt(e))
            ->plastic_newton_solve();
#endif

        record_stress_strain(); // Stress should drop to approx 0 here
        nsolve++;
        if (nsolve % 10 == 0)
        {
          output();
        }

        getStress(sigma);

      } while (sigma(0, 0) > 0);
    }
  }

  void solve_step()
  {
    if (steady_solve)
    {
      steady_newton_solve();
      time_pt()->time() += dtOne;
    }
    else
    {
      unsteady_newton_solve(dtOne);
    }
  }


public:
  void run()
  {
    setup();
    solve();
  }

  void set_steady_solve(bool steady)
  {
    steady_solve = steady;

    if (steady_solve)
      oomph_info << "Steady solve activated." << std::endl;
    else
      oomph_info << "Unsteady solve activated." << std::endl;
  }
};

int main(int argc, char* argv[])
{
  // Check, if a steady solve should be performed
  bool steady_analysis = false;

  for (int i = 1; i < argc; ++i)
  {
    std::string arg = argv[i];
    if (arg == "-s" || arg == "--steady")
    {
      steady_analysis = true;
    }
  }

// Select element based on DIM using simplified syntax
#ifdef PLASTIC
  PlasticUniaxialStrainTestProblem<QPlasticPVDElement<DIM, 2>> problem;
#else
  PlasticUniaxialStrainTestProblem<QPVDElement<DIM, 2>> problem;
#endif

  problem.set_steady_solve(steady_analysis);
  problem.run();
  return 0;
}