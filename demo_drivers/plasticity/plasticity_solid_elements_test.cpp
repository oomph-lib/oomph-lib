#include "constitutive/constitutive_laws.h"
#include "constitutive/plastic_constitutive_laws.h"

#include "solid/solid_elements.h"
#include "solid/solid_plastic_elements.h"

#include "solid_cubic_mesh.h"

#include "meshes/rectangular_quadmesh.h"

using namespace oomph;

namespace Traction
{
  double pressure = 0.0;
  void traction(const Vector<double>& xi,
                const Vector<double>& x,
                const Vector<double>& n,
                Vector<double>& traction)
  {
    for (unsigned i = 0; i < 2; i++)
    {
      traction[i] = n[i] * pressure;
    }
  }
}

namespace BodyForce
{
 double tMax = 2.5;
 double tZero = 10.0;
 double gMax = 1.5;
 void body_force(const double& t,
                 const Vector<double>& xi,
                 Vector<double>& b)
 {
   if(t<=tMax) b[0] = t * gMax/tMax;
   else b[0] = gMax * (1.0 - (t - tMax) / (tZero - tMax));

   if(b[0]>tZero) b[0] = 0.0;
  }
}

template<class ELEMENT>
class PlasticTestProblem : public Problem
{
private:
  SolidMesh* plasticMesh;

  const unsigned pin_boundary = 3;

  double nu = 0.3;
  double E = 10;

  double dt = 1e-2;

  ConstitutiveLaw* constitutive_law_pt;
  PlasticConstitutiveLaw* plastic_constitutive_law_pt;

  TimeStepper* my_time_stepper_pt;

  // Doc info.
  DocInfo doc_info;

  void setup()
  {
    my_time_stepper_pt = new Newmark<2>;
    add_time_stepper_pt(my_time_stepper_pt);
    plasticMesh = new ElasticRectangularQuadMesh<ELEMENT>(10, 10, 1, 1, time_stepper_pt());

    oomph_info << "Created a simple quad mesh with " << plasticMesh->nelement()
               << " elements" << std::endl;


    // Create elastic constitutive law
    constitutive_law_pt = new GeneralisedHookean(&nu, &E);

    plastic_constitutive_law_pt = new PlasticConstitutiveLaw();

    Newton_solver_tolerance = 1e-6;
    Max_newton_iterations = 100;

    // Apply the constitutive laws to the elements
    for (unsigned e = 0; e < plasticMesh->nelement(); e++)
    {
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(plasticMesh->element_pt(e));

      // Assign constitutive law
      el_pt->constitutive_law_pt() = constitutive_law_pt;

      el_pt->plastic_constitutive_law_pt() = plastic_constitutive_law_pt;

      el_pt->assign_plastic_timestepper(my_time_stepper_pt);

      el_pt->body_force_fct_pt() = &BodyForce::body_force;
    }

    add_sub_mesh(plasticMesh);
    build_global_mesh();

    initialise_dt(dt);
    assign_initial_values_impulsive(dt);


    doc_info.set_directory("PlasticTestProblemOut");

    setup_initial_condition();

    // Do equation numbering
    oomph_info << "# of dofs " << assign_eqn_numbers() << std::endl;
  }

  void output()
  {
    int currentOutputNumber = currentOutputNumber = doc_info.number()++;
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
    const unsigned nnode = plasticMesh->nboundary_node(pin_boundary);
    for (unsigned l = 0; l < nnode; l++)
    {
      SolidNode* node_pt = dynamic_cast<SolidNode*>(plasticMesh->boundary_node_pt(pin_boundary, l));
      node_pt->pin_position(0);
      node_pt->pin_position(1);
      node_pt->pin_position(2);
    }
  }


  void solve()
  {
    output();
    unsigned nsolve = 0;
    while (time_pt()->time() <= 15.0)
    {
      // Output the current state
      if (nsolve++ % 10 == 0)
      {
        output();
      }

      unsteady_newton_solve(dt);

      Traction::pressure = time_pt()->time();
    }
  }

public:
  void run()
  {
    setup();
    solve();
  }
};

int main()
{
  PlasticTestProblem<QPlasticPVDElement<2, 2>> problem;
  // PlasticTestProblem<QPVDElement<2, 2>> problem;

  problem.run();
  return 0;
}