// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2026 Matthias Heil and Andrew Hazel
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

// Driver code for benchmark of axisym non-linear solid equations
// in cylindrical coordinates: Thick spherical shell with external pressure
// Compare against analytical solution for an incompressible material

// Generic routines
#include "generic.h"

// Axisymmetric non-linear solid equations in cylindrical coords
#include "axisym_solid.h"

// Constitutive laws
#include "constitutive.h"

// The mesh: annular
#include "meshes/annular_mesh.h"
// Line mesh for auxiliary 1D problem
#include "meshes/one_d_mesh.h"

using namespace std;
using namespace oomph;

//=======start_namespace==========================================
/// Global physical variables
//================================================================
namespace GlobalPhysicalVariables
{
  // Poisson ratio
  double nu = 0.4;

  // Coefficient amplitudes for the manufactured solution
  double A0 = 0.05;
  double B0 = 0.05;

  // Oscillation frequency
  double omega = 5.0;

  // Get the principal stretches from the manufactured solution
  Vector<double> get_principal_stretches(const double& time, const double& r)
  {
    Vector<double> lambda(2);
    double A = A0 * sin(omega * time);
    double B = B0 * sin(omega * time);
    lambda[0] = 1.0 + A - 2.0 * B * pow(r, -3.0);
    lambda[1] = 1.0 + A + B * pow(r, -3.0);
    return lambda;
  }

  // Get derivatives of prinicipal stretches w.r.t. r
  Vector<double> get_principal_stretch_derivatives(const double& time,
                                                   const double& r)
  {
    Vector<double> dlambda_dr(2);
    double B = B0 * sin(omega * time);
    dlambda_dr[0] = 6.0 * B / pow(r, 4.0);
    dlambda_dr[1] = -3.0 * B / pow(r, 4.0);
    return dlambda_dr;
  }

  // Evaluate the deformed position
  void deformed_position_fct(const double& time,
                             const Vector<double>& xi,
                             Vector<double>& x)
  {
    double r = pow(xi[0] * xi[0] + xi[1] * xi[1], 0.5);
    double theta = atan2(xi[1], xi[0]);
    x[0] =
      xi[0] + (A0 * r + B0 * pow(r, -2.0)) * sin(omega * time) * cos(theta);
    x[1] =
      xi[1] + (A0 * r + B0 * pow(r, -2.0)) * sin(omega * time) * sin(theta);
  }

  // Evaluate the deformed position
  void deformed_velocity_fct(const double& time,
                             const Vector<double>& xi,
                             Vector<double>& x_dot)
  {
    double r = pow(xi[0] * xi[0] + xi[1] * xi[1], 0.5);
    double theta = atan2(xi[1], xi[0]);
    x_dot[0] =
      omega * (A0 * r + B0 * pow(r, -2.0)) * cos(omega * time) * cos(theta);
    x_dot[1] =
      omega * (A0 * r + B0 * pow(r, -2.0)) * cos(omega * time) * sin(theta);
  }

  // Evaluate the deformed position
  void deformed_accel_fct(const double& time,
                          const Vector<double>& xi,
                          Vector<double>& x_ddot)
  {
    double r = pow(xi[0] * xi[0] + xi[1] * xi[1], 0.5);
    double theta = atan2(xi[1], xi[0]);
    x_ddot[0] = -1.0 * omega * omega * (A0 * r + B0 * pow(r, -2.0)) *
                sin(omega * time) * cos(theta);
    x_ddot[1] = -1.0 * omega * omega * (A0 * r + B0 * pow(r, -2.0)) *
                sin(omega * time) * sin(theta);
  }

  // Body force function
  void body_force_fct(const double& time,
                      const Vector<double>& xi,
                      Vector<double>& b)
  {
    // Calculate the radial coordinate
    double r = pow(xi[0] * xi[0] + xi[1] * xi[1], 0.5);
    // Radial acceleration
    double accel =
      -1.0 * omega * omega * (A0 * r + B0 * pow(r, -2.0)) * sin(omega * time);
    // Principal stretches
    Vector<double> lambda = get_principal_stretches(time, r);
    // Jacobian
    double J = lambda[0] * lambda[1] * lambda[1];
    // Derivative of radial principal stretch and Jacobian w.r.t. r
    Vector<double> dlambda_dr = get_principal_stretch_derivatives(time, r);
    double dJdr = dlambda_dr[0] * lambda[1] * lambda[1] +
                  2.0 * lambda[0] * lambda[1] * dlambda_dr[1];

    // Dimensionless Lame parameters
    double mu_L = 1.0 / (2.0 * (1.0 + nu));
    double lambda_L = nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

    // Evaluate the body force
    double b_spher =
      accel - mu_L * (dlambda_dr[0] + (2.0 / r) * (lambda[0] - lambda[1])) -
      lambda_L * dJdr / (J * lambda[0]) -
      (dlambda_dr[0] / (lambda[0] * lambda[0]) +
       (2.0 / r) * (1.0 / lambda[1] - 1.0 / lambda[0])) *
        (mu_L - lambda_L * log(J));

    double theta = atan2(xi[1], xi[0]);
    b[0] = b_spher * cos(theta);
    b[1] = b_spher * sin(theta);
  }

  // External traction force function
  void traction_force_fct(const double& time,
                          const Vector<double>& xi,
                          const Vector<double>& x,
                          const Vector<double>& n,
                          Vector<double>& traction)
  {
    // Calculate the (undeformed) radial coordinate
    double r = pow(xi[0] * xi[0] + xi[1] * xi[1], 0.5);
    // Principal stretches
    Vector<double> lambda = get_principal_stretches(time, r);
    // Jacobian
    double J = lambda[0] * lambda[1] * lambda[1];

    // Dimensionless Lame parameters
    double mu_L = 1.0 / (2.0 * (1.0 + nu));
    double lambda_L = nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

    // Pressure required for the manufactured solution to work
    double ext_pressure =
      -1.0 * (mu_L * (lambda[0] * lambda[0] - 1.0) + lambda_L * log(J)) / J;

    // Return the traction
    traction[0] = -1.0 * ext_pressure * n[0];
    traction[1] = -1.0 * ext_pressure * n[1];
  }
} // namespace GlobalPhysicalVariables


//=======start_namespace==========================================
/// Global simulation settings
//================================================================
namespace GlobalSimSettings
{
  // Folder for output
  string result_folder = "RESLT";

  /// Flag for using a pressure formulation
  bool use_pressure_formulation = false;

  // Number of elements in radial direction
  unsigned n_radial = 14;

  // Number of elements in azimuthal direction
  unsigned n_azimuthal = 14;

  // Inner radius
  double a = 1.0;

  // Thickness
  double h = 0.5;

  // Number of timesteps
  unsigned Nsteps = 120;

  // Size of the timestep
  double dt;
} // namespace GlobalSimSettings

//======================start_mesh================================
/// Upgrade the mesh to a solidmesh.
/// Quarter-annulus in quadrant r>0, z>0
/// ntheta - no. elements in theta direction
/// nr - no. elements in the radial direction
/// a - Inner radius of annulus
/// h - thickness of annulus
//================================================================
template<class ELEMENT>
class ElasticTwoDAnnularMesh : public virtual TwoDAnnularMesh<ELEMENT>,
                               public virtual SolidMesh
{
public:
  /// Constructor: Build mesh and copy Eulerian coords to Lagrangian
  /// ones so that the initial configuration is the stress-free one.
  ElasticTwoDAnnularMesh<ELEMENT>(
    const unsigned& ntheta,
    const unsigned& nr,
    const double& a,
    const double& h,
    TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
    : RectangularQuadMesh<ELEMENT>(
        ntheta, nr, 1.0, 1.0, false, time_stepper_pt),
      TwoDAnnularMesh<ELEMENT>(false,
                               0.25,
                               ntheta,
                               nr,
                               a,
                               h,
                               MathematicalConstants::Pi,
                               time_stepper_pt)
  {
    // Make the current configuration the undeformed one by
    // setting the nodal Lagrangian coordinates to their current
    // Eulerian ones
    set_lagrangian_nodal_coordinates();
  }
};

//==============start_problem=========================================
/// Axisymmetric compression of a finite thickness shell in cylindrical
/// polar coordiantes
//====================================================================
template<class ELEMENT>
class ShellProblem : public Problem
{
public:
  /// Constructor:
  ShellProblem();

  /// Destructor
  ~ShellProblem() {}

  /// Doc the solution
  void doc_solution();

  /// Set the initial condition
  void set_initial_condition();

  /// Run through the timesteps
  void unsteady_run();

private:
  /// Complete the problem setup (attach function/physical parameter pointers)
  void complete_problem_setup();

  /// Create the boundary elements
  void create_boundary_elements();

  /// Delete boundary elements and wipe the associated mesh
  void delete_boundary_elements()
  {
    // How many surface elements are in the surface mesh
    unsigned n_element = Surface_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Kill surface element
      delete Surface_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Surface_mesh_pt->flush_element_and_node_storage();
  } // end of delete_boundary_elements

  /// Bulk mesh
  ElasticTwoDAnnularMesh<ELEMENT>* Bulk_mesh_pt;

  /// Surface mesh for traction force
  Mesh* Surface_mesh_pt;

  /// Pointer to constitutive law
  ConstitutiveLaw* Constitutive_law_pt;

  // DocInfo object
  DocInfo Doc_info;
};

//===============start_constructor========================================
/// Constructor for unstructured solid problem
//========================================================================
template<class ELEMENT>
ShellProblem<ELEMENT>::ShellProblem()
{
  // Set output directory
  Doc_info.set_directory(GlobalSimSettings::result_folder);

  // Initialise counter for solutions
  Doc_info.number() = 0;

  // Allocate the timestepper
  add_time_stepper_pt(new Newmark<1>);

  // Calculate the timestep size
  GlobalSimSettings::dt =
    2.0 * MathematicalConstants::Pi /
    (GlobalPhysicalVariables::omega * GlobalSimSettings::Nsteps);

  // Create the mesh: Quarter-annulus in the quadrant r>0, z>0.
  Bulk_mesh_pt =
    new ElasticTwoDAnnularMesh<ELEMENT>(GlobalSimSettings::n_azimuthal,
                                        GlobalSimSettings::n_radial,
                                        GlobalSimSettings::a,
                                        GlobalSimSettings::h,
                                        time_stepper_pt());

  // Create the surface mesh
  Surface_mesh_pt = new Mesh;

  // Create the boundary elements
  create_boundary_elements();

  // Add the two submeshes
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_pt);

  // Combine all sub-meshes into a single mesh
  build_global_mesh();

  // Define a constitutive law from a strain energy density function.
  // Compressible Neo-Hookean
  StrainEnergyFunction* strain_energy_function_pt =
    new NeoHookean(&GlobalPhysicalVariables::nu);
  Constitutive_law_pt =
    new IsotropicStrainEnergyFunctionConstitutiveLaw(strain_energy_function_pt);

  // Complete the problem set up
  complete_problem_setup();

  // Set the initial condition consistently
  set_initial_condition();

  // Assign eqn numbers
  oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
} // end constructor

//========start_of_complete_problem_setup==================================
/// Assign function/parameter pointers, set boundary conditions etc.
//=========================================================================
template<class ELEMENT>
void ShellProblem<ELEMENT>::complete_problem_setup()
{
  // Determine number of bulk elements in mesh
  const unsigned n_element_bulk = Bulk_mesh_pt->nelement();

  // Loop over the bulk elements
  for (unsigned e = 0; e < n_element_bulk; e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // Set the constitutive law
    el_pt->constitutive_law_pt() = Constitutive_law_pt;

    // Set the traction function for applying pressure
    el_pt->body_force_fct_pt() = &GlobalPhysicalVariables::body_force_fct;
  } // End of loop over bulk elements

  // Determine number of 1D interface elements in mesh
  const unsigned n_interface_element = Surface_mesh_pt->nelement();

  // Loop over the interface elements
  for (unsigned e = 0; e < n_interface_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    AxisymmetricCylindricalSolidTractionElement<ELEMENT>* el_pt =
      dynamic_cast<AxisymmetricCylindricalSolidTractionElement<ELEMENT>*>(
        Surface_mesh_pt->element_pt(e));

    // Set the traction function for applying pressure
    el_pt->traction_fct_pt() = &GlobalPhysicalVariables::traction_force_fct;

  } // End of loop over interface elements

  // Set Dirichlet boundary conditions
  unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(1);
  for (unsigned n = 0; n < n_boundary_node; n++)
  {
    // Pin radial position on symmetry axis only
    Bulk_mesh_pt->boundary_node_pt(1, n)->pin_position(1);
  }

  n_boundary_node = Bulk_mesh_pt->nboundary_node(3);
  // Loop over nodes on fixed boundary
  for (unsigned n = 0; n < n_boundary_node; n++)
  {
    // Pin radial position on symmetry axis only
    Bulk_mesh_pt->boundary_node_pt(3, n)->pin_position(0);
  }
} // end complete_problem_setup

//============start_of_create_boundary_elements==========================
/// Create boundary elements
//=======================================================================
template<class ELEMENT>
void ShellProblem<ELEMENT>::create_boundary_elements()
{
  // Loop over boundaries 0 and 2 (inner and outer surface of shell)
  for (unsigned ibound : {0, 2})
  {
    // Loop over elements on boundaries 0 and 2
    unsigned n_element = Bulk_mesh_pt->nboundary_element(ibound);
    for (unsigned e = 0; e < n_element; e++)
    {
      // Set a pointer to the bulk element we wish to our interface
      // element to
      ELEMENT* bulk_element_pt =
        dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(ibound, e));

      // Find the index of the face of element e along boundary b
      int face_index = Bulk_mesh_pt->face_index_at_boundary(ibound, e);

      // Create the interface element
      AxisymmetricCylindricalSolidTractionElement<ELEMENT>*
        interface_element_pt =
          new AxisymmetricCylindricalSolidTractionElement<ELEMENT>(
            bulk_element_pt, face_index);

      // Add the interface element to the surface mesh
      this->Surface_mesh_pt->add_element_pt(interface_element_pt);
      interface_element_pt->set_boundary_number_in_bulk_mesh(ibound);
    }
  }
} // end of create_boundary_elements

//============start_of_set_initial_condition==============================
/// Given the form of the manufactured solution, this really only needs
/// doing for the velocity
//========================================================================
template<class ELEMENT>
void ShellProblem<ELEMENT>::set_initial_condition()
{
  // First initialize dt and set the weights
  initialise_dt(GlobalSimSettings::dt);

  // Loop over nodes
  unsigned n_node = Bulk_mesh_pt->nnode();
  for (unsigned n = 0; n < n_node; n++)
  {
    // Container for position/velocity/acceleration
    Vector<double> value(2);

    // Cast to (solid)node
    SolidNode* nod_pt = dynamic_cast<SolidNode*>(Bulk_mesh_pt->node_pt(n));

    // Get the undeformed position
    Vector<double> xi(2), x(2);
    xi[0] = nod_pt->lagrangian_position(0);
    xi[1] = nod_pt->lagrangian_position(1);

    // Set values of deformed position at the present and previous times
    GlobalPhysicalVariables::deformed_position_fct(0.0, xi, value);
    nod_pt->x(0, 0) = value[0];
    nod_pt->x(0, 1) = value[1];
    GlobalPhysicalVariables::deformed_position_fct(
      -1.0 * GlobalSimSettings::dt, xi, value);
    nod_pt->x(1, 0) = value[0];
    nod_pt->x(1, 1) = value[1];

    // Set the velocity at the previous time
    GlobalPhysicalVariables::deformed_velocity_fct(
      -1.0 * GlobalSimSettings::dt, xi, value);
    nod_pt->x(2, 0) = value[0];
    nod_pt->x(2, 1) = value[1];

    // Set the acceleration at the previous time
    GlobalPhysicalVariables::deformed_accel_fct(
      -1.0 * GlobalSimSettings::dt, xi, value);
    nod_pt->x(3, 0) = value[0];
    nod_pt->x(3, 1) = value[1];
  }
} // end set_initial_condition

//========================================================================
/// Run through the timesteps
//========================================================================
template<class ELEMENT>
void ShellProblem<ELEMENT>::unsteady_run()
{
  // Loop over timesteps
  for (unsigned step = 0; step < GlobalSimSettings::Nsteps; step++)
  {
    // Perform the timestep
    unsteady_newton_solve(GlobalSimSettings::dt);

    // Doc the solution
    doc_solution();
  }
}

//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ShellProblem<ELEMENT>::doc_solution()
{
  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts = 3;

  sprintf(filename,
          "%s/bulk_soln%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file, npts);
  some_file.close();

  sprintf(filename,
          "%s/bulk_soln%i.vtu",
          Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output_paraview(some_file, npts);
  some_file.close();

  // Increment the counter
  Doc_info.number()++;
}

//===========start_main===================================================
/// Demonstrate how to solve an unstructured solid problem
//========================================================================
int main(int argc, char* argv[])
{
  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  string str_temp;
  CommandLineArgs::specify_command_line_flag(
    "--use_solid_pressure", &str_temp, "Use solid pressure");

  CommandLineArgs::specify_command_line_flag(
    "--poisson", &GlobalPhysicalVariables::nu, "Poisson");

  CommandLineArgs::specify_command_line_flag(
    "--Nsteps", &GlobalSimSettings::Nsteps, "Number of Steps");

  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  istringstream(str_temp) >> std::boolalpha >>
    GlobalSimSettings::use_pressure_formulation;

  // Create the problem
  if (GlobalSimSettings::use_pressure_formulation)
  {
    // Build the problem with an additional degree of freedom for pressure
    ShellProblem<QAxisymmetricCylindricalPVDWithPressureElement> problem;

    // Perform time integration and docs
    problem.unsteady_run();
  }
  else
  {
    // Standard build
    ShellProblem<QAxisymmetricCylindricalPVDElement<3>> problem;

    // Perform time integration and docs
    problem.unsteady_run();
  }
} // end main
