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

// Elements used to solve an auxillary problem for benchmarking
#include "benchmark_elements.h"

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
  // Poisson ratio - only used in the compressible case
  // (nu = 0.5 for incompressible material)
  double nu = 0.4;

  // The external pressure
  double ext_pressure = 0.05;

  // Traction function for appyling a pressure on the outer surface
  void external_pressure_fct(const double& time,
                             const Vector<double>& xi,
                             const Vector<double>& x,
                             const Vector<double>& n,
                             Vector<double>& traction)
  {
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

  /// If using a pressure formulation, is the solid incompressible?
  bool incompressible = false;

  // Number of elements in radial direction
  unsigned n_radial = 10;

  // Number of elements in azimuthal direction
  unsigned n_azimuthal = 10;

  // Inner radius
  double a = 1.0;

  // Thickness
  double h = 0.5;
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
}; // End mesh

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

private:
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

  /// Complete the problem setup (attach function/physical parameter pointers)
  void complete_problem_setup();

  /// Bulk mesh
  ElasticTwoDAnnularMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to the "surface" mesh
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

  // Create the mesh: Quarter-annulus in the quadrant r>0, z>0.
  Bulk_mesh_pt =
    new ElasticTwoDAnnularMesh<ELEMENT>(GlobalSimSettings::n_azimuthal,
                                        GlobalSimSettings::n_radial,
                                        GlobalSimSettings::a,
                                        GlobalSimSettings::h);

  // Create the "surface mesh" that will contain only the interface
  // elements. The constructor just creates the mesh without giving
  // it any elements, nodes, etc.
  Surface_mesh_pt = new Mesh;

  // Create the boundary elements
  create_boundary_elements();

  // Add the two submeshes
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_pt);

  // Combine all sub-meshes into a single mesh
  build_global_mesh();

  // Define a constitutive law from a strain energy density function.
  // Version depends on the incompressibility flag
  StrainEnergyFunction* strain_energy_function_pt;
  if (GlobalSimSettings::incompressible)
  {
    strain_energy_function_pt = new IncompressibleNeoHookean();
  }
  else
  {
    strain_energy_function_pt = new NeoHookean(&GlobalPhysicalVariables::nu);
  }
  Constitutive_law_pt =
    new IsotropicStrainEnergyFunctionConstitutiveLaw(strain_energy_function_pt);

  // Complete the problem set up
  complete_problem_setup();

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

    // Cast to element with pressure dofs
    QAxisymmetricCylindricalPVDWithPressureElement* el_with_pres_pt =
      dynamic_cast<QAxisymmetricCylindricalPVDWithPressureElement*>(el_pt);

    // Can set compressible/incompressible condition here
    if (el_with_pres_pt != 0)
    {
      if (GlobalSimSettings::incompressible)
      {
        el_with_pres_pt->set_incompressible();
      }
      else
      {
        el_with_pres_pt->set_compressible();
      }
    }
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
    el_pt->traction_fct_pt() = &GlobalPhysicalVariables::external_pressure_fct;

  } // End of loop over interface elements

  // Set the Dirichlet boundary conditions for this problem
  // Find number of nodes on boundary 3
  unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(3);
  // Loop over nodes on fixed boundary
  for (unsigned n = 0; n < n_boundary_node; n++)
  {
    // Pin radial position on symmetry axis only
    Bulk_mesh_pt->boundary_node_pt(3, n)->pin_position(0);
  }
  n_boundary_node = Bulk_mesh_pt->nboundary_node(1);
  // Loop over nodes on fixed boundary
  for (unsigned n = 0; n < n_boundary_node; n++)
  {
    // Pin radial position on symmetry axis only
    Bulk_mesh_pt->boundary_node_pt(1, n)->pin_position(1);
  }
} // End complete_problem_setup

//============start_of_create_boundary_elements==========================
/// Create boundary elements
//=======================================================================
template<class ELEMENT>
void ShellProblem<ELEMENT>::create_boundary_elements()
{
  // Loop over elements on boundary 2
  unsigned n_element = Bulk_mesh_pt->nboundary_element(2);
  for (unsigned e = 0; e < n_element; e++)
  {
    // Set a pointer to the bulk element we wish to our interface
    // element to
    ELEMENT* bulk_element_pt =
      dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(2, e));

    // Find the index of the face of element e along boundary b
    int face_index = Bulk_mesh_pt->face_index_at_boundary(2, e);

    // Create the interface element
    AxisymmetricCylindricalSolidTractionElement<ELEMENT>* interface_element_pt =
      new AxisymmetricCylindricalSolidTractionElement<ELEMENT>(bulk_element_pt,
                                                               face_index);

    // Add the interface element to the surface mesh
    this->Surface_mesh_pt->add_element_pt(interface_element_pt);
    interface_element_pt->set_boundary_number_in_bulk_mesh(2);
  }
} // end of create_boundary_elements


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
          "%s/bulk_soln_nu_%g_pres_%g.dat",
          Doc_info.directory().c_str(),
          GlobalPhysicalVariables::nu,
          GlobalPhysicalVariables::ext_pressure);
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file, npts);
  some_file.close();

  sprintf(filename,
          "%s/bulk_soln_nu_%g_pres_%g.vtu",
          Doc_info.directory().c_str(),
          GlobalPhysicalVariables::nu,
          GlobalPhysicalVariables::ext_pressure);
  some_file.open(filename);
  Bulk_mesh_pt->output_paraview(some_file, npts);
  some_file.close();

  // Increment the counter
  Doc_info.number()++;
}


//==============start_auxillary_problem===============================
/// Axisymmetric compression of a finite thickness shell in cylindrical
/// polar coordiantes
//====================================================================
class AuxiliaryShellProblem : public Problem
{
public:
  /// Constructor:
  AuxiliaryShellProblem();

  /// Destructor
  ~AuxiliaryShellProblem() {}

  /// Doc the solution
  void doc_solution();

  /// Solve the auxiliary problem
  void solve_auxiliary_problem();

private:
  /// Bulk mesh
  OneDMesh<ShellBenchmarkElement<3>>* Bulk_mesh_pt;

  /// Pointer to the "surface" mesh
  Mesh* Surface_mesh_pt;

  // DocInfo object
  DocInfo Doc_info;
};


//===============start_constructor========================================
/// Constructor for auxiliary 1D problem
//========================================================================
AuxiliaryShellProblem::AuxiliaryShellProblem()
{
  // Set output directory
  Doc_info.set_directory(GlobalSimSettings::result_folder);

  // Initialise counter for solutions
  Doc_info.number() = 0;

  // Create the mesh: use n_radial elements to match the 2D problem
  double outer_coord = GlobalSimSettings::a + GlobalSimSettings::h;
  Bulk_mesh_pt = new OneDMesh<ShellBenchmarkElement<3>>(
    GlobalSimSettings::n_radial, GlobalSimSettings::a, outer_coord);

  // Only complete the setup if the material is compressible. Otherwise,
  // an analytical solution will be used.
  if (!GlobalSimSettings::incompressible)
  {
    // Set the Poisson ratio in each element
    unsigned n_el = Bulk_mesh_pt->nelement();
    for (unsigned e = 0; e < n_el; e++)
    {
      // Cast to element
      ShellBenchmarkElement<3>* el_pt =
        dynamic_cast<ShellBenchmarkElement<3>*>(Bulk_mesh_pt->element_pt(e));

      // Set the Poisson ratio
      el_pt->nu_pt() = &GlobalPhysicalVariables::nu;
    }

    // Loop over nodes and set the initial guess in each node: r=R
    unsigned n_node = Bulk_mesh_pt->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      double r = Bulk_mesh_pt->node_pt(n)->x(0);
      Bulk_mesh_pt->node_pt(n)->set_value(0, r);
    }

    // Create the "surface mesh" for application of boundary pressure
    Surface_mesh_pt = new Mesh;

    // Cast to the right-most element in the mesh
    ShellBenchmarkElement<3>* bulk_element_pt =
      dynamic_cast<ShellBenchmarkElement<3>*>(
        Bulk_mesh_pt->element_pt(GlobalSimSettings::n_radial - 1));

    // Find the index of the face of element e (0) along boundary b (1)
    int face_index = 1; // Bulk_mesh_pt->face_index_at_boundary(1, 0);

    // Create the interface element
    ShellBenchmarkPressureElement* interface_element_pt =
      new ShellBenchmarkPressureElement(bulk_element_pt, face_index);

    // Add the interface element to the surface mesh
    this->Surface_mesh_pt->add_element_pt(interface_element_pt);
    interface_element_pt->set_boundary_number_in_bulk_mesh(1);

    // Set the pressure in the element
    interface_element_pt->pressure_pt() =
      &GlobalPhysicalVariables::ext_pressure;

    // Add the two submeshes
    add_sub_mesh(Bulk_mesh_pt);
    add_sub_mesh(Surface_mesh_pt);

    // Combine all sub-meshes into a single mesh
    build_global_mesh();

    // Assign eqn numbers
    oomph_info << "Number of equations for auxiliary problem: "
               << assign_eqn_numbers() << std::endl;
  }
} // end constructor

//========================================================================
/// Solve the auxiliary problem
//========================================================================
void AuxiliaryShellProblem::solve_auxiliary_problem()
{
  // If the material is incompressible, the problem can be solved
  // analytically, but yields an algebriac equation which must still be
  // solved numerically. This is done using a manual newton solve.
  if (GlobalSimSettings::incompressible)
  {
    // Storage for residual and Jacobian
    double residual, jacobian;

    // Containers for undeformed and deformed inner/outer radii.
    Vector<double> xi(2), r(2);
    xi[0] = GlobalSimSettings::a;
    xi[1] = GlobalSimSettings::a + GlobalSimSettings::h;

    // Set the initial guess
    r[0] = xi[0];
    r[1] = xi[1];

    // Derivative of outer coord. w.r.t. inner coord.
    double dr1_dr0;

    // Poisson ratio in incompressible case
    const double nu = 0.5;

    // Residual from the initial guess
    r[1] = pow(pow(r[0], 3.0) + pow(xi[1], 3.0) - pow(xi[0], 3.0), 1.0 / 3.0);
    residual =
      (1.0 / (1.0 + nu)) * (0.25 * pow(xi[0] / r[0], 4.0) + xi[0] / r[0] -
                            0.25 * pow(xi[1] / r[1], 4.0) - xi[1] / r[1]) -
      GlobalPhysicalVariables::ext_pressure;

    // Newton iteration loop
    while (std::fabs(residual) > 1.0e-8)
    {
      // Calculate the Jacobian
      dr1_dr0 =
        pow(pow(r[0], 3.0) + pow(xi[1], 3.0) - pow(xi[0], 3.0), -2.0 / 3.0) *
        r[0] * r[0];
      jacobian = (1.0 / (1.0 + nu)) *
                 (-1.0 * (pow(xi[0] / r[0], 4.0) + xi[0] / r[0]) / r[0] +
                  (pow(xi[1] / r[1], 4.0) + xi[1] / r[1]) * dr1_dr0 / r[1]);
      // Take a Newton step
      r[0] = r[0] - residual / jacobian;

      // Calculate the new residual
      r[1] = pow(pow(r[0], 3.0) + pow(xi[1], 3.0) - pow(xi[0], 3.0), 1.0 / 3.0);
      residual =
        (1 / (1.0 + nu)) * (0.25 * pow(xi[0] / r[0], 4.0) + xi[0] / r[0] -
                            0.25 * pow(xi[1] / r[1], 4.0) - xi[1] / r[1]) -
        GlobalPhysicalVariables::ext_pressure;
    }

    // With the solution obtained, set nodal values from the analytic solution
    unsigned n_node = Bulk_mesh_pt->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      double xi_tmp = Bulk_mesh_pt->node_pt(n)->x(0);
      Bulk_mesh_pt->node_pt(n)->set_value(
        0, pow(pow(r[0], 3.0) + pow(xi_tmp, 3.0) - pow(xi[0], 3.0), 1.0 / 3.0));
    }
  }
  // Otherwise, perform a standard Newton solve using the benchmark elements.
  else
  {
    newton_solve();
  }
}

//========================================================================
/// Doc the solution
//========================================================================
void AuxiliaryShellProblem::doc_solution()
{
  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts = 3;

  sprintf(filename,
          "%s/aux_soln_nu_%g_pres_%g.dat",
          Doc_info.directory().c_str(),
          GlobalPhysicalVariables::nu,
          GlobalPhysicalVariables::ext_pressure);
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file, npts);
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

  string str_temp_solid_pressure;
  CommandLineArgs::specify_command_line_flag(
    "--use_solid_pressure", &str_temp_solid_pressure, "Use solid pressure");

  string str_temp_incompressible;
  CommandLineArgs::specify_command_line_flag(
    "--incompressible", &str_temp_incompressible, "Incompressible");

  CommandLineArgs::specify_command_line_flag(
    "--poisson", &GlobalPhysicalVariables::nu, "Poisson");

  CommandLineArgs::specify_command_line_flag(
    "--pressure", &GlobalPhysicalVariables::ext_pressure, "External pressure");

  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  istringstream(str_temp_solid_pressure) >> std::boolalpha >>
    GlobalSimSettings::use_pressure_formulation;
  istringstream(str_temp_incompressible) >> std::boolalpha >>
    GlobalSimSettings::incompressible;

  // If the Poisson ratio has been set to 0.5, but the
  // use_pressure_formulation and/or incompressible flags are false, set them to
  // true
  if (std::fabs(GlobalPhysicalVariables::nu - 0.5) < 1.0e-10 &&
      !(GlobalSimSettings::incompressible &&
        GlobalSimSettings::use_pressure_formulation))
  {
    std::string warning_message =
      "Inconsistency: Poisson ratio has been set to 0.5 but the  \n"
      "use_pressure_formulation and/or incompressible flags are false. \n"
      "Setting both of these to true to proceed.";
    OomphLibWarning(
      warning_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    GlobalSimSettings::use_pressure_formulation = true;
    GlobalSimSettings::incompressible = true;
  }

  // If the incompressibility flag has been set, we must use the solid
  // pressure formulation. The Poisson ratio will not be used in this case.
  if (GlobalSimSettings::incompressible &&
      !GlobalSimSettings::use_pressure_formulation)
  {
    std::string warning_message =
      "Inconsistency: The incompressibility flag is true but we are not using "
      "\n"
      "the solid pressure formulation. Setting use_solid_pressure to true to "
      "proceed. \n";
    OomphLibWarning(
      warning_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    GlobalSimSettings::use_pressure_formulation = true;
  }

  // Set the Poisson ratio to 0.5 if incompressibility has been specified. This
  // is solely for the purpose of output
  if (GlobalSimSettings::incompressible)
  {
    GlobalPhysicalVariables::nu = 0.5;
  }

  // End consistency checks

  // Create the problem
  if (GlobalSimSettings::use_pressure_formulation)
  {
    // Build the problem with an additional degree of freedom for pressure
    ShellProblem<QAxisymmetricCylindricalPVDWithPressureElement> problem;

    problem.newton_solve();
    problem.doc_solution();
  }
  else
  {
    // Standard build
    ShellProblem<QAxisymmetricCylindricalPVDElement<3>> problem;

    problem.newton_solve();
    problem.doc_solution();
  }

  // Build and solve the auxiliary problem
  AuxiliaryShellProblem aux_problem;
  aux_problem.solve_auxiliary_problem();
  aux_problem.doc_solution();
} // end main
