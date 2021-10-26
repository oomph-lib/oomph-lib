//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
// Generic oomph-lib stuff
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// Include the headers required to generate the space-time mesh
#include "meshes/rectangle_with_moving_cylinder_mesh.h"
#include "meshes/extruded_cube_mesh_from_quad_mesh_with_macro_elements.h"

using namespace oomph;

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//======start_of_OscillatingCylinder_class================================
/// Oscillating cylinder class
//========================================================================
class OscillatingCylinder : public GeomObject
{
public:
  /// Constructor: Pass in the radius, the amplitude of the cylinder
  /// motion, the simulation Strouhal number and a pointer to time object.
  OscillatingCylinder(double* radius_pt,
                      double* amplitude_pt,
                      Time* time_pt) :
    GeomObject(1,2),
    Radius_pt(radius_pt),
    Amplitude_pt(amplitude_pt),
    Time_pt(time_pt)
  {}

  /// Destructor: Empty
  virtual ~OscillatingCylinder() {}

  /// Access function for the amplitude (lvalue)
  double& amplitude()
  {
    // Return the value of the amplitude
    return *Amplitude_pt;
  } // End of amplitude

  /// Access function for the Time pointer
  Time* time_pt()
  {
    // Return the Time pointer
    return Time_pt;
  } // End of time_pt

  /// Current position vector to material point at Lagrangian
  /// coordinate xi (steady version)
  void position(const Vector<double>& xi,
                Vector<double>& r) const
  {
    // X-coordinate
    r[0]=(*Radius_pt)*cos(xi[0]);

    // Y-coordinate
    r[1]=(*Radius_pt)*sin(xi[0]);
  } // End of position

  /// Current position vector to material point at Lagrangian
  /// coordinate xi (unsteady version). Implementation includes a
  /// transition phase where the cylinder oscillates to a smaller
  /// amplitude than the target value. Used to ensure that the solution
  /// isn't drastically different to that at the next time step. This
  /// can be disabled by setting Use_transition_phase to false.
  void position(const unsigned& t,
                const Vector<double>& xi,
                Vector<double>& r) const
  {
    // We can only use this if t=0
    if (t!=0)
    {
      // Create an output stream
      std::ostringstream error_message_stream;

      // Create an error message
      error_message_stream << "Trying to access the position at time history "
                           << "value, t: " << t << ".\nThis doesn't make sense "
                           << "in the space-time mesh!" << std::endl;

      // Throw the error message
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Calculate the coordinate before translation
    position(xi,r);

    // Get current time
    double time=Time_pt->time(t);

    // Scaling factor
    double arg=2.0*MathematicalConstants::Pi;

    // Calculate the translation
    double translation=(*Amplitude_pt)*sin(arg*time);

    // Update the y-coordinate
    r[1]+=translation;
  } // End of position

  /// Parametrised position on object: r(zeta). Evaluated at
  /// the continuous time value, t.
  virtual void position(const double& t,
                        const Vector<double>& xi,
                        Vector<double>& r) const
  {
    // Calculate the coordinate before translation
    position(xi,r);

    // Scaling factor
    double arg=2.0*MathematicalConstants::Pi;

    // Calculate the translation
    double translation=(*Amplitude_pt)*sin(arg*t);

    // Update the y-coordinate
    r[1]+=translation;
  } // End of position

  /// Velocity at any given point on the rigid cylinder at time, t
  virtual void velocity(const double& t, Vector<double>& u) const
  {
    // Scaling factor
    double arg=2.0*MathematicalConstants::Pi;

    // Zero velocity component in the x-direction
    u[0]=0.0;

    // Calculate the (non-zero) velocity component in the y-direction
    u[1]=arg*(*Amplitude_pt)*cos(arg*t);
  } // End of velocity

private:

  /// Radius of the cylinder
  double* Radius_pt;

  /// Non-dimensionalised amplitude of the cylinder motion
  double* Amplitude_pt;

  /// Pointer to the current time in the problem
  Time* Time_pt;
}; // End of OscillatingCylinder class

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//===== start_of_namespace=============================================
/// Namespace for exact solution for Poisson equation with "sharp step"
//=====================================================================
namespace TanhSolnForPoisson
{
  /// Parameter for steepness of step
  double Alpha=1.0;

  /// Orientation of non-normalised x-component of vector in plane direction
  double N_x=-1.0;

  /// Orientation of non-normalised y-component of vector in plane direction
  double N_y=-1.0;

  /// Orientation of non-normalised z-component of vector in plane direction
  double N_z=1.0;

  /// Orientation (x-coordinate of step plane)
  double X_0=0.0;

  /// Orientation (y-coordinate of step plane)
  double Y_0=0.0;

  /// Orientation (z-coordinate of step plane)
  double Z_0=0.0;

  /// Exact solution as a Vector
  void get_exact_u(const Vector<double>& x, Vector<double>& u)
  {
    double d=std::sqrt(N_x*N_x+N_y*N_y+N_z*N_z);
    u[0]=std::tanh(Alpha*((x[0]-X_0)*N_x+
                          (x[1]-Y_0)*N_y+
                          (x[2]-Z_0)*N_z)/d);
  } // End of get_exact_u

  /// Exact solution as a scalar
  void get_exact_u(const Vector<double>& x, double& u)
  {
    double d=std::sqrt(N_x*N_x+N_y*N_y+N_z*N_z);
    u=std::tanh(Alpha*((x[0]-X_0)*N_x+
                       (x[1]-Y_0)*N_y+
                       (x[2]-Z_0)*N_z)/d);
  } // End of get_exact_u

  /// Source function to make it an exact solution
  void get_source(const Vector<double>& x, double& source)
  {
    // Storage for the solution
    double u=0.0;

    // Calculate the solution
    get_exact_u(x,u);

    // Calculate the source value
    source=-2.0*u*(1.0-pow(u,2.0))*Alpha*Alpha;
  } // End of get_source
} // End of TanhSolnForPoisson namespace

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//======start_of_GlobalParameters_namespace===============================
/// Global parameters for the problem
//========================================================================
namespace GlobalParameters
{
  /// -------------------------Cylinder Properties-------------------------
  /// Pointer to the cylinder
  OscillatingCylinder* Cylinder_pt=0;

  /// Radius of the cylinder
  double Radius=0.5;

  /// Amplitude of the cylinder motion
  double Amplitude=0.50;
  /// -------------------------Cylinder Properties-------------------------

  /// -------------------------Domain Properties---------------------------
  /// Length of square central box domain
  double Length_of_box=10.0;

  /// The length of the extruded mesh in the z-direction
  double Length_z=1.0;

  /// Number of elements in the z-direction (in the extruded mesh)
  unsigned N_element_z=10;

  /// Number of uniform refinements before any solve
  unsigned N_uniform_refinement_before_solve=2;

  /// The radius of the annular region surrounding the cylinder
  double Annular_region_radius=1.0;

  /// Update parameters
  void update_parameters()
  {
    // Update the radius of the annular mesh region surrounding the cylinder
    Annular_region_radius=
      Radius+std::min(2.0*Radius,0.5*((0.5*Length_of_box)-Radius));
  } // End of update_parameters
  /// -------------------------Domain Properties-------------------------

  /// -----------------------Documentation Helpers-----------------------
  /// Doc info object
  DocInfo Doc_info;
  /// -----------------------Documentation Helpers-----------------------
} // End of GlobalParameters namespace

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

//=====start_of_ExtrudedMovingCylinderProblem=========================
/// Poisson problem in an extruded 3D domain. Used to illustrate the mesh
/// extrusion machinery with objects that move "over time" through their
/// macro-element representation. We're using a steady Poisson problem to
/// show that we don't need genuine time dependence to make use of the
/// extrusion toolset.
//========================================================================
template<class TWO_D_ELEMENT,class THREE_D_ELEMENT>
class ExtrudedMovingCylinderProblem : public Problem
{
public:
  /// Constructor: pass a pointer to the source function
  ExtrudedMovingCylinderProblem(
    const PoissonEquations<3>::PoissonSourceFctPt& source_fct_pt);

  /// Destructor: Empty
  ~ExtrudedMovingCylinderProblem() {}

  /// Update the problem specs before solve
  void actions_before_newton_solve() {}

  /// Update the problem specs after solve (empty)
  void actions_after_newton_solve()
  {
    /// Apply boundary conditions
    apply_boundary_conditions();
  }

  /// Update the problem specs after adaptation
  void actions_after_adapt()
  {
    /// Apply boundary conditions
    apply_boundary_conditions();

    // Complete the build of the elements; pass pointers to physical variables
    complete_element_setup();
  } // End of actions_after_adapt

  /// Create the 2D mesh and extrude it to create the 3D mesh
  void create_extruded_mesh();

  /// Apply boundary conditions
  void apply_boundary_conditions();

  // Complete the build of the elements; pass pointers to physical variables
  void complete_element_setup();

  /// Document the solution
  void doc_solution();

private:
  /// Pointer to the extruded messh
  RefineableExtrudedCubeMeshFromQuadMesh<THREE_D_ELEMENT>* Bulk_mesh_pt;

  /// Pointer to source function
  PoissonEquations<3>::PoissonSourceFctPt Source_fct_pt;
};

//=====start_of_ExtrudedMovingCylinderProblem=========================
/// Constructor
//========================================================================
template<class TWO_D_ELEMENT,class THREE_D_ELEMENT>
ExtrudedMovingCylinderProblem<TWO_D_ELEMENT,THREE_D_ELEMENT>::
ExtrudedMovingCylinderProblem(
  const PoissonEquations<3>::PoissonSourceFctPt& source_fct_pt) :
  Bulk_mesh_pt(0),
  Source_fct_pt(0)
{
  // Setup the steepness of step in the exact tanh solution
  TanhSolnForPoisson::Alpha=1.0;

  // Generate the 2D mesh and extrude it to create the 3D mesh we want
  create_extruded_mesh();

  // Apply the chosen Dirichlet boundary conditions
  apply_boundary_conditions();

  // Complete the build of the elements; pass pointers to physical variables
  complete_element_setup();

  // Attach the boundary conditions to the mesh
  oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
} // End of ExtrudedMovingCylinderProblem


//=====start_of_create_extruded_mesh======================================
/// Build the extruded mesh
//========================================================================
template<class TWO_D_ELEMENT,class THREE_D_ELEMENT>
void ExtrudedMovingCylinderProblem<TWO_D_ELEMENT,THREE_D_ELEMENT>::
create_extruded_mesh()
{
  // Start the clock
  double timer_s=TimingHelpers::timer();

  //-------------------
  // Create the 2D mesh
  //-------------------
  // Use BDF2
  add_time_stepper_pt(new Steady<0>);

  // Make a new cylinder
  // NOTE: This needs to be kept alive as long as updates are going to be
  // made to the mesh as it will be needed to calculate the positions of
  // nodes through their extruded macro-element local coordinates
  GlobalParameters::Cylinder_pt=
    new OscillatingCylinder(&GlobalParameters::Radius,
                            &GlobalParameters::Amplitude,
                            time_pt());

  // Build central mesh
  RefineableRectangleWithHoleAndAnnularRegionMesh<TWO_D_ELEMENT>* two_d_mesh_pt=
    new RefineableRectangleWithHoleAndAnnularRegionMesh<TWO_D_ELEMENT>(
    GlobalParameters::Cylinder_pt,
    GlobalParameters::Annular_region_radius,
    GlobalParameters::Length_of_box,
    time_stepper_pt(0));

  /// Number of uniform refinements before any solve
  for (unsigned i=0; i<GlobalParameters::N_uniform_refinement_before_solve; i++)
  {
    // Refine the pre-extruded mesh
    two_d_mesh_pt->refine_uniformly();
  }

  //--------------------------------
  // Now generate the extruded mesh!
  //--------------------------------
  // The number of elements in the z-direction
  unsigned n_z=GlobalParameters::N_element_z;

  // The length of the mesh in the z-direction
  double l_z=GlobalParameters::Length_z;

  // Create the extruded mesh
  Bulk_mesh_pt=new RefineableExtrudedCubeMeshFromQuadMesh<THREE_D_ELEMENT>(
    two_d_mesh_pt,n_z,l_z);

  // Set the error estimator
  Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

  // Set limits on the error
  Bulk_mesh_pt->max_permitted_error()=1.0e-02;
  Bulk_mesh_pt->min_permitted_error()=1.0e-04;

  // Store the mesh pointer
  Problem::mesh_pt()=Bulk_mesh_pt;

  // Document the setup time
  oomph_info << "\nTime taken for mesh extrusion/setup [sec]: "
             << TimingHelpers::timer()-timer_s << std::endl;
} // End of create_extruded_mesh


//=====start_of_apply_boundary_conditions=================================
/// Apply the chosen Dirichlet boundary conditions
//========================================================================
template<class TWO_D_ELEMENT,class THREE_D_ELEMENT>
void ExtrudedMovingCylinderProblem<TWO_D_ELEMENT,THREE_D_ELEMENT>::
apply_boundary_conditions()
{
  // How many boundaries does this mesh have?
  unsigned n_boundary=mesh_pt()->nboundary();

  // Loop over the boundaries
  for (unsigned b=0; b<n_boundary; b++)
  {
    // How many nodes are there on this boundary?
    unsigned n_boundary_node=mesh_pt()->nboundary_node(b);

    // Loop over the nodes on boundary
    for (unsigned n=0; n<n_boundary_node; n++)
    {
      // Grab a pointer to the n-th node on the b-th boundary
      Node* nod_pt=mesh_pt()->boundary_node_pt(b,n);

      Vector<double> x(3);
      x[0]=nod_pt->x(0);
      x[1]=nod_pt->x(1);
      x[2]=nod_pt->x(2);

      // Storage for the solution
      double u=0.0;

      // What is the exact solution at this Node?
      TanhSolnForPoisson::get_exact_u(x,u);

      // Pin this BC node
      nod_pt->pin(0);

      // Assign the solution
      nod_pt->set_value(0,u);
    }
  } // for (unsigned b=0; b<n_boundary; b++)
} // End of apply_boundary_conditions


//=====start_of_complete_element_setup====================================
/// Complete problem setup: pass pointers to physical variables.
//========================================================================
template<class TWO_D_ELEMENT,class THREE_D_ELEMENT>
void ExtrudedMovingCylinderProblem<TWO_D_ELEMENT,THREE_D_ELEMENT>::
complete_element_setup()
{
  // Complete the build of all elements so they are fully functional

  // The number of elements in the mesh
  unsigned n_element=Bulk_mesh_pt->nelement();

  // Loop over the elements to set up element-specific things that cannot
  // be handled by the (argument-free!) element constructor: Pass pointer
  // to source function
  for (unsigned i=0; i<n_element; i++)
  {
    // Upcast from GeneralsedElement to the present element
    THREE_D_ELEMENT *el_pt=
      dynamic_cast<THREE_D_ELEMENT*>(Bulk_mesh_pt->element_pt(i));

    // Set the source function pointer
    el_pt->source_fct_pt()=Source_fct_pt;
  }
} // End of complete_element_setup


//=====start_of_doc_solution==============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class TWO_D_ELEMENT,class THREE_D_ELEMENT>
void ExtrudedMovingCylinderProblem<TWO_D_ELEMENT,THREE_D_ELEMENT>::
doc_solution()
{
  // Start the clock
  double timer_s=TimingHelpers::timer();

  // Make an ofstream object to output the solution
  std::ofstream some_file;

  // Storage for the filename
  char filename[200];

  // Number of plot points to use for the big space-time solution
  unsigned n_plot_point=2;

  //-----------------
  // Output solution:
  //-----------------
  // Create the filename
  sprintf(filename,"%s/soln%i.dat",
          GlobalParameters::Doc_info.directory().c_str(),
          GlobalParameters::Doc_info.number());

  // Open a file with the constructed filename
  some_file.open(filename);

  // Output the (numerically) approximated solution
  Bulk_mesh_pt->output(some_file,n_plot_point);

  // We're done; close the file
  some_file.close();

  //-----------------------
  // Output exact solution:
  //-----------------------
  sprintf(filename,"%s/exact_soln%i.dat",
          GlobalParameters::Doc_info.directory().c_str(),
          GlobalParameters::Doc_info.number());

  // Open a file with the constructed filename
  some_file.open(filename);

  // Output the exact solution
  mesh_pt()->output_fct(some_file,n_plot_point,TanhSolnForPoisson::get_exact_u);

  // We're done; close the file
  some_file.close();

  //----------------------------------
  // Output the error in the solution:
  //----------------------------------
  // Storage for the error and solution norm
  double error=0.0, norm=0.0;

  // Create the filename
  sprintf(filename,"%s/error%i.dat",
          GlobalParameters::Doc_info.directory().c_str(),
          GlobalParameters::Doc_info.number());

  // Open a file with the constructed filename
  some_file.open(filename);

  // Output the (numerically) approximated solution
  Bulk_mesh_pt->compute_error(some_file,
                              TanhSolnForPoisson::get_exact_u,
                              error,norm);

  // We're done; close the file
  some_file.close();

  //--------------------
  // Document the error:
  //--------------------
  // Document the error
  oomph_info << "Solution norm: " << sqrt(norm)
             << "\nError norm: " << sqrt(error)
             << "\nRelative error norm: " << sqrt(error)/sqrt(norm)
             << std::endl;

  // Finally, output the time taken
  oomph_info << "\nTotal time for documentation [sec]: "
             << TimingHelpers::timer()-timer_s << std::endl;

  // Increment counter for solutions
  GlobalParameters::Doc_info.number()++;
} // End of doc_solution


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

//=====start_of_main======================================================
/// Constructor
//========================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc,argv);

  // Switch off output modifier
  oomph_info.output_modifier_pt()=&default_output_modifier;

  // Switch off oomph_info output for all processors but rank 0
  if (MPI_Helpers::communicator_pt()->my_rank()!=0)
  {
    oomph_info.stream_pt()=&oomph_nullstream;
    OomphLibWarning::set_stream_pt(&oomph_nullstream);
    OomphLibError::set_stream_pt(&oomph_nullstream);
  }
  else
  {
    oomph_info << "\n=====================================================\n"
               << "Number of processors: "
               << MPI_Helpers::communicator_pt()->nproc()
               << "\n=====================================================\n"
               << std::endl;
  }
#endif
  // Start the clock
  double timer_s=TimingHelpers::timer();

  // Set the documentation directory
  GlobalParameters::Doc_info.set_directory("RESLT");

  // Typedef the element type
  typedef RefineableQPoissonElement<2,2> ELEMENT;

  // Typedef the extruded element type
  typedef RefineableQPoissonElement<3,2> EXTRUDED_ELEMENT;

  // Create Problem
  ExtrudedMovingCylinderProblem<ELEMENT,EXTRUDED_ELEMENT> problem(
    TanhSolnForPoisson::get_source);

  // Document the mesh
  problem.newton_solve();

  // Document the mesh
  problem.doc_solution();

  // Maximum number of adaptations
  unsigned max_adapt=2;

  /// Number of uniform refinements before any solve
  for (unsigned i=0; i<max_adapt; i++)
  {
    // Document the mesh
    problem.adapt();

    // Document the mesh
    problem.newton_solve();

    // Document the mesh
    problem.doc_solution();
  }

  // Tell the user we're done
  oomph_info << "\nSimulation complete!\nTotal time for simulation [sec]: "
             << TimingHelpers::timer()-timer_s << "\n" << std::endl;
}
