//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
// Driver for 2D moving block

// Generic includes
#include "generic.h"
#include "navier_stokes.h"

// Mesh headers
#include "meshes/quad_from_triangle_mesh.h"
#include "meshes/extruded_cube_mesh_from_quad_mesh_with_macro_elements.h"

using namespace oomph;

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

//======start_of_GlobalParameters_namespace===============================
/// Global parameters for the problem
//========================================================================
namespace GlobalParameters
{
  /// Reynolds number
  double Re=100.0;

  /// Length of the mesh in the z-direction
  double L_z=3.0;

  /// Number of elements in the z-direction
  unsigned N_z=5;

  /// Helper for documenting
  DocInfo Doc_info;
} // End of GlobalParameters

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

//======start_of_problem_class============================================
/// Moving block problem
//========================================================================
template<class ELEMENT>
class MovingBlockProblem : public Problem
{
public:
  /// Constructor
  MovingBlockProblem();

  /// Destructor (empty)
  ~MovingBlockProblem() {}

  /// Update the problem specs before solve. Reset velocity boundary
  /// conditions just to be on the safe side...
  void actions_before_newton_solve()
  {
    // Call the helper function to set the boundary conditions
    apply_boundary_conditions();
  } // End of actions_before_newton_solve

  /// Update the after solve (empty)
  void actions_after_newton_solve() {}

  /// After adaptation: Unpin pressure and pin redudant pressure dofs.
  void actions_after_adapt();

  /// Doc the solution
  void doc_solution();

private:
  // Make an alias for the 3D cube element
  typedef RefineableQTaylorHoodElement<3> EXTRUDED_ELEMENT;

  /// Pointer to the 3D "extruded" mesh
  RefineableExtrudedCubeMeshFromQuadMesh<EXTRUDED_ELEMENT>* Extruded_mesh_pt;

  /// Generate a 2D mesh and "extrude" it to create a 3D mesh
  void create_extruded_mesh();

  /// Set the chosen boundary conditions
  void apply_boundary_conditions();

  /// Finish the build of the elements: assign their problem parameter pointers
  void complete_element_build();

  /// Fix pressure in element e at pressure dof pdof and set to pvalue
  void fix_pressure(const unsigned &e,
                    const unsigned &p_dof,
                    const double &p_value)
  {
    // Cast to full element type and fix the pressure at that element
    dynamic_cast<EXTRUDED_ELEMENT*>(mesh_pt()->element_pt(e))->
    fix_pressure(p_dof,p_value);
  } // End of fix_pressure

  // Enumeration of the boundaries of the 3D mesh
  enum
  {
    Lower_boundary_id=0,
    Right_boundary_id=1,
    Upper_boundary_id=2,
    Left_boundary_id=3,
    Cylinder_surface_boundary_id=4,
    Initial_time_boundary_id=5,
    Final_time_boundary_id=6
  };
}; // End of MovingBlockProblem class

//==start_of_constructor==================================================
/// Constructor for MovingBlock problem
//========================================================================
template<class ELEMENT>
MovingBlockProblem<ELEMENT>::MovingBlockProblem()
{
  // Set the maximum residuals value
  Problem::Max_residuals=100.0;

  // Create the extruded 3D mesh
  create_extruded_mesh();

  // Set the chosen boundary conditions
  apply_boundary_conditions();

  // Finish the build of the elements: assign their problem parameter pointers
  complete_element_build();

  // We're using QTaylorHood elements so we need to pin the redundant nodal
  // pressure dofs and because we're not using traction BCs, the pressure is
  // only determined up to a constant so need to add a constraint to make the
  // solution unique. For this, we set the first pressure value in element 0
  // to 0.0. This is all implemented in actions_after_adapt().
  actions_after_adapt();

  // Setup equation numbering scheme
  std::cout << "Number of equations: " << assign_eqn_numbers() << std::endl;
} // End of MovingBlockProblem


//======start_of_create_extruded_mesh======================================
/// Generate the extruded 3D mesh
//=========================================================================
template<class ELEMENT>
void MovingBlockProblem<ELEMENT>::create_extruded_mesh()
{
  // Convert arguments to strings that specify the input file names
  const string node_file_name("triangle_meshes/box_hole.1.node");
  const string element_file_name("triangle_meshes/box_hole.1.ele");
  const string poly_file_name("triangle_meshes/box_hole.1.poly");

  // Record the start time
  double start_t=TimingHelpers::timer();

  // Create the bulk mesh
  RefineableQuadFromTriangleMesh<ELEMENT>* twod_mesh_pt=
    new RefineableQuadFromTriangleMesh<ELEMENT>(
    node_file_name,element_file_name,poly_file_name);

  // Record the end time and compute the mesh setup time
  oomph_info << "Time to generate quad mesh [sec]: "
             << TimingHelpers::timer()-start_t << std::endl;

  // Output the number of elements
  oomph_info << "\nNumber of quad mesh elements: " << twod_mesh_pt->nelement()
             << "\nNumber of quad mesh nodes: " << twod_mesh_pt->nnode()
             << std::endl;

  // Record the start time
  start_t=TimingHelpers::timer();

  // Create the extruded mesh
  Extruded_mesh_pt=new RefineableExtrudedCubeMeshFromQuadMesh<EXTRUDED_ELEMENT>(
    twod_mesh_pt,GlobalParameters::N_z,GlobalParameters::L_z);

  // Record the end time and compute the mesh setup time
  oomph_info << "Time to generate extruded mesh [sec]: "
             << TimingHelpers::timer()-start_t << std::endl;

  // Create/set error estimator
  Extruded_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

  // Choose error tolerances to force some uniform refinement
  Extruded_mesh_pt->min_permitted_error()=1.0e-03;
  Extruded_mesh_pt->max_permitted_error()=3.0e-02;

  // Create the main mesh
  mesh_pt()=Extruded_mesh_pt;

  // Output the number of elements
  oomph_info << "\nNumber of extruded mesh elements: " << mesh_pt()->nelement()
             << "\nNumber of extruded mesh nodes: " << mesh_pt()->nnode()
             << std::endl;
} // End of create_extruded_mesh


//=====start_of_apply_boundary_conditions===================================
/// Set the chosen boundary conditions
//==========================================================================
template<class ELEMENT>
void MovingBlockProblem<ELEMENT>::apply_boundary_conditions()
{
  // Get the start time
  double start_t=TimingHelpers::timer();

  // Set the boundary conditions for this problem: All nodes are free by
  // default -- just pin the ones that have Dirichlet conditions here.
  // Ignore the lower boundary (b=0) so we can do it separately afterwards
  unsigned n_boundary=mesh_pt()->nboundary();
  for (unsigned b=1; b<n_boundary; b++)
  {
    unsigned n_node=mesh_pt()->nboundary_node(b);
    for (unsigned n=0; n<n_node; n++)
    {
      Node* boundary_node_pt=mesh_pt()->boundary_node_pt(b,n);

      unsigned n_value=boundary_node_pt->nvalue();

      // Just pin the velocity dofs
      for (unsigned i=0; i<n_value-1; i++)
      {
        // Pin all boundary nodes; we're applying Dirichlet BCs
        boundary_node_pt->pin(i);

        // No-slip BCs on the side walls
        boundary_node_pt->set_value(i,0.0);
      }
    } // for (unsigned n=0; n<n_node; n++)
  } // for (unsigned b=0; b<n_boundary; b++)

  // Set the boundary conditions for this problem: All nodes are free by
  // default -- just pin the ones that have Dirichlet conditions here.
  unsigned n_node=mesh_pt()->nboundary_node(Lower_boundary_id);
  for (unsigned n=0; n<n_node; n++)
  {
    // Pointer to the n-th node on the lower boundary
    Node* boundary_node_pt=mesh_pt()->boundary_node_pt(Lower_boundary_id,n);

    // The number of N.St. unknowns
    unsigned n_value=boundary_node_pt->nvalue();

    // Just pin the velocity dofs
    for (unsigned i=0; i<n_value-1; i++)
    {
      // Pin all boundary nodes; we're applying Dirichlet BCs
      boundary_node_pt->pin(i);
    }

    // Horizontal flow
    boundary_node_pt->set_value(0,1.0);
    boundary_node_pt->set_value(1,0.0);
    boundary_node_pt->set_value(2,0.0);
  } // for (unsigned n=0; n<n_node; n++)

  // Output the setup time to the screen
  oomph_info << "\nTime taken for application of boundary conditions [sec]: "
             << TimingHelpers::timer()-start_t << std::endl;
} // End of apply_boundary_conditions


//=====start_of_complete_element_build======================================
/// Complete problem setup: pass pointers to physical variables.
//==========================================================================
template<class ELEMENT>
void MovingBlockProblem<ELEMENT>::complete_element_build()
{
  // Get the start time
  double start_t=TimingHelpers::timer();

  // Find number of elements in mesh
  unsigned n_element=mesh_pt()->nelement();

  // Loop over the elements to set up element-specific things that cannot be
  // handled by constructor
  for (unsigned e=0; e<n_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    EXTRUDED_ELEMENT* el_pt=
      dynamic_cast<EXTRUDED_ELEMENT*>(mesh_pt()->element_pt(e));

    // Set the Reynolds number
    el_pt->re_pt()=&GlobalParameters::Re;
  }

  // Output the setup time to the screen
  oomph_info << "\nTime taken to complete build of elements [sec]: "
             << TimingHelpers::timer()-start_t << std::endl;
} // End of complete_element_build


//======start_of_actions_after_adapt========================================
/// After adaptation: Unpin pressure and pin redudant pressure dofs.
//==========================================================================
template<class ELEMENT>
void MovingBlockProblem<ELEMENT>::actions_after_adapt()
{
  // Unpin all pressure dofs
  RefineableNavierStokesEquations<3>::
  unpin_all_pressure_dofs(mesh_pt()->element_pt());

  // Pin redundant pressure dofs
  RefineableNavierStokesEquations<3>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());

  // Now set the first pressure dof in the first element to 0.0
  fix_pressure(0,0,0.0);
} // End of actions_after_adapt


//======start_of_doc_solution=============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void MovingBlockProblem<ELEMENT>::doc_solution()
{
  // Stream to output the data to the output file
  std::ofstream some_file;

  // Storage for the filename
  char filename[200];

  // Number of plot points
  unsigned n_plot_point=2;

  // Output solution
  sprintf(filename,"%s/soln%i.dat",
          GlobalParameters::Doc_info.directory().c_str(),
          GlobalParameters::Doc_info.number());

  // Document the solution
  some_file.open(filename);
  mesh_pt()->output(some_file,n_plot_point);
  some_file.close();

  // Increment the counter
  GlobalParameters::Doc_info.number()++;
} // End of doc_solution

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//==start_of_main======================================================
/// Driver for MovingBlock test problem -- test drive
/// with two different types of element.
//=====================================================================
int main(int argc, char *argv[])
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

  // Record the start time
  double start_t=TimingHelpers::timer();

  // Set output directory
  GlobalParameters::Doc_info.set_directory("RESLT");

  // Step number
  GlobalParameters::Doc_info.number()=0;

  // Typedef the ELEMENT and MESH type
  typedef RefineableQTaylorHoodElement<2> ELEMENT;

  // Build the problem with QCrouzeixRaviartElements
  MovingBlockProblem<ELEMENT> problem;

  // Tell the user
  oomph_info << "Using RefineableQTaylorHoodElement<2>" << std::endl;

  // Solve the problem
  problem.newton_solve();

  // Output the solution
  problem.doc_solution();

  // Solve the problem
  problem.adapt();

  // Solve the problem
  problem.newton_solve();

  // Output the solution
  problem.doc_solution();

  // Record the end time and compute the mesh setup time
  oomph_info << "Simulation complete.\n\nTime taken for simulation [sec]: "
             << TimingHelpers::timer()-start_t << std::endl;

#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::finalize();
#endif
} // end_of_main

