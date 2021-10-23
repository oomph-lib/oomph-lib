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
// Generic routines
#include "generic.h"

// The space-time Unsteady heat equations
#include "space_time_unsteady_heat_equal_order_galerkin_petrov.h"

// The mesh
#include "meshes/simple_cubic_mesh.h"

// Add in the block preconditioning machinery
#include "space_time_block_preconditioner.h"

// Include oomph namespace
using namespace oomph;

using namespace std;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//======start_of_GlobalParameters==========================================
/// Global parameters for the problem
//=========================================================================
namespace GlobalParameters
{
  ///---------------------Unsteady Heat Parameters------------------------
  /// Alpha value (thermal inertia)
  double Alpha=1.0;

  /// Beta value (thermal conductivity)
  double Beta=1.0;
  ///---------------------Unsteady Heat Parameters------------------------

  ///--------------------------Mesh Properties----------------------------
  /// Length of the mesh in the x-direction
  double L_x=1.0;

  /// Length of the mesh in the y-direction
  double L_y=1.0;

  /// Length of the mesh in the t-direction
  double L_t=1.0;

  /// Number of elements in the x-direction
  unsigned N_x=5;

  /// Number of elements in the y-direction
  unsigned N_y=5;

  /// Number of elements in the t-direction
  unsigned N_t=5;

  /// Should we apply time-periodic boundary conditions?
  bool Apply_time_periodic_boundary_conditions=true;
  ///--------------------------Mesh Properties----------------------------

  ///----------------------------Solver Info------------------------------
  /// Variable to choose which preconditioner to use. The actual
  /// preconditioner we choose to use is defined by the enumeration class
  /// implemented in the problem
  unsigned Preconditioner=1;

  /// Storage for the number of dof types in the mesh. Will be
  /// assigned in the function assign_time_slice_id()
  unsigned N_dof_type=0;

  /// Helper function which sets up the mapping between DOF types
  /// and which block they should be assigned to. This relies on the concept
  /// of "time slabs" in the space-time formulation. All dofs in a given
  /// time slab will be aggregrated together
  void set_up_dof_to_block_mapping(Vector<unsigned>& dof_to_block_map)
  {
    // Resize the vector
    dof_to_block_map.resize(N_dof_type);

    // Loop over the dofs
    for (unsigned i=0; i<N_dof_type; i++)
    {
      // Each time slice should correspond to its own block
      dof_to_block_map[i]=i;
    }
  } // End of set_up_dof_to_block_mapping
  ///----------------------------Solver Info------------------------------

  ///-----------------------Documentation Helpers-------------------------
  // DocInfo object for documentation
  DocInfo Doc_info;
  ///-----------------------Documentation Helpers-------------------------

  ///---------------------------Miscellaneous-----------------------------
  /// Function to round a double to the nearest integral value
  double round(const double& d)
  {
    // Round it
    return std::floor(d+0.5);
  } // End of round
  ///---------------------------Miscellaneous-----------------------------
} // End of GlobalParameters

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//======start_of_SinSolution===============================================
/// Namespace for forced exact solution for UnsteadyHeat equation
//=========================================================================
namespace SinSolution
{
  // Shift the solution from the origin
  double U_shift=8.0;

  /// Exact solution as a vector
  void get_exact_u(const double& t,
                   const Vector<double>& x,
                   Vector<double>& u)
  {
    // Cache the value of pi
    double pi=MathematicalConstants::Pi;

    // The full solution
    u[0]=sin(pi*x[0])*sin(pi*x[1])*cos(2.0*pi*t)+U_shift;
  } // End of get_exact_u

  /// Exact solution as a scalar
  void get_exact_u(const double& t,
                   const Vector<double>& x,
                   double& u)
  {
    // Cache the value of pi
    double pi=MathematicalConstants::Pi;

    // The full solution
    u=sin(pi*x[0])*sin(pi*x[1])*cos(2.0*pi*t)+U_shift;
  } // End of get_exact_u


  /// Source function to make it an exact solution
  void get_source(const double& t,
                  const Vector<double>& x,
                  double& source)
  {
    // Cache the value of pi
    double pi=MathematicalConstants::Pi;

    // The temporal derivative of the solution
    double du_dt=-2.0*pi*sin(pi*x[0])*sin(pi*x[1])*sin(2.0*pi*t);

    // Second derivative of spatial contribution w.r.t. x
    double d2u_dx2=-pow(pi,2)*sin(pi*x[0])*sin(pi*x[1])*cos(2.0*pi*t);

    // Second derivative of spatial contribution w.r.t. y
    double d2u_dy2=-pow(pi,2)*sin(pi*x[0])*sin(pi*x[1])*cos(2.0*pi*t);

    // Assign the source function value
    source=-GlobalParameters::Alpha*du_dt+
           GlobalParameters::Beta*(d2u_dx2+d2u_dy2);
  } // End of get_source
} // End of SinSolution

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//======start_of_UnsteadyHeatProblem========================================
/// UnsteadyHeat problem
//=========================================================================
template<class ELEMENT>
class UnsteadyHeatProblem : public virtual Problem
{
  // Typedef the source function pointer type to keep things tidy
  typedef SpaceTimeUnsteadyHeatEquations<2>::
  SpaceTimeUnsteadyHeatSourceFctPt SpaceTimeFunctionPt;

public:
  /// Constructor
  UnsteadyHeatProblem(const SpaceTimeFunctionPt& source_fct_pt);

  /// Destructor (empty)
  ~UnsteadyHeatProblem();

  /// Update the problem specs before solve (empty)
  void actions_before_newton_solve() {}

  /// Update the problem specs after solve (empty)
  void actions_after_newton_solve() {}

  /// Update the problem specs before solve (empty)
  void actions_before_adapt() {}

  /// Update the problem specs after solve
  void actions_after_adapt()
  {
    // Assign the appropriate boundary conditions
    apply_boundary_conditions();

    // Re-assign time-slab IDs and wipe the old solver/preconditioner and
    // create a new one based on the (possibly) updated dof-to-block map
    update_block_preconditioner_after_refinement();

    // Complete the setup of the elements
    complete_problem_setup();
  } // End of actions_after_adapt

  /// A wrapper around the Problem's refine_uniformly(...) function;
  /// update the number of space-time slabs after a uniform refinement so
  /// that the preconditioner solve the system more efficiently.
  /// NOTE: We don't implement this for
  void refine_uniformly()
  {
    // The number of space-time slabs will double after this refinement
    N_space_time_slab*=2;

    // Call the Problem's implementation of refine_uniformly(...)
    Problem::refine_uniformly();
  } // End of refine_uniformly

  /// A wrapper around the Problem's adapt(...) function; issue a
  /// warning that adaptation doesn't work for these elements
  void adapt()
  {
    // Create an output stream
    std::ostringstream warning_message_stream;

    // Create an error message
    warning_message_stream << "Hey, adaptation doesn't work properly with "
                           << "Petrov-Galerkin elements!\nI'm not going to "
                           << "stop the simulation because you may already "
                           << "know this\nbut if you didn't, consider "
                           << "yourself warned!" << std::endl;

    // Throw an error
    OomphLibWarning(warning_message_stream.str(),
                    OOMPH_CURRENT_FUNCTION,
                    OOMPH_EXCEPTION_LOCATION);

    // Call the Problem's implementation of refine_uniformly(...)
    Problem::adapt();
  } // End of adapt

  /// Create the space-time mesh with the chosen number of elements in
  /// the time direction and the chosen spatial resolution
  void create_spacetime_mesh();

  /// Assign the Dirichlet and maybe time-periodic BCs
  void apply_boundary_conditions();

  /// Assign the appropriate boundary conditions and enforce periodicity
  /// in the time direction
  void enforce_time_periodic_boundary_conditions();

  /// Function to set the periodicity between two octrees and assign
  /// the up and right equivalents to both elements
  void set_neighbour_periodic_and_up_right_equivalents(FiniteElement* el0_pt,
      FiniteElement* el1_pt,
      const int& direction0);

  // /// Assign the dof type to each dof according to the temporal position
  // /// of their associated Node. This allows the block preconditioner to group
  // /// dofs in the same time-slice together.
  // /// NOTE: The boundary conditions need to be applied before this function is
  // /// invoked to ensure that the correct number of time slide IDs are assigned.
  // void assign_time_slice_id();

  /// Helper function when space-time block preconditioning is being used
  void assign_time_slab_id();

  /// Helper function to update the block preconditioner after, what
  /// seems like, a uniform refinement
  void update_block_preconditioner_after_refinement();

  /// Complete problem setup; make all the elements fully functional
  /// by passing pointers to all physical parameters
  void complete_problem_setup();

  /// Assign the chosen solver (and preconditioner if so desired)
  void set_up_spacetime_solver();

  /// Doc the solution
  void doc_solution(const bool& doc_error=false);

private:
  /// Oomph-lib iterative linear solver
  IterativeLinearSolver* Solver_pt;

  /// Preconditioner
  Preconditioner* Prec_pt;

  // Enumeration of the preconditioners
  enum
  {
    Diagonal_preconditioner=0,
    Lower_triangular_preconditioner=1
  };

  /// Pointer to the mesh
  RefineableSimpleCubicMesh<ELEMENT>* Bulk_mesh_pt;

  /// The number of space-time slabs in the mesh; this is essentially
  /// a copy of N_t in the GlobalParameters namespace but will be updated
  /// when we complete a uniform refinement
  unsigned N_space_time_slab;

  /// Pointer to source function
  SpaceTimeFunctionPt Source_fct_pt;

  // Boolean variable to check if periodicity has been set up
  bool Periodicity_has_been_enforced;

  // Enumeration of the boundaries of the space-time mesh
  enum
  {
    Initial_time_boundary_id=0,
    Lower_spatial_boundary_id=1,
    Right_spatial_boundary_id=2,
    Upper_spatial_boundary_id=3,
    Left_spatial_boundary_id=4,
    Final_time_boundary_id=5
  };
}; // End of UnsteadyHeatProblem class


//======start_of_constructor=============================================
/// Constructor for UnsteadyHeat problem in cubic domain
//=========================================================================
template<class ELEMENT>
UnsteadyHeatProblem<ELEMENT>::UnsteadyHeatProblem(
  const SpaceTimeFunctionPt& source_fct_pt) :
  Solver_pt(0),
  Prec_pt(0),
  Bulk_mesh_pt(0),
  Source_fct_pt(source_fct_pt),
  Periodicity_has_been_enforced(false)
{
  // Generate the space-time mesh (stored as Bulk_mesh_pt)
  create_spacetime_mesh();

  // Assign the appropriate boundary conditions
  apply_boundary_conditions();

  // Assign the time-slab IDs (for the block preconditioner)
  assign_time_slab_id();
  // assign_time_slice_id();

  // Complete the setup of the elements
  complete_problem_setup();

  // Call the auxiliary solver setup function
  set_up_spacetime_solver();

  // Do equation numbering
  oomph_info << "\nNumber of equations: " << assign_eqn_numbers() << std::endl;
} // End of UnsteadyHeatProblem


//======start_of_destructor==============================================
/// Destructor for UnsteadyHeat problem in cubic domain
//=========================================================================
template<class ELEMENT>
UnsteadyHeatProblem<ELEMENT>::~UnsteadyHeatProblem()
{
  // If we're not using the default linear solver
  if (Prec_pt!=0)
  {
    // Kill the preconditioner
    delete Prec_pt;

    // Make it a null pointer
    Prec_pt=0;
  }

  // If we're not using the default linear solver
  if (Solver_pt!=0)
  {
    // Kill the linear solver
    delete Solver_pt;

    // Make it a null pointer
    Solver_pt=0;

    // Make the Problem's linear solver pointer null too
    linear_solver_pt()=0;
  }

  // Delete the error estimator
  delete Bulk_mesh_pt->spatial_error_estimator_pt();

  // Make it a null pointer
  Bulk_mesh_pt->spatial_error_estimator_pt()=0;

  // Delete the mesh
  delete Bulk_mesh_pt;

  // Set the pointer to null
  Bulk_mesh_pt=0;
} // End of ~UnsteadyHeatProblem


//======start_of_create_spacetime_mesh======================================
/// Helper function to create the space-time mesh (to be assigned to
/// Bulk_mesh_pt) with the chosen number of elements in the time direction
/// and an appropriate spatial resolution (to capture the time-periodic
/// solution properly).
//=========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::create_spacetime_mesh()
{
  // The number of space-time slabs in the mesh
  N_space_time_slab=GlobalParameters::N_t;

  // Build the refineable mesh
  Bulk_mesh_pt=new RefineableSimpleCubicMesh<ELEMENT>
  (GlobalParameters::N_x,GlobalParameters::N_y,GlobalParameters::N_t,
   GlobalParameters::L_x,GlobalParameters::L_y,GlobalParameters::L_t);

  // Assign the mesh pointer
  mesh_pt()=Bulk_mesh_pt;

  // Add an error estimator
  Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

  // Maximum permitted errors
  Bulk_mesh_pt->max_permitted_error()=1.0e-02;

  // Minimum permitted errors
  Bulk_mesh_pt->min_permitted_error()=1.0e-03;
} // End of create_spacetime_mesh


//=====start_of_assign_time_slab_id=========================================
/// Assign the time slice IDs to each BlockPreconditionable element
//==========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::assign_time_slab_id()
{
  // Number of dimensions
  const unsigned n_dim=3;

  // The ID of the coordinate we want, i.e. the ID of the time-direction
  const unsigned t_index=n_dim-1;

  // Space for the local coordinates (at the center of the element)
  Vector<double> s(n_dim,0.0);

  // Space for the (Eulerian) coordinates
  double t=0.0;

  // Get the number of elements in the mesh
  const unsigned n_element=Bulk_mesh_pt->nelement();

  // Loop over the elements
  for (unsigned i=0; i<n_element; i++)
  {
    // Upcast the element
    ELEMENT* const el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

    // Get the interpolated time value at the center of the element
    t=el_pt->interpolated_x(s,t_index);

    // The idea is to group all dofs in each space-time slab together, by
    // making sure that elements in the same space-time slab have the same ID

    // To ensure this is done correctly, we provide an offset. Note, however
    unsigned id=std::floor((double(N_space_time_slab)*t)/
                           GlobalParameters::L_t);

    // Upcast and assign the time slice ID
    el_pt->set_time_slab_id(id);

    // We're grouping dofs in each space-time slab together and there are only
    // "N_space_time_slab" space-time slabs
    GlobalParameters::N_dof_type=N_space_time_slab;

    // Finally, tell it how many dof types there are in the mesh
    el_pt->set_ndof_types(GlobalParameters::N_dof_type);
  } // for (unsigned i=0;i<n_element;i++)
} // End of assign_time_slab_id


//========start_of_update_block_preconditioner_after_refinement===========
/// Helper function to update the block preconditioner after, what
/// seems like, a uniform refinement
//========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::
update_block_preconditioner_after_refinement()
{
  // Make sure we're even using the preconditioner
  if (Prec_pt!=0)
  {
    // Set the time-slab id
    assign_time_slab_id();

    // Don't waste time, just kill the old solver
    // and set up a new one. God knows why it's
    // struggling with me trying to replace the
    // current dof-to-block map...

    // Kill the preconditioner
    delete Prec_pt;

    // Kill the solver
    delete Solver_pt;

    // Make the reference to the preconditioner a null pointer
    Prec_pt=0;

    // Make the reference to the solver a null pointer
    Solver_pt=0;

    // Set up a new solver and preconditioner
    set_up_spacetime_solver();
  }
} // End of update_block_preconditioner_after_refinement


//======start_of_set_up_spacetime_solver===================================
/// Set up the solver for this Problem object
//=========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::set_up_spacetime_solver()
{
  // Create oomph-lib iterative linear solver
  Solver_pt=new GMRES<CRDoubleMatrix>;

  // Use LHS preconditioning
  dynamic_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_RHS();

  // Set the tolerance
  Solver_pt->tolerance()=1.0e-12;

  // Maximum number of iterations
  Solver_pt->max_iter()=200;

  // Set linear solver
  linear_solver_pt()=Solver_pt;

  //-----------------------------------------
  // Create the master-level preconditioners:
  //-----------------------------------------
  // Solve the diagonal blocks associated with each time-slice separately
  if (GlobalParameters::Preconditioner==Diagonal_preconditioner)
  {
    // Create a new instance of the space-time preconditioner
    Prec_pt=new BlockDiagonalPreconditioner<CRDoubleMatrix>;
  }
  // Solve the block lower-triangular part of the system matrix
  else if (GlobalParameters::Preconditioner==Lower_triangular_preconditioner)
  {
    // Create a new instance of the space-time preconditioner
    Prec_pt=new BandedBlockTriangularPreconditioner<CRDoubleMatrix>;

    // Create a new instance of the space-time preconditioner
    BandedBlockTriangularPreconditioner<CRDoubleMatrix>* st_block_prec_pt=
      dynamic_cast<BandedBlockTriangularPreconditioner<CRDoubleMatrix>*>(
        Prec_pt);

    // Indicate that we're using a block lower triangular solve
    st_block_prec_pt->lower_triangular();

    // The block bandwidth for the preconditioner; because we group dofs by
    // their space-time-slab, the bandwidth will be one (we only have diagonal
    // blocks and subdiagonal blocks in the block-organised system matrix for
    // the backwards-looking mixed order discretisation).
    // NOTE: Don't delete this as this drastically reduces the time taken to
    // to extract blocks during the preconditioner setup.
    bool bandwidth=1;

    // Provide the bandwidth; only subdiagonal block entries
    st_block_prec_pt->set_block_bandwidth(bandwidth);

    // Do we want to document the memory usage?
    bool document_memory_usage=true;

    // If we want to document the memory usage
    if (document_memory_usage)
    {
      // Tell the block preconditioner that we want it to doc. the memory usage
      st_block_prec_pt->enable_doc_memory_usage();
    }
  }
  // If the user provided an invalid input
  else
  {
    // Throw an error
    throw OomphLibError("Invalid choice of preconditioner.",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  } // if (GlobalParameters::Preconditioner==Diagonal_preconditioner)

  // Allocate space for the DOF to block map; this tells the master block
  // preconditioner which DOF types to aggregate
  Vector<unsigned> dof_to_block_map;

  // Call the auxiliary function which sets up the mapping
  GlobalParameters::set_up_dof_to_block_mapping(dof_to_block_map);

  // Create an upcasted pointer to the master preconditioner
  GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* upcasted_master_prec_pt=
    dynamic_cast<GeneralPurposeBlockPreconditioner<CRDoubleMatrix>*>(Prec_pt);

  // Build silently!
  upcasted_master_prec_pt->enable_silent_preconditioner_setup();

  // Pass the DOF-to-block map to the preconditioner
  upcasted_master_prec_pt->set_dof_to_block_map(dof_to_block_map);

  // Pass a pointer to the (space-time) mesh
  upcasted_master_prec_pt->add_mesh(Bulk_mesh_pt);

  // Now assign the preconditioner to the linear solver
  Solver_pt->preconditioner_pt()=Prec_pt;
} // End of set_up_spacetime_solver


//======start_of_complete_problem_setup====================================
/// Complete problem setup: pass pointers to physical variables.
//=========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::complete_problem_setup()
{
  // Find number of elements in mesh
  unsigned n_element=Bulk_mesh_pt->nelement();

  // Loop over the elements to set up element-specific things that
  // cannot be handled by the constructor
  for (unsigned i=0; i<n_element; i++)
  {
    // Upcast from FiniteElement to the present element
    ELEMENT* el_pt=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

    // Set the source function pointer
    el_pt->source_fct_pt()=Source_fct_pt;

    // Set the alpha value pointer
    el_pt->alpha_pt()=&GlobalParameters::Alpha;

    // Set the beta value pointer
    el_pt->beta_pt()=&GlobalParameters::Beta;
  }
} // End of complete_problem_setup


//======start_of_apply_boundary_conditions=================================
/// Apply the Dirichlet conditions on the spatial boundaries and
/// the initial time boundary, if specified, otherwise apply time-periodic
/// boundary conditions on the time boundaries. More explicitly:
///             Boundary 0 (t=0) -- Dirichlet or time-periodic BCs
///             Boundary 1 (y=0) -- Dirichlet
///             Boundary 2 (x=1) -- Dirichlet
///             Boundary 3 (y=1) -- Dirichlet
///             Boundary 4 (x=0) -- Dirichlet
///             Boundary 5 (t=1) -- Do nothing!
//=========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::apply_boundary_conditions()
{
  // Number of spatial dimensions
  unsigned n_dim=2;

  // Storage for the time
  double time=0.0;

  // Storage for the spatial coordinates
  Vector<double> spatial_coordinates(n_dim,0.0);

  // Storage for the solution
  double u_exact=0.0;

  // Get the number of boundaries in the mesh
  unsigned n_bound=Bulk_mesh_pt->nboundary();

  // Loop over the boundaries
  for (unsigned b=0; b<n_bound; b++)
  {
    // Get the number of nodes on the b-th boundary
    unsigned n_node=Bulk_mesh_pt->nboundary_node(b);

    // Don't do anything on the final time boundary
    if (b!=Final_time_boundary_id)
    {
      // If we're applying time-periodic boundary conditions and we're
      if ((GlobalParameters::Apply_time_periodic_boundary_conditions)&&
          (b==Initial_time_boundary_id))
      {
        // We need to apply time-periodic BCs
        enforce_time_periodic_boundary_conditions();
      }
      // Otherwise apply Dirichlet boundary conditions
      else
      {
        // Loop over the nodes on the b-th boundary
        for (unsigned n=0; n<n_node; n++)
        {
          // Get a pointer to the n-th node on the b-th boundary
          Node* node_pt=Bulk_mesh_pt->boundary_node_pt(b,n);

          // Pin the (one and only) dof at this node!
          node_pt->pin(0);

          // Loop over the coordinates
          for (unsigned i=0; i<n_dim; i++)
          {
            // Get the i-th coordinate of the node
            spatial_coordinates[i]=node_pt->x(i);
          }

          // Get the current time
          time=node_pt->x(n_dim);

          // Get the exact solution at this node
          SinSolution::get_exact_u(time,spatial_coordinates,u_exact);

          // Set the value of the solution
          node_pt->set_value(0,u_exact);
        }
      } // if ((GlobalParameters::Apply_time_periodic_boundary_conditions)&&
    } // if (b!=Final_time_boundary_id)
  } // for (unsigned b=0;b<n_bound;b++)
} // End of apply_boundary_conditions


//======start_of_enforce_time_periodic_boundary_conditions=================
/// Assign the appropriate boundary conditions, i.e. periodicity
/// in the t-direction. In the x and y-direction apply Dirichlet boundary
/// conditions. In summary:
///             Boundary 0 (t=0) -- Periodic in time (w.r.t. boundary 5)
///             Boundary 1 (y=0) -- Dirichlet
///             Boundary 2 (x=1) -- Dirichlet
///             Boundary 3 (y=1) -- Dirichlet
///             Boundary 4 (x=0) -- Dirichlet
///             Boundary 5 (t=1) -- Periodic in time (w.r.t. boundary 0)
//=========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::enforce_time_periodic_boundary_conditions()
{
  // If we need to set up periodicity
  if (!Periodicity_has_been_enforced)
  {
    // Number of dimensions
    unsigned n_dim=3;

    // Index of the t=0 boundary
    unsigned boundary0_index=0;

    // Index of the t=1 boundary
    unsigned boundary1_index=5;

    // Number of nodes on t=0 boundary
    unsigned n_boundary0_node=Bulk_mesh_pt->nboundary_node(boundary0_index);

    // Number of nodes on t=1 boundary
    unsigned n_boundary1_node=Bulk_mesh_pt->nboundary_node(boundary1_index);

    //----------------------------
    // Establish nodal periodicity
    //----------------------------
    // Make sure there are as many nodes on boundary 0 as there are on
    // boundary (n_boundary-1)
    if (n_boundary0_node!=n_boundary1_node)
    {
      // Create an output stream
      std::ofstream output_file;

      // Open a file
      output_file.open("RESLT/nodes_b0.csv");

      // Loop over the nodes on the t=0 boundary
      for (unsigned i=0; i<n_boundary0_node; i++)
      {
        // Output the coordinates of the i-th node on the t=0 boundary
        Bulk_mesh_pt->boundary_node_pt(boundary0_index,i)->output(output_file);
      }

      // Close the file
      output_file.close();

      // Open a file
      output_file.open("RESLT/nodes_b1.csv");

      // Loop over the nodes on the t=0 boundary
      for (unsigned i=0; i<n_boundary1_node; i++)
      {
        // Output the coordinates of the i-th node on the t=1 boundary
        Bulk_mesh_pt->boundary_node_pt(boundary1_index,i)->output(output_file);
      }

      // Close the file
      output_file.close();

      // Create an output stream
      std::ostringstream error_message_stream;

      // Create an error message
      error_message_stream << "Different number of nodes on t=0 and t=1 "
                           << "boundary!\nThere are " << n_boundary0_node
                           << " nodes on boundary " << boundary0_index
                           << " and " << n_boundary1_node
                           << " nodes on boundary " << boundary1_index
                           << "!" << std::endl;

      // Throw an error
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Loop over the nodes on the t=0 boundary
    for (unsigned i=0; i<n_boundary0_node; i++)
    {
      // Get the pointer to the associated node
      Node* node0_pt=Bulk_mesh_pt->boundary_node_pt(boundary0_index,i);

      // Boolean to indicate whether or not the neighbour has been found
      bool has_neighbour_node_been_found=false;

      // Loop over the nodes on the t=1 boundary
      for (unsigned j=0; i<n_boundary1_node; j++)
      {
        // Get the pointer to the associated node
        Node* node1_pt=Bulk_mesh_pt->boundary_node_pt(boundary1_index,j);

        // Distance value
        double distance=0.0;

        // Loop over the entries of x
        for (unsigned k=0; k<n_dim-1; k++)
        {
          // Update the distance (2 norm)
          distance+=pow(((node0_pt->x(k))-(node1_pt->x(k))),2.0);
        }

        // Square root it
        distance=std::sqrt(distance);

        // Check if it matches to within a reasonable tolerance
        if (std::fabs(distance)<Tree::max_neighbour_finding_tolerance())
        {
          // Make the nodes periodic; the node on the t=0 boundary now points
          // to the node on the t=1 boundary.
          node0_pt->make_periodic(node1_pt);

          // We've found the neighbouring node
          has_neighbour_node_been_found=true;

          // We're done; break out!
          break;
        }
      } // for (unsigned i=0;i<n_boundary0_node;i++)

      // If we get here and we haven't found the neighbouring node, something's
      // wrong so throw an error
      if (!has_neighbour_node_been_found)
      {
        // Throw an error
        throw OomphLibError("Couldn't find neighbouring node!",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    } // for (unsigned i=0;i<n_boundary0_node;i++)

    //--------------------------------
    // Establish elemental periodicity
    //--------------------------------
    // Number of elements on t=0 boundary
    unsigned n_boundary0_element=
      Bulk_mesh_pt->nboundary_element(boundary0_index);

    // Number of elements on t=1 boundary
    unsigned n_boundary1_element=
      Bulk_mesh_pt->nboundary_element(boundary1_index);

    // Make sure there are as many nodes on boundary 0 as there are on
    // boundary (n_boundary-1)
    if (n_boundary0_element!=n_boundary1_element)
    {
      // Throw an error
      throw OomphLibError("Different number of elements on time boundaries!",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Storage for the local coordinates of the centre of the face of the
    // element on the t=0 boundary
    Vector<double> s_back(n_dim,0.0);

    // Storage for the local coordinates of the centre of the face of the
    // element on the t=1 boundary
    Vector<double> s_front(n_dim,0.0);

    // The t=0 boundary corresponds to the s[n_dim-1]=-1 (in the
    // appropriate element)
    s_back[n_dim-1]=-1.0;

    // The t=1 boundary corresponds to the s[n_dim-1]=1 (in the
    // appropriate element)
    s_front[n_dim-1]=1.0;

    // Loop over the elements on t=0 boundary
    for (unsigned i=0; i<n_boundary0_element; i++)
    {
      // Get a pointer to the associated element
      FiniteElement* el0_pt=Bulk_mesh_pt->boundary_element_pt(boundary0_index,i);

      // Boolean to indicate whether or not the neighbour has been found
      bool has_neighbour_element_been_found=false;

      // Loop over the elements on t=1 boundary
      for (unsigned j=0; j<n_boundary1_element; j++)
      {
        // Get a pointer to the associated element
        FiniteElement* el1_pt=Bulk_mesh_pt->boundary_element_pt(boundary1_index,j);

        // Distance value
        double distance=0.0;

        // Loop over the entries of x
        for (unsigned k=0; k<n_dim-1; k++)
        {
          // Update the distance (2 norm)
          distance+=pow((el0_pt->interpolated_x(s_back,k))-
                        (el1_pt->interpolated_x(s_front,k)),2.0);
        }

        // Square root it
        distance=std::sqrt(distance);

        // Check if it matches to within a reasonable tolerance
        if (std::fabs(distance)<Tree::max_neighbour_finding_tolerance())
        {
          //-----------------------------------------------------------
          //                   ...FACE NEIGHBOURS...
          //-----------------------------------------------------------
          //-----------------------------------------------------------
          // First task: Set the periodicity in the face direction
          // and assign the neighbour pointers.
          //-----------------------------------------------------------
          // Variable to hold the face direction of the neighbouring element
          // from the perspective of the current element (t=1 boundary element)
          int direction1=OcTreeNames::F;

          // Set the periodicity between the neighbouring OcTree objects and
          // set the up/right equivalents of the elements (it is assumed that
          // the elements are both oriented in the same way -- seems reasonable
          // for an extruded mesh)
          set_neighbour_periodic_and_up_right_equivalents(el1_pt,el0_pt,
              direction1);

          //-----------------------------------------------------------
          //                   ...EDGE NEIGHBOURS...
          //-----------------------------------------------------------
          //-----------------------------------------------------------
          // First task: Find an element on the t=0 boundary and its
          // corresponding element on the t=1 boundary in a given edge
          // direction. We already have a pointer to the element on the
          // back face (t=0 boundary) pointed to by el0_pt so we don't
          // need to calculate that. To get the edge neighbours (on the
          // t=1 boundary) we can simply use the OcTree data structure.
          //-----------------------------------------------------------
          // Pointer to the tree root associated with the element on
          // the t=0 boundary
          OcTreeRoot* octree0_pt=
            dynamic_cast<OcTreeRoot*>
            (dynamic_cast<ELEMENT*>(el0_pt)->tree_pt()->root_pt());

          // First edge direction
          direction1=OcTreeNames::LF;

          // If there's an face neighbour in the given direction
          if (octree0_pt->neighbour_pt(OcTreeNames::L)!=0)
          {
            // Get a pointer to the LF edge neighbour of the t=1 boundary element
            el0_pt=octree0_pt->neighbour_pt(OcTreeNames::L)->object_pt();

            // Set the periodicity between the neighbouring OcTree objects and
            // set the up/right equivalents of the elements (it is assumed that
            // the elements are both oriented in the same way -- seems reasonable
            // for an extruded mesh)
            set_neighbour_periodic_and_up_right_equivalents(
              el1_pt,el0_pt,direction1);
          }

          // Second edge direction
          direction1=OcTreeNames::RF;

          // If there's an face neighbour in the given direction
          if (octree0_pt->neighbour_pt(OcTreeNames::R)!=0)
          {
            // Get a pointer to the RF edge neighbour of the t=1 boundary element
            el0_pt=octree0_pt->neighbour_pt(OcTreeNames::R)->object_pt();

            // Set the periodicity between the neighbouring OcTree objects and
            // set the up/right equivalents of the elements
            set_neighbour_periodic_and_up_right_equivalents(
              el1_pt,el0_pt,direction1);
          }

          // Third edge direction
          direction1=OcTreeNames::DF;

          // If there's an face neighbour in the given direction
          if (octree0_pt->neighbour_pt(OcTreeNames::D)!=0)
          {
            // Get a pointer to the DF edge neighbour of the t=1 boundary element
            el0_pt=octree0_pt->neighbour_pt(OcTreeNames::D)->object_pt();

            // Set the periodicity between the neighbouring OcTree objects and
            // set the up/right equivalents of the elements
            set_neighbour_periodic_and_up_right_equivalents(
              el1_pt,el0_pt,direction1);
          }

          // Fourth edge direction
          direction1=OcTreeNames::UF;

          // If there's an face neighbour in the given direction
          if (octree0_pt->neighbour_pt(OcTreeNames::U)!=0)
          {
            // Get a pointer to the UF edge neighbour of the t=1 boundary element
            el0_pt=octree0_pt->neighbour_pt(OcTreeNames::U)->object_pt();

            // Set the periodicity between the neighbouring OcTree objects and
            // set the up/right equivalents of the elements
            set_neighbour_periodic_and_up_right_equivalents(
              el1_pt,el0_pt,direction1);
          }

          // We've found the neighbouring node
          has_neighbour_element_been_found=true;

          // We're done; break out!
          break;
        }
      } // for (unsigned j=0;j<n_boundary1_element;j++)

      // If we get here and we haven't found the neighbouring element,
      // something's wrong
      if (!has_neighbour_element_been_found)
      {
        // Throw an error
        throw OomphLibError("Couldn't find neighbouring element!",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    } // for (unsigned i=0;i<n_boundary0_element;i++)

    // Notify user
    oomph_info << "\nFinished enforcing periodicity!" << std::endl;

    // If we've got here then the periodicity has been set up
    Periodicity_has_been_enforced=true;
  } // if (!Periodicity_has_been_enforced)
} // End of enforce_time_periodic_boundary_conditions


//======start_of_set_neighbour_periodic_and_up_right_equivalents==========
/// Function to set the periodicity between two octrees and assign
/// the up and right equivalents to both elements. The arguments are the
/// current element, the (i,j,k) coordinates of the neighbouring element
/// and the direction to the neighbouring element.
//=========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::
set_neighbour_periodic_and_up_right_equivalents(FiniteElement* el0_pt,
    FiniteElement* el1_pt,
    const int& direction0)
{
  // The direction from the element on the other face/edge
  int direction1=0;

  // If we're on a face
  if ((direction0>=OcTreeNames::L)&&(direction0<=OcTreeNames::F))
  {
    // Calculate the reflected face direction
    direction1=OcTree::Reflect_face[direction0];
  }
  // If we're on a edge
  else if ((direction0>=OcTreeNames::LB)&&(direction0<=OcTreeNames::UF))
  {
    // Calculate the reflected edge direction
    direction1=OcTree::Reflect_edge[direction0];
  }
  // We should never reach here
  else
  {
    // Throw an error
    throw OomphLibError("Not looking for a face or edge, something's wrong!",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }

  // Get the tree root of the second element
  OcTreeRoot* octree0_pt=dynamic_cast<OcTreeRoot*>
                         (dynamic_cast<ELEMENT*>(el0_pt)->tree_pt()->root_pt());

  // Get the tree root of the second element
  OcTreeRoot* octree1_pt=dynamic_cast<OcTreeRoot*>
                         (dynamic_cast<ELEMENT*>(el1_pt)->tree_pt()->root_pt());

  //--------------------------------------------------------------------
  // Set the periodicity between the neighbouring OcTree objects and
  // set the up/right equivalents of the elements (it is assumed that
  // the elements are both oriented in the same way -- seems reasonable
  // for a uniform cube mesh)
  //--------------------------------------------------------------------
  // The first element (at t=0) will be periodic in the direction direction0
  octree0_pt->set_neighbour_periodic(direction0);

  // The second element (at t=1) will be periodic in the direction direction1
  octree1_pt->set_neighbour_periodic(direction1);

  // If we're on a face
  if ((direction0>=OcTreeNames::L)&&(direction0<=OcTreeNames::F))
  {
    // Pass the tree root of the neighbouring element
    octree0_pt->neighbour_pt(direction0)=octree1_pt;

    // Pass the tree root of the neighbouring element
    octree1_pt->neighbour_pt(direction1)=octree0_pt;
  }
  // If we're on an edge
  else if ((direction0>=OcTreeNames::LB)&&(direction0<=OcTreeNames::UF))
  {
    // Pass the tree root of the neighbouring element
    octree0_pt->add_edge_neighbour_pt(octree1_pt,direction0);

    // Pass the tree root of the neighbouring element
    octree1_pt->add_edge_neighbour_pt(octree0_pt,direction1);
  }

  //----------------------------------------------------------------
  // Now the appropriate OcTree objects are aware of each other,
  // they need to know their relative orientations. Specifically,
  // they need to know their up and right direction relative to each
  // other. This is painful to do by calling the function
  //               construct_up_right_equivalents()
  // from the Forest pointer of the mesh so just hardcode it (hacky,
  // I know...).
  //----------------------------------------------------------------
  // Set the up equivalent
  octree0_pt->set_up_equivalent(octree1_pt,OcTreeNames::U);

  // Set the up equivalent
  octree1_pt->set_up_equivalent(octree0_pt,OcTreeNames::U);

  // Set the right equivalent
  octree0_pt->set_right_equivalent(octree1_pt,OcTreeNames::R);

  // Set the right equivalent
  octree1_pt->set_right_equivalent(octree0_pt,OcTreeNames::R);
} // End of set_neighbour_periodic_and_up_right_equivalents


//======start_of_doc_solution=============================================
/// Document the solution
//=========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::doc_solution(const bool& doc_error)
{
  // Make an ofstream object to output the solution
  ofstream some_file;

  // Storage for the filename
  char filename[100];

  // The number of plot points
  unsigned n_plot_point=2;

  // The filename suffix added depending on the BCs we're using
  std::string filename_suffix="";

  // If we're solving the time-periodic problem
  if (GlobalParameters::Apply_time_periodic_boundary_conditions)
  {
    // Set the suffix
    filename_suffix="_pbc";
  }
  // If we're solving the initial-value problem
  else
  {
    // Set the suffix
    filename_suffix="_ic";
  }

  //----------------
  // Output solution
  //----------------
  // If we're solving the time-periodic problem
  if (GlobalParameters::Apply_time_periodic_boundary_conditions)
  {
    // Output the (numerically) approximated solution
    sprintf(filename,"%s/soln%s%i.dat",
            GlobalParameters::Doc_info.directory().c_str(),
            filename_suffix.c_str(),
            GlobalParameters::Doc_info.number());
  }
  // If we're solving the initial-value problem
  else
  {
    // Output the (numerically) approximated solution
    sprintf(filename,"%s/soln_ic%i.dat",
            GlobalParameters::Doc_info.directory().c_str(),
            GlobalParameters::Doc_info.number());
  }

  // Open the file
  some_file.open(filename);

  // Output the solution
  Bulk_mesh_pt->output(some_file,n_plot_point);

  // Now close the file
  some_file.close();

  //----------------------
  // Output exact solution
  //----------------------
  // Output the exact solution
  sprintf(filename,"%s/exact_soln%i.dat",
          GlobalParameters::Doc_info.directory().c_str(),
          GlobalParameters::Doc_info.number());

  // Open the file
  some_file.open(filename);

  // Dummy time value
  double dummy_time=0.0;

  // Output the exact solution
  Bulk_mesh_pt->output_fct(some_file,
                           n_plot_point,
                           dummy_time,
                           SinSolution::get_exact_u);

  // Now close the file
  some_file.close();

  //-------------------
  // Document the error
  //-------------------
  // Storage for the norm and the solution and the error
  double norm=0.0, error=0.0;

  // Output the error
  sprintf(filename,"%s/error%s%i.dat",
          GlobalParameters::Doc_info.directory().c_str(),
          filename_suffix.c_str(),
          GlobalParameters::Doc_info.number());

  // Open the file
  some_file.open(filename);

  // Compute the error in the solution
  Bulk_mesh_pt->compute_error(some_file,
                              SinSolution::get_exact_u,
                              dummy_time,
                              error,norm);

  // Now close the file
  some_file.close();

  //------------------------
  // Doc. solution and error
  //------------------------
  // Calculate the norm value we want
  norm=std::sqrt(norm);

  // Square root the error value
  error=std::sqrt(error);

  // Output the solution norm and error values
  oomph_info << "Solution norm : " << norm << std::endl;
  oomph_info << "Absolute error: " << error << std::endl;
  oomph_info << "Relative error: " << error/norm << std::endl;

  // Increment counter for solutions
  GlobalParameters::Doc_info.number()++;
} // End of doc_solution

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//======start_of_main=====================================================
/// Driver code for unsteady heat equation
//=========================================================================
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

  // Output directory
  GlobalParameters::Doc_info.set_directory("RESLT");

  // Start the timer
  double timer_s=TimingHelpers::timer();

  // Typedef the block element type
  typedef BlockPrecRefineableQUnsteadyHeatSpaceTimeElement<2,3> ELEMENT;

  // -------
  // TEST 1:
  // -------
  // Solve the problem using time-periodic boundary conditions.
  // NOTE: We can't just re-apply the boundary conditions for TEST 2 (where
  // we solve the initial-value problem) as applying periodic boundary
  // conditions involves copying nodes so we'd also have to regenerate
  // the mesh. It's easier just to create another Problem...
  {
    // Tell the user
    oomph_info << ANSIEscapeCode::Green
               << "\nTest 1: Solving with time-periodic boundary conditions!"
               << ANSIEscapeCode::Reset << std::endl;

    // Use time-periodic boundary conditions
    GlobalParameters::Apply_time_periodic_boundary_conditions=true;

    // Create a Problem pointer
    UnsteadyHeatProblem<ELEMENT> problem(&SinSolution::get_source);

    // Solve the space-time problem
    problem.newton_solve();

    // Document the solution
    problem.doc_solution();

    // Uniform mesh refinement
    problem.refine_uniformly();

    // Solve the space-time problem
    problem.newton_solve();

    // Document the solution
    problem.doc_solution();
  } // End of TEST 1

  // -------
  // TEST 2:
  // -------
  // Solve the problem using Dirichlet boundary conditions and an initial
  // condition.
  {
    // Tell the user
    oomph_info << ANSIEscapeCode::Green
               << "\nTest 2: Solving with an initial condition!"
               << ANSIEscapeCode::Reset << std::endl;

    // Reset the documentation counter
    GlobalParameters::Doc_info.number()=0;

    // Use time-periodic boundary conditions
    GlobalParameters::Apply_time_periodic_boundary_conditions=false;

    // Create a Problem pointer
    UnsteadyHeatProblem<ELEMENT> problem(&SinSolution::get_source);

    // Solve the space-time problem
    problem.newton_solve();

    // Document the solution
    problem.doc_solution();

    // Uniform mesh refinement
    problem.refine_uniformly();

    // Solve the space-time problem
    problem.newton_solve();

    // Document the solution
    problem.doc_solution();
  } // End of TEST 2

  // Tell the user we're done
  oomph_info << "\n3D space-time simulation complete!"
             << "\nTotal time for simulation [sec]: "
             << TimingHelpers::timer()-timer_s << "\n" << std::endl;
} // End of main
