// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
#include "lagrange_enforced_flow_preconditioner.h"

namespace oomph
{
  //==========================================================================
  /// Namespace for subsidiary preconditioner creation helper functions
  //==========================================================================
  namespace Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper
  {
    /// CG with diagonal preconditioner for W-block subsidiary linear
    /// systems.
    Preconditioner* get_w_cg_preconditioner()
    {
#ifdef OOMPH_HAS_TRILINOS
      InnerIterationPreconditioner<TrilinosAztecOOSolver,
                                   MatrixBasedDiagPreconditioner>* prec_pt =
        new InnerIterationPreconditioner<TrilinosAztecOOSolver,
                                         MatrixBasedDiagPreconditioner>;

      // Note: This makes CG a proper "inner iteration" for
      // which GMRES (may) no longer converge. We should really
      // use FGMRES or GMRESR for this. However, here the solver
      // is so good that it'll converge very quickly anyway
      // so there isn't much to be gained by limiting the number
      // of iterations...
      prec_pt->max_iter() = 4;
      prec_pt->solver_pt()->solver_type() = TrilinosAztecOOSolver::CG;
      prec_pt->solver_pt()->disable_doc_time();
      return prec_pt;
#else
      std::ostringstream err_msg;
      err_msg << "Inner CG preconditioner is unavailable.\n"
              << "Please install Trilinos.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    } // function get_w_cg_preconditioner

    /// Hypre Boomer AMG setting for the augmented momentum block
    /// of a 2D Navier-Stokes problem using the simple form of the viscous
    /// term (for serial code).
    Preconditioner* boomer_amg_for_2D_momentum_simple_visc()
    {
#ifdef OOMPH_HAS_HYPRE
      // Create a new HyprePreconditioner
      HyprePreconditioner* hypre_preconditioner_pt = new HyprePreconditioner;

      // Coarsening strategy
      // 1 = classical RS with no boundary treatment (not recommended in
      // parallel)
      hypre_preconditioner_pt->amg_coarsening() = 1;

      // Strength of dependence = 0.25
      hypre_preconditioner_pt->amg_strength() = 0.25;


      // Set the smoothers
      //   1 = Gauss-Seidel, sequential (very slow in parallel!)
      hypre_preconditioner_pt->amg_simple_smoother() = 1;

      // Set smoother damping (not required, so set to -1)
      hypre_preconditioner_pt->amg_damping() = -1;


      // Set number of cycles to 1xV(2,2)
      hypre_preconditioner_pt->set_amg_iterations(1);
      hypre_preconditioner_pt->amg_smoother_iterations() = 2;

      return hypre_preconditioner_pt;
#else
      std::ostringstream err_msg;
      err_msg << "hypre preconditioner is not available.\n"
              << "Please install Hypre.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    } // function boomer_amg_for_2D_momentum_simple_visc

    /// Hypre Boomer AMG setting for the augmented momentum block
    /// of a 2D Navier-Stokes problem using the stress divergence form of the
    /// viscous term (for serial code).
    Preconditioner* boomer_amg_for_2D_momentum_stressdiv_visc()
    {
#ifdef OOMPH_HAS_HYPRE
      // Create a new HyprePreconditioner
      HyprePreconditioner* hypre_preconditioner_pt = new HyprePreconditioner;

      // Coarsening strategy
      // 1 = classical RS with no boundary treatment (not recommended in
      // parallel)
      hypre_preconditioner_pt->amg_coarsening() = 1;

      // Strength of dependence = 0.668
      hypre_preconditioner_pt->amg_strength() = 0.668;


      // Set the smoothers
      //   1 = Gauss-Seidel, sequential (very slow in parallel!)
      hypre_preconditioner_pt->amg_simple_smoother() = 1;

      // Set smoother damping (not required, so set to -1)
      hypre_preconditioner_pt->amg_damping() = -1;


      // Set number of cycles to 1xV(2,2)
      hypre_preconditioner_pt->set_amg_iterations(1);
      hypre_preconditioner_pt->amg_smoother_iterations() = 2;

      return hypre_preconditioner_pt;
#else
      std::ostringstream err_msg;
      err_msg << "hypre preconditioner is not available.\n"
              << "Please install Hypre.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    } // function boomer_amg_for_2D_momentum_stressdiv_visc

    /// Hypre Boomer AMG setting for the augmented momentum block
    /// of a 3D Navier-Stokes problem (for serial code).
    Preconditioner* boomer_amg_for_3D_momentum()
    {
#ifdef OOMPH_HAS_HYPRE
      // Create a new HyprePreconditioner
      HyprePreconditioner* hypre_preconditioner_pt = new HyprePreconditioner;

      // Coarsening strategy
      // 1 = classical RS with no boundary treatment (not recommended in
      // parallel)
      hypre_preconditioner_pt->amg_coarsening() = 1;

      // Strength of dependence = 0.668
      hypre_preconditioner_pt->amg_strength() = 0.8;


      // Set the smoothers
      //   1 = Gauss-Seidel, sequential (very slow in parallel!)
      hypre_preconditioner_pt->amg_simple_smoother() = 1;

      // Set smoother damping (not required, so set to -1)
      hypre_preconditioner_pt->amg_damping() = -1;


      // Set number of cycles to 1xV(2,2)
      hypre_preconditioner_pt->set_amg_iterations(1);
      hypre_preconditioner_pt->amg_smoother_iterations() = 2;

      return hypre_preconditioner_pt;
#else
      std::ostringstream err_msg;
      err_msg << "hypre preconditioner is not available.\n"
              << "Please install Hypre.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    } // function boomer_amg_for_3D_momentum

    /// Hypre Boomer AMG setting for the augmented momentum block
    /// of a 3D Navier-Stokes problem (for serial code).
    Preconditioner* boomer_amg2v22_for_3D_momentum()
    {
#ifdef OOMPH_HAS_HYPRE
      // Create a new HyprePreconditioner
      HyprePreconditioner* hypre_preconditioner_pt = new HyprePreconditioner;

      // Coarsening strategy
      // 1 = classical RS with no boundary treatment (not recommended in
      // parallel)
      hypre_preconditioner_pt->amg_coarsening() = 1;

      // Strength of dependence = 0.668
      hypre_preconditioner_pt->amg_strength() = 0.8;


      // Set the smoothers
      //   1 = Gauss-Seidel, sequential (very slow in parallel!)
      hypre_preconditioner_pt->amg_simple_smoother() = 1;

      // Set smoother damping (not required, so set to -1)
      hypre_preconditioner_pt->amg_damping() = -1;


      // Set number of cycles to 1xV(2,2)
      hypre_preconditioner_pt->set_amg_iterations(2);
      hypre_preconditioner_pt->amg_smoother_iterations() = 2;

      return hypre_preconditioner_pt;
#else
      std::ostringstream err_msg;
      err_msg << "hypre preconditioner is not available.\n"
              << "Please install Hypre.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    } // function boomer_amg_for_3D_momentum


    /// Hypre Boomer AMG setting for the 2D Poisson problem
    /// (for serial code).
    Preconditioner* boomer_amg_for_2D_poisson_problem()
    {
#ifdef OOMPH_HAS_HYPRE
      // Create a new HyprePreconditioner
      HyprePreconditioner* hypre_preconditioner_pt = new HyprePreconditioner;

      // Coarsening strategy
      // 1 = classical RS with no boundary treatment (not recommended in
      // parallel)
      hypre_preconditioner_pt->amg_coarsening() = 1;

      // Strength of dependence = 0.25
      hypre_preconditioner_pt->amg_strength() = 0.25;


      // Set the smoothers
      //   0 = Jacobi
      hypre_preconditioner_pt->amg_simple_smoother() = 0;

      // Set Jacobi damping = 2/3
      hypre_preconditioner_pt->amg_damping() = 0.668;


      // Set number of cycles to 1xV(2,2)
      hypre_preconditioner_pt->set_amg_iterations(2);
      hypre_preconditioner_pt->amg_smoother_iterations() = 1;

      return hypre_preconditioner_pt;
#else
      std::ostringstream err_msg;
      err_msg << "hypre preconditioner is not available.\n"
              << "Please install Hypre.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    } // function boomer_amg_for_2D_poisson_problem

    /// Hypre Boomer AMG setting for the 3D Poisson problem
    /// (for serial code).
    Preconditioner* boomer_amg_for_3D_poisson_problem()
    {
#ifdef OOMPH_HAS_HYPRE
      // Create a new HyprePreconditioner
      HyprePreconditioner* hypre_preconditioner_pt = new HyprePreconditioner;

      // Coarsening strategy
      // 1 = classical RS with no boundary treatment (not recommended in
      // parallel)
      hypre_preconditioner_pt->amg_coarsening() = 1;

      // Strength of dependence = 0.7
      hypre_preconditioner_pt->amg_strength() = 0.7;


      // Set the smoothers
      //   0 = Jacobi
      hypre_preconditioner_pt->amg_simple_smoother() = 0;

      // Set smoother damping = 2/3
      hypre_preconditioner_pt->amg_damping() = 0.668;


      // Set number of cycles to 2xV(1,1)
      hypre_preconditioner_pt->set_amg_iterations(2);
      hypre_preconditioner_pt->amg_smoother_iterations() = 1;

      return hypre_preconditioner_pt;
#else
      std::ostringstream err_msg;
      err_msg << "hypre preconditioner is not available.\n"
              << "Please install Hypre.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    } // function boomer_amg_for_3D_poisson_problem

  } // namespace
    // Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper

  /// /////////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////////

  /// Apply the preconditioner.
  /// r is the residual (rhs), z will contain the solution.
  void LagrangeEnforcedFlowPreconditioner::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
  {
#ifdef PARANOID
    if (Preconditioner_has_been_setup == false)
    {
      std::ostringstream error_message;
      error_message << "setup() must be called before using "
                    << "preconditioner_solve()";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (z.built())
    {
      if (z.nrow() != r.nrow())
      {
        std::ostringstream error_message;
        error_message << "The vectors z and r must have the same number of "
                      << "of global rows";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // if z is not setup then give it the same distribution
    if (!z.distribution_pt()->built())
    {
      z.build(r.distribution_pt(), 0.0);
    }

    // Working vectors.
    DoubleVector temp_vec;
    DoubleVector another_temp_vec;
    DoubleVector yet_another_temp_vec;


    // -----------------------------------------------------------------------
    // Step 1 - apply approximate W block inverse to Lagrange multiplier
    // unknowns
    // -----------------------------------------------------------------------
    // For each subsystem associated with each Lagrange multiplier, we loop
    // through and:
    // 1) extract the block entries from r
    // 2) apply the inverse
    // 3) return the entries to z.
    for (unsigned l_i = 0; l_i < N_lagrange_doftypes; l_i++)
    {
      // The Lagrange multiplier block type.
      const unsigned l_ii = N_fluid_doftypes + l_i;

      // Extract the block
      this->get_block_vector(l_ii, r, temp_vec);

      // Apply the inverse.
      const unsigned vec_nrow_local = temp_vec.nrow_local();
      double* vec_values_pt = temp_vec.values_pt();
      for (unsigned i = 0; i < vec_nrow_local; i++)
      {
        vec_values_pt[i] = vec_values_pt[i] * Inv_w_diag_values[l_i][i];
      } // for

      // Return the unknowns
      this->return_block_vector(l_ii, temp_vec, z);

      // Clear vectors.
      temp_vec.clear();
    } // for

    // -----------------------------------------------------------------------
    // Step 2 - apply the augmented Navier-Stokes matrix inverse to the
    // velocity and pressure unknowns
    // -----------------------------------------------------------------------

    // At this point, all vectors are cleared.
    if (Using_superlu_ns_preconditioner)
    {
      // Which block types corresponds to the fluid block types.
      Vector<unsigned> fluid_block_indices(N_fluid_doftypes, 0);
      for (unsigned b = 0; b < N_fluid_doftypes; b++)
      {
        fluid_block_indices[b] = b;
      }

      this->get_concatenated_block_vector(fluid_block_indices, r, temp_vec);

      // temp_vec contains the (concatenated) fluid rhs.
      Navier_stokes_preconditioner_pt->preconditioner_solve(temp_vec,
                                                            another_temp_vec);

      temp_vec.clear();

      // Return it to the unknowns.
      this->return_concatenated_block_vector(
        fluid_block_indices, another_temp_vec, z);

      another_temp_vec.clear();
    }
    else
    {
      // The Navier-Stokes preconditioner is a block preconditioner.
      // Thus is handles all of the block vector extraction and returns.
      Navier_stokes_preconditioner_pt->preconditioner_solve(r, z);
    }
  } // end of preconditioner_solve

  /// Set the meshes,
  /// the first mesh in the vector must be the bulk mesh.
  void LagrangeEnforcedFlowPreconditioner::set_meshes(
    const Vector<Mesh*>& mesh_pt)
  {
    // There should be at least two meshes passed to this preconditioner.
    const unsigned nmesh = mesh_pt.size();

#ifdef PARANOID
    if (nmesh < 2)
    {
      std::ostringstream err_msg;
      err_msg << "There should be at least two meshes.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check that all pointers are not null
    for (unsigned mesh_i = 0; mesh_i < nmesh; mesh_i++)
    {
      if (mesh_pt[mesh_i] == 0)
      {
        std::ostringstream err_msg;
        err_msg << "The pointer mesh_pt[" << mesh_i << "] is null.\n";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }

    // We assume that the first mesh is the Navier-Stokes "bulk" mesh.
    // To check this, the elemental dimension must be the same as the
    // nodal (spatial) dimension.
    //
    // We store the elemental dimension i.e. the number of local coordinates
    // required to parametrise its geometry.
    const unsigned elemental_dim = mesh_pt[0]->elemental_dimension();

    // The dimension of the nodes in the first element in the (supposedly)
    // bulk mesh.
    const unsigned nodal_dim = mesh_pt[0]->nodal_dimension();

    // Check if the first mesh is the "bulk" mesh.
    // Here we assume only one mesh contains "bulk" elements.
    // All subsequent meshes contain block preconditionable elements which
    // re-classify the bulk velocity DOFs to constrained velocity DOFs.
    if (elemental_dim != nodal_dim)
    {
      std::ostringstream err_msg;
      err_msg << "In the first mesh, the elements have elemental dimension "
              << "of " << elemental_dim << ",\n"
              << "with a nodal dimension of " << nodal_dim << ".\n"
              << "The first mesh does not contain 'bulk' elements.\n"
              << "Please re-order your mesh_pt vector.\n";

      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Set the number of meshes
    this->set_nmesh(nmesh);

    // Set the meshes
    for (unsigned mesh_i = 0; mesh_i < nmesh; mesh_i++)
    {
      this->set_mesh(mesh_i, mesh_pt[mesh_i]);
    }

    // We also store the meshes and number of meshes locally in this class.
    // This may seem slightly redundant, since we always have all the meshes
    // stored in the upper most master block preconditioner.
    // But at the moment there is no mapping/look up scheme between
    // master and subsidiary block preconditioners for meshes.
    //
    // If this is a subsidiary block preconditioner, we don't know which of
    // the master's meshes belong to us. We need this information to set up
    // look up lists in the function setup(...).
    // Thus we store them local to this class.
    My_mesh_pt = mesh_pt;
    My_nmesh = nmesh;
  } // function set_meshes


  //==========================================================================
  /// Setup the Lagrange enforced flow preconditioner. This
  /// extracts blocks corresponding to the velocity and Lagrange multiplier
  /// unknowns, creates the matrices actually needed in the application of the
  /// preconditioner and deletes what can be deleted... Note that
  /// this preconditioner needs a CRDoubleMatrix.
  //==========================================================================
  void LagrangeEnforcedFlowPreconditioner::setup()
  {
    // clean
    this->clean_up_memory();

#ifdef PARANOID
    // Paranoid check that meshes have been set. In this preconditioner, we
    // always have to set meshes.
    if (My_nmesh == 0)
    {
      std::ostringstream err_msg;
      err_msg << "There are no meshes set. Please call set_meshes(...)\n"
              << "with at least two mesh.";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // -----------------------------------------------------------------------
    // Step 1 - Construct the dof_to_block_map vector.
    // -----------------------------------------------------------------------
    // Assumption: The first mesh is always the "bulk" mesh
    // (Navier-Stokes mesh), which contains block preconditionable
    // Navier-Stokes elements. Thus the bulk elements classify the velocity
    // and pressure degrees of freedom (ndof_types = 3(4) in 2(3)D).
    // All subsequent meshes contain the constrained velocity DOF types,
    // then the Lagrange multiplier DOF types.
    //
    // Thus, a general ordering of DOF types (in 3D) follows the ordering:
    //
    //  0 1 2 3   4 5 6 7  8  ..x   x+0 x+1 x+2 x+3 x+4
    // [u v w p] [u v w l1 l2 ...] [u   v   w   l1  l2 ...] ...
    //
    // where the square brackets [] represent the DOF types in each mesh.
    //
    // Example:
    // Consider the case of imposing parallel outflow (3 constrained velocity
    // DOF types and 2 Lagrange multiplier DOF types) and tangential flow (3
    // constrained velocity DOF types and 1 Lagrange multiplier DOF type)
    // along two different boundaries in 3D. The resulting natural ordering of
    // the DOF types is:
    //
    // [0 1 2 3] [4  5  6   7   8 ] [9  10 11 12 ]
    // [u v w p] [up vp wp Lp1 Lp2] [ut vt wt Lt1]
    //
    //
    // In our implementation, the desired block structure is:
    // | u v w | up vp wp | ut vt wt | p | Lp1 Lp2 Lt1 |
    //
    // The dof_to_block_map should have the following construction:
    //
    //    dof_to_block_map[$dof_number] = $block_number
    //
    // Thus, the dof_to_block_map for the example above should be:
    //
    // To achieve this, we use the dof_to_block_map:
    // dof_to_block_map = [0 1 2 9 3  4  5  10  11  6  7  8  12]
    //
    // To generalise the construction of the dof_to_block_map vector
    // (to work for any number of meshes), we first require some auxiliary
    // variables to aid us in this endeavour.
    // -----------------------------------------------------------------------

    // Set up the My_ndof_types_in_mesh vector.
    // If this is already constructed, we reuse it instead.
    if (My_ndof_types_in_mesh.size() == 0)
    {
      for (unsigned mesh_i = 0; mesh_i < My_nmesh; mesh_i++)
      {
        My_ndof_types_in_mesh.push_back(My_mesh_pt[mesh_i]->ndof_types());
      }
    }

    // Get the spatial dimension of the problem.
    unsigned spatial_dim = My_mesh_pt[0]->nodal_dimension();

    // Get the number of DOF types.
    unsigned n_dof_types = ndof_types();

#ifdef PARANOID
    // We cannot check which DOF types are "correct" in the sense that there
    // is a distinction between bulk and constrained velocity DOF tyoes.
    // But we can at least check if the ndof_types matches the total number
    // of DOF types from the meshes.
    //
    // This check is not necessary for the upper most master block
    // preconditioner since the ndof_types() is calculated by looping
    // through the meshes!
    //
    // This check is useful if this is a subsidiary block preconditioner and
    // incorrect DOF type merging has taken place.
    if (is_subsidiary_block_preconditioner())
    {
      unsigned tmp_ndof_types = 0;
      for (unsigned mesh_i = 0; mesh_i < My_nmesh; mesh_i++)
      {
        tmp_ndof_types += My_ndof_types_in_mesh[mesh_i];
      }

      if (tmp_ndof_types != n_dof_types)
      {
        std::ostringstream err_msg;
        err_msg << "The number of DOF types are incorrect.\n"
                << "The total DOF types from the meshes is: " << tmp_ndof_types
                << ".\n"
                << "The number of DOF types from "
                << "BlockPreconditioner::ndof_types() is " << n_dof_types
                << ".\n";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // The number of velocity DOF types: We assume that all meshes classify
    // at least some of the velocity DOF types (bulk/constrained), thus the
    // total number of velocity DOF types is the spatial dimension multiplied
    // by the number of meshes.
    N_velocity_doftypes = My_nmesh * spatial_dim;

    // Fluid has +1 for the pressure.
    N_fluid_doftypes = N_velocity_doftypes + 1;

    // The rest are Lagrange multiplier DOF types.
    N_lagrange_doftypes = n_dof_types - N_fluid_doftypes;

    /// ///////////////////////////////////////////////////////////
    //   Now construct the DOF to block map for block_setup()   //
    /// ///////////////////////////////////////////////////////////

    // Now that we have
    //
    // (1) spatial dimension of the problem and
    // (2) the number of DOF types in each of the meshes.
    //
    // We observe that the problem dimension is 3 and
    // My_ndof_type_in_mesh = [4, 5, 4].
    //
    // With these information we can construct the desired block structure:
    // | u v w | up vp wp | ut vt wt | p | Lp1 Lp2 Lt1 |
    //
    // The block structure is determined by the vector dof_to_block_map we
    // give to the function block_setup(...).

    // This preconditioner permutes the DOF numbers to get the block numbers.
    // I.e. nblock_types = ndof_types, but they are re-ordered
    // so that we have (in order):
    // 1) bulk velocity DOF types
    // 2) constrained velocity DOF types
    // 3) pressure DOF type
    // 4) Lagrange multiplier DOF types
    //
    // Recall that the natural ordering of the DOF types are ordered first by
    // their meshes, then the elemental DOF type as described by the
    // function get_dof_numbers_for_unknowns().
    //
    // Also recall that (for this problem), we assume/require that every mesh
    // (re-)classify at least some of the velocity DOF types, furthermore,
    // the velocity DOF type classification comes before the
    // pressure / Lagrange multiplier DOF types.
    //
    // Consider the same example as shown previously, repeated here for your
    // convenience:
    //
    // Natural DOF type ordering:
    //  0 1 2 3   4  5  6   7   8    9  10 11 12   <- DOF number.
    // [u v w p] [up vp wp Lp1 Lp2] [ut vt wt Lt1] <- DOF type.
    //
    // Desired block structure:
    //  0 1 2   3  4  5    6  7  8     9   10  11  12   <- block number
    // [u v w | up vp wp | ut vt wt ] [p | Lp1 Lp2 Lt1] <- DOF type
    //
    // dof_to_block_map[$dof_number] = $block_number
    //
    // Required dof_to_block_map:
    //  0 1 2 9 3  4  5  10  11  6  7  8  12   <- block number
    // [u v w p up vp wp Lp1 Lp2 ut vt wt Lt1] <- DOF type
    //
    // Consider the 4th entry of the dof_to_block_map, this represents the
    // pressure DOF type, from the desired block structure we see this takes
    // the block number 9.
    //
    // One way to generalise the construction of the dof_to_block_map  is to
    // lump the first $spatial_dimension number of DOF types
    // from each mesh together, then lump the remaining DOF types together.
    //
    // Notice that the values in the velocity DOF types of the
    // dof_to_block_map vector increases sequentially, from 0 to 8. The values
    // of the Lagrange multiplier DOF types follow the same pattern, from 9 to
    // 12. We follow this construction.

    // Storage for the dof_to_block_map vector.
    Vector<unsigned> dof_to_block_map(n_dof_types, 0);

    // Index for the dof_to_block_map vector.
    unsigned temp_index = 0;

    // Value for the velocity DOF type.
    unsigned velocity_number = 0;

    // Value for the pressure/Lagrange multiplier DOF type.
    unsigned pres_lgr_number = N_velocity_doftypes;

    // Loop through the number of meshes.
    for (unsigned mesh_i = 0; mesh_i < My_nmesh; mesh_i++)
    {
      // Fill in the velocity DOF types.
      for (unsigned dim_i = 0; dim_i < spatial_dim; dim_i++)
      {
        dof_to_block_map[temp_index++] = velocity_number++;
      } // for

      // Loop through the pressure/Lagrange multiplier DOF types.
      unsigned ndof_type_in_mesh_i = My_ndof_types_in_mesh[mesh_i];
      for (unsigned doftype_i = spatial_dim; doftype_i < ndof_type_in_mesh_i;
           doftype_i++)
      {
        dof_to_block_map[temp_index++] = pres_lgr_number++;
      } // for
    } // for

    // Call the block setup
    this->block_setup(dof_to_block_map);


    // -----------------------------------------------------------------------
    // Step 2 - Get the infinite norm of Navier-Stokes F block.
    // -----------------------------------------------------------------------

    // Extract the velocity blocks.
    DenseMatrix<CRDoubleMatrix*> v_aug_pt(
      N_velocity_doftypes, N_velocity_doftypes, 0);

    for (unsigned row_i = 0; row_i < N_velocity_doftypes; row_i++)
    {
      for (unsigned col_i = 0; col_i < N_velocity_doftypes; col_i++)
      {
        v_aug_pt(row_i, col_i) = new CRDoubleMatrix;
        this->get_block(row_i, col_i, *v_aug_pt(row_i, col_i));
      } // for
    } // for

    // Now get the infinite norm.
    if (Use_norm_f_for_scaling_sigma)
    {
      Scaling_sigma = -CRDoubleMatrixHelpers::inf_norm(v_aug_pt);
    } // if

#ifdef PARANOID
    // Warning for division by zero.
    if (Scaling_sigma == 0.0)
    {
      std::ostringstream warning_stream;
      warning_stream << "WARNING: " << std::endl
                     << "The scaling (Scaling_sigma) is " << Scaling_sigma
                     << ".\n"
                     << "Division by 0 will occur."
                     << "\n";
      OomphLibWarning(
        warning_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (Scaling_sigma > 0.0)
    {
      std::ostringstream warning_stream;
      warning_stream << "WARNING: " << std::endl
                     << "The scaling (Scaling_sigma) is positive: "
                     << Scaling_sigma << std::endl
                     << "Performance may be degraded." << std::endl;
      OomphLibWarning(
        warning_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // -----------------------------------------------------------------------
    // Step 3 - Augment the constrained fluid blocks.
    // -----------------------------------------------------------------------
    // Loop through the Lagrange multipliers and do three things:
    // For each Lagrange block:
    //   3.1) Extract the mass matrices and store the location of non-zero mass
    //        matrices.
    //
    //   3.2) a) Create inv_w (for the augmentation)
    //        b) Store the diagonal values of inv_w (for preconditioner solve)
    //
    //   3.3) Perform the augmentation: v_aug + m_i * inv(w_i) * m_j
    //
    //   3.4) Clean up memory.

    // Storage for the inv w diag values.
    Inv_w_diag_values.clear();
    for (unsigned l_i = 0; l_i < N_lagrange_doftypes; l_i++)
    {
      // Step 3.1: Location and extraction of non-zero mass matrices.

      // Storage for the current Lagrange block mass matrices.
      Vector<CRDoubleMatrix*> mm_pt;
      Vector<CRDoubleMatrix*> mmt_pt;

      // Block type for the Lagrange multiplier.
      const unsigned lgr_block_num = N_fluid_doftypes + l_i;

      // Store the mass matrix locations for the current Lagrange block.
      Vector<unsigned> mm_locations;

      // Store the number of mass matrices.
      unsigned n_mm = 0;

      // Go along the block columns for the current Lagrange block ROW.
      for (unsigned col_i = 0; col_i < N_velocity_doftypes; col_i++)
      {
        // Get the block matrix for this block column.
        CRDoubleMatrix* mm_temp_pt = new CRDoubleMatrix;
        this->get_block(lgr_block_num, col_i, *mm_temp_pt);

        // Check if this is non-zero
        if (mm_temp_pt->nnz() > 0)
        {
          mm_locations.push_back(col_i);
          mm_pt.push_back(mm_temp_pt);
          n_mm++;
        }
        else
        {
          // This is just an empty matrix. No need to keep this.
          delete mm_temp_pt;
        }
      } // loop through the columns of the Lagrange row.

#ifdef PARANOID
      if (n_mm == 0)
      {
        std::ostringstream warning_stream;
        warning_stream << "WARNING:\n"
                       << "There are no mass matrices on Lagrange block row "
                       << l_i << ".\n"
                       << "Perhaps the problem setup is incorrect."
                       << std::endl;
        OomphLibWarning(warning_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Get the transpose of the mass matrices.
      for (unsigned mm_i = 0; mm_i < n_mm; mm_i++)
      {
        // Get the block matrix for this block column.
        CRDoubleMatrix* mm_temp_pt = new CRDoubleMatrix;
        this->get_block(mm_locations[mm_i], lgr_block_num, *mm_temp_pt);

        if (mm_temp_pt->nnz() > 0)
        {
          mmt_pt.push_back(mm_temp_pt);
        }
        else
        {
          // There should be a non-zero mass matrix here, since L=(L^T)^T
#ifdef PARANOID
          {
            std::ostringstream warning_stream;
            warning_stream << "WARNING:\n"
                           << "The mass matrix block " << mm_locations[mm_i]
                           << " in L^T block " << l_i << " is zero.\n"
                           << "Perhaps the problem setup is incorrect."
                           << std::endl;
            OomphLibWarning(warning_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
          }
#endif
        }
      } // loop through the ROW of the Lagrange COLUMN.

      //  Step 3.2: Create inv_w and store its diagonal entries.

      // Storage for inv_w
      CRDoubleMatrix* inv_w_pt =
        new CRDoubleMatrix(this->Block_distribution_pt[lgr_block_num]);

      // Get the number of local rows for this Lagrange block.
      unsigned long l_i_nrow_local =
        this->Block_distribution_pt[lgr_block_num]->nrow_local();

      // The first row, for the column offset (required in parallel).
      unsigned l_i_first_row =
        this->Block_distribution_pt[lgr_block_num]->first_row();

      // A vector to contain the results of mass matrices squared.
      Vector<double> w_i_diag_values(l_i_nrow_local, 0);

      // Create mm*mm^T, and component-wise add the mass matrices
      for (unsigned m_i = 0; m_i < n_mm; m_i++)
      {
        // Create mm*mm^T
        CRDoubleMatrix temp_mm_sqrd;
        temp_mm_sqrd.build(mm_pt[m_i]->distribution_pt());
        mm_pt[m_i]->multiply(*mmt_pt[m_i], temp_mm_sqrd);

        // Extract diagonal entries
        Vector<double> m_diag = temp_mm_sqrd.diagonal_entries();

        // Loop through the entries, add them.
        for (unsigned long row_i = 0; row_i < l_i_nrow_local; row_i++)
        {
          w_i_diag_values[row_i] += m_diag[row_i];
        }
      }

      // Storage for inv_w matrix vectors
      Vector<double> invw_i_diag_values(l_i_nrow_local, 0);
      Vector<int> w_i_column_indices(l_i_nrow_local);
      Vector<int> w_i_row_start(l_i_nrow_local + 1);

      // Divide by Scaling_sigma and create the inverse of w.
      for (unsigned long row_i = 0; row_i < l_i_nrow_local; row_i++)
      {
        invw_i_diag_values[row_i] = Scaling_sigma / w_i_diag_values[row_i];

        w_i_column_indices[row_i] = row_i + l_i_first_row;
        w_i_row_start[row_i] = row_i;
      }

      w_i_row_start[l_i_nrow_local] = l_i_nrow_local;

      // Theses are square matrices. So we can use the l_i_nrow_global for the
      // number of columns.
      unsigned long l_i_nrow_global =
        this->Block_distribution_pt[lgr_block_num]->nrow();
      inv_w_pt->build(
        l_i_nrow_global, invw_i_diag_values, w_i_column_indices, w_i_row_start);

      Inv_w_diag_values.push_back(invw_i_diag_values);


      // Step 3.3: Perform the augmentation: v_aug + m_i * inv(w_i) * m_j

      /// /////////////////////////////////////////////////////////////////////
      // Now we create the augmented matrix in v_aug_pt.
      // v_aug_pt is already re-ordered
      // Loop through the mm_locations
      for (unsigned ii = 0; ii < n_mm; ii++)
      {
        // Location of the mass matrix
        unsigned aug_i = mm_locations[ii];

        for (unsigned jj = 0; jj < n_mm; jj++)
        {
          // Location of the mass matrix
          unsigned aug_j = mm_locations[jj];

          // Storage for intermediate results.
          CRDoubleMatrix aug_block;
          CRDoubleMatrix another_aug_block;

          // aug_block = mmt*inv_w
          mmt_pt[ii]->multiply((*inv_w_pt), (aug_block));

          // another_aug_block = aug_block*mm = mmt*inv_w*mm
          aug_block.multiply(*mm_pt[jj], another_aug_block);

          // Add the augmentation.
          v_aug_pt(aug_i, aug_j)
            ->add(another_aug_block, *v_aug_pt(aug_i, aug_j));
        } // loop jj
      } // loop ii

      // Step 3.4: Clean up memory.
      delete inv_w_pt;
      inv_w_pt = 0;
      for (unsigned m_i = 0; m_i < n_mm; m_i++)
      {
        delete mm_pt[m_i];
        mm_pt[m_i] = 0;
        delete mmt_pt[m_i];
        mmt_pt[m_i] = 0;
      }
    } // loop through Lagrange multipliers.

    // -----------------------------------------------------------------------
    // Step 4 - Replace all the velocity blocks
    // -----------------------------------------------------------------------
    // When we replace (DOF) blocks, the indices have to respect the DOF type
    // ordering. This is because only DOF type information is passed between
    // preconditioners. So we need to create the inverse of the
    // dof_to_block_map... a block_to_dof_map!
    //
    // Note: Once the DOF blocks have been replaced, further calls to
    // get_block at this preconditioning hierarchy level and/or lower will use
    // the (nearest) replaced blocks.
    Vector<unsigned> block_to_dof_map(n_dof_types, 0);
    for (unsigned dof_i = 0; dof_i < n_dof_types; dof_i++)
    {
      block_to_dof_map[dof_to_block_map[dof_i]] = dof_i;
    }

    // Now do the replacement of all blocks in v_aug_pt
    for (unsigned block_row_i = 0; block_row_i < N_velocity_doftypes;
         block_row_i++)
    {
      unsigned temp_dof_row_i = block_to_dof_map[block_row_i];
      for (unsigned block_col_i = 0; block_col_i < N_velocity_doftypes;
           block_col_i++)
      {
        unsigned temp_dof_col_i = block_to_dof_map[block_col_i];
        this->set_replacement_dof_block(
          temp_dof_row_i, temp_dof_col_i, v_aug_pt(block_row_i, block_col_i));
      }
    }

    // -----------------------------------------------------------------------
    // Step 5 - Set up Navier-Stokes preconditioner
    // If the subsidiary fluid preconditioner is a block preconditioner:
    //   5.1) Set up the dof_number_in_master_map
    //   5.2) Set up the doftype_coarsen_map
    // otherwise:
    //   5.3) Concatenate the fluid matrices into one big matrix and pass it
    //        to the solver.
    // -----------------------------------------------------------------------

    // First we determine if we're using a block preconditioner or not.
    BlockPreconditioner<CRDoubleMatrix>* navier_stokes_block_preconditioner_pt =
      dynamic_cast<BlockPreconditioner<CRDoubleMatrix>*>(
        Navier_stokes_preconditioner_pt);
    Navier_stokes_preconditioner_is_block_preconditioner = true;
    if (navier_stokes_block_preconditioner_pt == 0)
    {
      Navier_stokes_preconditioner_is_block_preconditioner = false;
    }

    // If the Navier-Stokes preconditioner is a block preconditioner, then we
    // need to turn it into a subsidiary block preconditioner.
    if (Navier_stokes_preconditioner_is_block_preconditioner)
    {
      // Assumption: All Navier-Stokes block preconditioners take $dim
      // number of velocity DOF types and 1 pressure DOF type.

      // Step 5.1: Set up the dof_number_in_master_map

      // First we create the dof_number_in_master_map vector to pass to the
      // subsidiary preconditioner's
      // turn_into_subsidiary_block_preconditioner(...) function.
      //
      // This vector maps the subsidiary DOF numbers and the master DOF
      // numbers and has the construction:
      //
      //   dof_number_in_master_map[subsidiary DOF number] = master DOF number
      //
      // Example: Using the example above, our problem has the natural
      // DOF type ordering:
      //
      //  0 1 2 3   4  5  6   7   8    9  10 11 12   <- DOF number in master
      // [u v w p] [up vp wp Lp1 Lp2] [ut vt wt Lt1] <- DOF type
      //
      // For now, we ignore the fact that this preconditioner's number of
      // DOF types may be more fine grain than what is assumed by the
      // subsidiary block preconditioner. Normally, the order of the values in
      // dof_number_in_master_map matters, since the indices must match the
      // DOF type in the subsidiary block preconditioner (see the assumption
      // above). For example, if the (subsidiary) LSC block preconditioner
      // requires the DOF type ordering:
      //
      // [0 1 2 3]
      //  u v w p
      //
      // Then the dof_number_in_master_map vector must match the u velocity
      // DOF type in the subsidiary preconditioner with the u velocity in the
      // master preconditioner, etc...
      //
      // However, we shall see (later) that it does not matter in this
      // instance because the DOF type coarsening feature overrides the
      // ordering provided here.
      //
      // For now, we only need to ensure that the subsidiary preconditioner
      // knows about 10 master DOF types (9 velocity and 1 pressure), the
      // ordering does not matter.
      //
      // We pass to the subsidiary block preconditioner the following
      // dof_number_in_master_map:
      // [0 1 2 4  5  6  9  10 11 3] <- DOF number in master
      //  u v w up vp wp ut vt wt p  <- corresponding DOF type

      // Variable to keep track of the DOF number.
      unsigned temp_dof_number = 0;

      // Storage for the dof_number_in_master_map
      Vector<unsigned> dof_number_in_master_map;
      for (unsigned mesh_i = 0; mesh_i < My_nmesh; mesh_i++)
      {
        // Store the velocity dof types.
        for (unsigned dim_i = 0; dim_i < spatial_dim; dim_i++)
        {
          dof_number_in_master_map.push_back(temp_dof_number + dim_i);
        } // for spatial_dim

        // Update the DOF index
        temp_dof_number += My_ndof_types_in_mesh[mesh_i];
      } // for My_nmesh

      // Push back the pressure DOF type
      dof_number_in_master_map.push_back(spatial_dim);


      // Step 5.2 DOF type coarsening.

      // Since this preconditioner works with more fine grained DOF types than
      // the subsidiary block preconditioner, we need to tell the subsidiary
      // preconditioner which DOF types to coarsen together.
      //
      // The LSC preconditioner expects 2(3) velocity dof types, and 1
      // pressure DOF types. Thus, we give it this list:
      // u [0, 3, 6]
      // v [1, 4, 7]
      // w [2, 5, 8]
      // p [9]
      //
      // See how the ordering of dof_number_in_master_map (constructed above)
      // does not matter as long as we construct the doftype_coarsen_map
      // correctly.

      // Storage for which subsidiary DOF types to coarsen
      Vector<Vector<unsigned>> doftype_coarsen_map;
      for (unsigned direction = 0; direction < spatial_dim; direction++)
      {
        Vector<unsigned> dir_doftypes_vec(My_nmesh, 0);
        for (unsigned mesh_i = 0; mesh_i < My_nmesh; mesh_i++)
        {
          dir_doftypes_vec[mesh_i] = spatial_dim * mesh_i + direction;
        }
        doftype_coarsen_map.push_back(dir_doftypes_vec);
      }

      Vector<unsigned> ns_p_vec(1, 0);

      // This is simply the number of velocity dof types,
      ns_p_vec[0] = My_nmesh * spatial_dim;

      doftype_coarsen_map.push_back(ns_p_vec);

      // Turn the Navier-Stokes block preconditioner into a subsidiary block
      // preconditioner.
      navier_stokes_block_preconditioner_pt
        ->turn_into_subsidiary_block_preconditioner(
          this, dof_number_in_master_map, doftype_coarsen_map);

      // Call the setup function
      navier_stokes_block_preconditioner_pt->setup(matrix_pt());
    }
    else
    {
      // Step 5.3: This is not a block preconditioner, thus we need to
      // concatenate all the fluid matrices and pass them to the solver.

      // Select all the fluid blocks (velocity and pressure)
      VectorMatrix<BlockSelector> f_aug_blocks(N_fluid_doftypes,
                                               N_fluid_doftypes);
      for (unsigned block_i = 0; block_i < N_fluid_doftypes; block_i++)
      {
        for (unsigned block_j = 0; block_j < N_fluid_doftypes; block_j++)
        {
          f_aug_blocks[block_i][block_j].select_block(block_i, block_j, true);
        }
      }

      // Note: This will use the replaced blocks.
      CRDoubleMatrix f_aug_block = this->get_concatenated_block(f_aug_blocks);

      if (Using_superlu_ns_preconditioner)
      {
        if (Navier_stokes_preconditioner_pt == 0)
        {
          Navier_stokes_preconditioner_pt = new SuperLUPreconditioner;
        }
      }
      else
      {
        if (Navier_stokes_preconditioner_pt == 0)
        {
          std::ostringstream err_msg;
          err_msg << "Not using SuperLUPreconditioner for NS block,\n"
                  << "but the Navier_stokes_preconditioner_pt is null.\n";
          throw OomphLibError(
            err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Setup the solver.
      Navier_stokes_preconditioner_pt->setup(&f_aug_block);

    } // else

    // Clean up memory
    const unsigned v_aug_nrow = v_aug_pt.nrow();
    const unsigned v_aug_ncol = v_aug_pt.ncol();
    for (unsigned v_row = 0; v_row < v_aug_nrow; v_row++)
    {
      for (unsigned v_col = 0; v_col < v_aug_ncol; v_col++)
      {
        delete v_aug_pt(v_row, v_col);
        v_aug_pt(v_row, v_col) = 0;
      }
    }

    Preconditioner_has_been_setup = true;
  } // func setup

  /// Function to set a new momentum matrix preconditioner
  /// (inexact solver)
  void LagrangeEnforcedFlowPreconditioner::set_navier_stokes_preconditioner(
    Preconditioner* new_ns_preconditioner_pt)
  {
    // Check if pointer is non-zero.
    if (new_ns_preconditioner_pt == 0)
    {
      std::ostringstream warning_stream;
      warning_stream << "WARNING: \n"
                     << "The LSC preconditioner point is null.\n"
                     << "Using default (SuperLU) preconditioner.\n"
                     << std::endl;
      OomphLibWarning(
        warning_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);

      Navier_stokes_preconditioner_pt = 0;
      Using_superlu_ns_preconditioner = true;
    }
    else
    {
      // If the default SuperLU preconditioner has been used
      // clean it up now...
      if (Using_superlu_ns_preconditioner &&
          Navier_stokes_preconditioner_pt != 0)
      {
        delete Navier_stokes_preconditioner_pt;
        Navier_stokes_preconditioner_pt = 0;
      }

      Navier_stokes_preconditioner_pt = new_ns_preconditioner_pt;
      Using_superlu_ns_preconditioner = false;
    }
  } // func set_navier_stokes_preconditioner

  //========================================================================
  /// Clears the memory.
  //========================================================================
  void LagrangeEnforcedFlowPreconditioner::clean_up_memory()
  {
    // clean the block preconditioner base class memory
    this->clear_block_preconditioner_base();

    // Delete the Navier-Stokes preconditioner pointer if we have created it.
    if (Using_superlu_ns_preconditioner && Navier_stokes_preconditioner_pt != 0)
    {
      delete Navier_stokes_preconditioner_pt;
      Navier_stokes_preconditioner_pt = 0;
    }

    Preconditioner_has_been_setup = false;
  } // func LagrangeEnforcedFlowPreconditioner::clean_up_memory

} // namespace oomph
