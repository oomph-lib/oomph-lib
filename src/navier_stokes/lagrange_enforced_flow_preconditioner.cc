//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision$
//LIC//
//LIC// $LastChangedDate$
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
//LIC//=====================================================================
#include "lagrange_enforced_flow_preconditioner.h"

namespace oomph
{
//==========================================================================
/// Namespace for subsidiary preconditioner creation helper functions
//==========================================================================
namespace Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper
{
  /// \short CG with diagonal preconditioner for W-block subsidiary linear 
  /// systems.
  Preconditioner* get_w_cg_preconditioner()
  {
#ifdef OOMPH_HAS_TRILINOS
    InnerIterationPreconditioner
      <TrilinosAztecOOSolver,MatrixBasedDiagPreconditioner>* prec_pt = 
        new InnerIterationPreconditioner
          <TrilinosAztecOOSolver,MatrixBasedDiagPreconditioner>;

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
    throw OomphLibError(err_msg.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
#endif
  } // function get_w_cg_preconditioner
} // namespace

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

/// \short Apply the preconditioner.
/// r is the residual (rhs), z will contain the solution.
void LagrangeEnforcedFlowPreconditioner::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
{
#ifdef PARANOID
  if (Preconditioner_has_been_setup==false)
   {
    std::ostringstream error_message;
    error_message << "setup() must be called before using "
                  << "preconditioner_solve()";
    throw OomphLibError(
     error_message.str(),
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }
  if (z.built())
   {
    if (z.nrow() != r.nrow())
     {
      std::ostringstream error_message;
      error_message << "The vectors z and r must have the same number of "
                    << "of global rows";
      throw OomphLibError(
       error_message.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);      
     }
   }
#endif

  // if z is not setup then give it the same distribution
  if (!z.distribution_pt()->built())
   {
    z.build(r.distribution_pt(),0.0);
   }

  // Working vectors.
  DoubleVector temp_vec;
  DoubleVector another_temp_vec;
  DoubleVector yet_another_temp_vec;


  // -----------------------------------------------------------------------
  // Step 1 - apply approximate W block inverse to Lagrange multiplier
  // unknowns
  // -----------------------------------------------------------------------

  // If the W preconditioner is a block preconditioner, then the extraction
  // associated unknowns and RHS is handled by the (subsidiary) block 
  // preconditioner.
  if(W_preconditioner_is_block_preconditioner)
  {
    // Loop through the solve the subsystems associated with the W block.
    for (unsigned l_i = 0; l_i < N_lagrange_doftypes; l_i++) 
    {
      W_preconditioner_pt[l_i]->preconditioner_solve(r,z);
    }
  }
  else
  {
    // For each subsystem associated with each Lagrange multiplier, we loop
    // through and:
    // 1) extract the block entries from r
    // 2) apply the inverse
    // 3) return the entries to z.
    for(unsigned l_i = 0; l_i < N_lagrange_doftypes; l_i++)
    {
      // The Lagrange multiplier block type.
      const unsigned l_ii = N_fluid_doftypes + l_i;

      // Extract the block
      this->get_block_vector(l_ii,r,temp_vec);

      // Solve
      W_preconditioner_pt[l_i]
        ->preconditioner_solve(temp_vec,another_temp_vec);

//      // Apply the scaling sigma. The reason why we apply the scaling sigma
//      // here instead of within the W_preconditioner_pt is to reduce error 
//      // propagation. If sigma is small - then 1/sigma is large, we want 
//      // avoid large growth in values.
//      const unsigned vec_nrow_local = another_temp_vec.nrow_local();
//      double* vec_values_pt = another_temp_vec.values_pt();
//      for (unsigned i = 0; i < vec_nrow_local; i++) 
//      {
//        vec_values_pt[i] = vec_values_pt[i]*Scaling_sigma;
//      } // for

      // Return the unknowns
      this->return_block_vector(l_ii,another_temp_vec,z);

      // Clear vectors.
      temp_vec.clear();
      another_temp_vec.clear();
    } // for
  } // else

  // -----------------------------------------------------------------------
  // Step 2 - apply the augmented Navier-Stokes matrix inverse to the 
  // velocity and pressure unknowns
  // -----------------------------------------------------------------------

  // At this point, all vectors are cleared.
  if(Using_superlu_ns_preconditioner)
  {
    // Which block types corresponds to the fluid block types.
    Vector<unsigned> fluid_block_indices(N_fluid_doftypes,0);
    for (unsigned b = 0; b < N_fluid_doftypes; b++) 
    {
      fluid_block_indices[b] = b;
    }

    this->get_concatenated_block_vector(fluid_block_indices,r,temp_vec);

    // temp_vec contains the (concatenated) fluid rhs.
    Navier_stokes_preconditioner_pt
      ->preconditioner_solve(temp_vec,another_temp_vec);

    temp_vec.clear();

    // Now return it.
    this->return_concatenated_block_vector(fluid_block_indices,
        another_temp_vec,z);

    another_temp_vec.clear();
  }
  else
  {
    // This is a BlockPreconditioner
    Navier_stokes_preconditioner_pt->preconditioner_solve(r,z);
  }
} // end of preconditioner_solve

/// \short Set the meshes, the first mesh must be the bulk fluid mesh
void LagrangeEnforcedFlowPreconditioner::set_meshes(
    const Vector<Mesh*> &mesh_pt)
{
  // There should be at least two meshes for this preconditioner.
  const unsigned nmesh = mesh_pt.size();

#ifdef PARANOID
  if(nmesh < 2)
  {
    std::ostringstream err_msg;
    err_msg << "There should be at least two meshes.\n";
    throw OomphLibError(err_msg.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }

  // Check that all pointers are not null
  for(unsigned mesh_i = 0; mesh_i < nmesh; mesh_i++)
  {
    if (mesh_pt[mesh_i]==0)
    {
      std::ostringstream err_msg;
      err_msg << "The pointer mesh_pt[" << mesh_i << "] is null.\n";
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
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
    err_msg << "In the first mesh, the elements have elemental dimension of "
      << elemental_dim << ", with a nodal dimension of "
      << nodal_dim << ".\n"
      << "The first mesh does not contain 'bulk' elements.\n"
      << "Please re-order your mesh_pt vector.\n";

    throw OomphLibError(err_msg.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }
#endif

  // Set the number of meshes 
  this->set_nmesh(nmesh);

  // Set the meshes
  for(unsigned mesh_i = 0; mesh_i < nmesh; mesh_i++)
  {
    this->set_mesh(mesh_i,mesh_pt[mesh_i]);
  }

  // We also store the meshes and number of meshes locally in this class.
  // This is slightly redundant, since we always have it stored in the upper
  // most master block preconditioner. But at the moment there is no 
  // mapping/look up scheme between master and subsidiary block 
  // preconditioners for meshes.
  // So if this is a subsidiary block preconditioner, we don't know which of
  // the master's meshes belong to us. We need this information to set up
  // look up lists in the function setup(...).
  // So we store them local to this class.
  My_mesh_pt = mesh_pt;
  My_nmesh = nmesh;
} // EoFunc set_meshes


//===========================================================================
/// Setup the Lagrange enforced flow preconditioner. This
/// extracts blocks corresponding to the velocity and pressure unknowns,
/// creates the matrices actually needed in the application of the
/// preconditioner and deletes what can be deleted... Note that
/// this preconditioner needs a CRDoubleMatrix.
//============================================================================

//========================================================================
/// Setup method for the LagrangeEnforcedFlowPreconditioner.
//========================================================================
void LagrangeEnforcedFlowPreconditioner::setup()
{
  // clean
  this->clean_up_memory();

#ifdef PARANOID
  // Paranoid check that meshes have been set.
  if(My_nmesh == 0)
  {
    std::ostringstream err_msg;
    err_msg << "There are no meshes set. Please call set_meshes(...)\n"
      << "with at least two mesh.";
    throw OomphLibError(err_msg.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }
#endif


  // Set up the My_ndof_types_in_mesh vector.
  // - this can be optimised by creating a lookup list to the master
  // block preconditioner, so then we get the number of DOF types in each mesh
  // from BlockPreconditioner::ndof_types_in_mesh(...). The function
  //  Mesh::ndof_types() requires communication, so we store this for now.
  //
  //  We try to be a bit efficient by not setting up the list twice.
  if(My_ndof_types_in_mesh.size() == 0)
  {
    for (unsigned mesh_i = 0; mesh_i < My_nmesh; mesh_i++) 
    {
      My_ndof_types_in_mesh.push_back(My_mesh_pt[mesh_i]->ndof_types());
    }
  }


  // To construct the desired block structure for this preconditioner, with
  // assumption on the natural ordering of the DOF types, we require only the
  //
  // (1) spatial dimension of the problem and
  // (2) the number of DOF types in each of the meshes.
  //
  // Assumption: We assume the first mesh is always the "bulk" mesh 
  // (the Navier-Stokes mesh), which contains block preconditionable 
  // Navier-Stokes elements classifies the velocity and pressure 
  // degrees of freedom (ndof_types = 3(4) in 2(3)D). All subsequent meshes
  // contains the constrained velocity DOF types, then the Lagrange multiplier
  // DOF types.
  //
  // Example:
  // The general ordering for the DOF types (for this program, in 3D) is
  //
  //  0 1 2 3   4 5 6 7  8  ..x   x+0 x+1 x+2 x+3 x+4
  // [u v w p] [u v w l1 l2 ...] [u   v   w   l1  l2 ...] ...
  //
  // where then square brackets [] represents the DOF types in each mesh.
  //
  // Thus we need the mesh pointers to determine how many Lagrange multiplier
  // DOF types there are in each mesh, which is used to determine how many
  // Lagrange multiplier DOF types there are overall.
  //
  // Consider the case of imposing parallel outflow (3 constrained velocity
  // DOF types and 2 Lagrange multiplier DOF types) and tangential flow (3
  // constrained velocity DOF types and 1 Lagrange multiplier DOF type)
  // along two different boundaries in 3D. The resulting natural ordering of
  // the DOF types, determined by the by the two-level ordering, (first mesh 
  // order, then the elemental DOF type ordering) is:
  //
  // [0 1 2 3] [4  5  6   7   8 ] [9  10 11 12 ]
  // [u v w p] [up vp wp Lp1 Lp2] [ut vt wt Lt1]
  //
  // We observe that the problem dimension is 3 and 
  // N_doftype_in_mesh = [4, 5, 4].
  //
  // With these information we can construct the desired block structure:
  // | u v w | up vp wp | ut vt wt | p | Lp1 Lp2 Lt1 |
  //
  // The block structure is determined by the vector dof_to_block_map we 
  // give to the function block_setup(...).
  // 
  // First we work out a few variables to aid us in this endeavour.


  // Get the spatial dimension of the problem.
  // This is used to create look up lists later.
  // Note: Mesh::nodal_dimension() requires communication, use it sparingly!
  unsigned spatial_dim = My_mesh_pt[0]->nodal_dimension();

  // Get the number of DOF types.
  unsigned n_dof_types = ndof_types();

  // Check if the number of DOF types make sense.
#ifdef PARANOID
  // The Navier-Stokes mesh should (classify $spatial_dim + 1) number of 
  // DOF types (velocities and pressure).
  // For the meshes classifying constrained DOF types, there will always be
  // $spatial_dim number of constrained velocity DOF types first, followed by 
  // the Lagrange multiplier DOF types.
  //
  // Of course, we cannot check which DOF type is which, and thus cannot check 
  // if they are in the correct "assumed/required" order.
  //
  // We can at least check if the n_dof_types is the same as the total number
  // of DOF types in the meshes. This check is not necessary for the upper most
  // master block preconditioner since it gets the ndof_types() by looping 
  // through the meshes!
  if(is_subsidiary_block_preconditioner())
  {
    // Calculate the ndof types from the meshes passed to this block 
    // preconditioner. This should be equal to the value returned by 
    // ndof_types()
    unsigned tmp_ndof_types = 0;

    for (unsigned mesh_i = 0; mesh_i < My_nmesh; mesh_i++) 
    {
      tmp_ndof_types += My_ndof_types_in_mesh[mesh_i];
    }

    if(tmp_ndof_types != n_dof_types)
    {
      std::ostringstream err_msg;
      err_msg << "The number of DOF types are incorrect.\n"
        << "The total DOF types from the meshes is: " 
        << tmp_ndof_types << ".\n"
        << "The number of DOF types from "
        << "BlockPreconditioner::ndof_types() is " << n_dof_types <<".\n";
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif

  // The number of velocity DOF types: We assume that all meshes classify
  // at least some of the velocity DOF types (bulk/constrained), thus the 
  // total number of velocity DOF types is the spatial dimension multiplied
  // by the number of meshes.
  N_velocity_doftypes = My_nmesh*spatial_dim;
  //    std::cout << "N_velocity_doftypes: " << N_velocity_doftypes << std::endl; 

  // Fluid has +1 for the pressure.
  N_fluid_doftypes = N_velocity_doftypes + 1;
  //    std::cout << "N_fluid_doftypes: " << N_fluid_doftypes << std::endl; 

  // The rest are Lagrange multiplier DOF types.
  N_lagrange_doftypes = n_dof_types - N_fluid_doftypes;
  //    std::cout << "N_lagrange_doftypes: " << N_lagrange_doftypes << std::endl; 

  ///////////////////////////////////////////////////////////
  //   Now create the DOF to block map for block_setup()   //
  ///////////////////////////////////////////////////////////

  // This preconditioner permutes the DOF types to get the block types.
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
  // pressure / Lagrange multiplier
  // (in the case of the bulk mesh)/(in the case of the subsequent meshes).
  //
  // Consider the same example as shown previously, the natural DOF type 
  // ordering is:
  //  0 1 2 3   4  5  6   7   8    9  10 11 12   <- Natural DOF type ordering.
  // [u v w p] [up vp wp Lp1 Lp2] [ut vt wt Lt1] <- DOF type.
  //
  // We want:
  // [u v w | up vp wp | ut vt wt ] [p | Lp1 Lp2 Lt1]
  //
  // We can see that we lump the first $spatial_dimension number of DOF types
  // from each mesh together, and put the rest at the end.
  //
  // To achieve this, we use the dof_to_block_map:
  //  0 1 2 9 3  4  5  10  11  6  7  8  12   <- dof_to_block_map
  //
  // In general, we want:
  //
  //    dof_to_block_map[$dof_number] = $block_number
  //
  //
  // We put this all together for your convenience:
  //
  // Natural DOF type ordering:
  //  0 1 2 3   4  5  6   7   8    9  10 11 12   <- Natural DOF type ordering.
  // [u v w p] [up vp wp Lp1 Lp2] [ut vt wt Lt1] <- DOF type.
  //
  // Desired block structure:
  //  0 1 2   3  4  5    6  7  8     9   10  11  12   <- block index
  // [u v w | up vp wp | ut vt wt ] [p | Lp1 Lp2 Lt1] <- DOF type
  //
  // dof_to_block_map:
  //  0 1 2 9 3  4  5  10  11  6  7  8  12   <- dof_to_block_map
  // [u v w p up vp wp Lp1 Lp2 ut vt wt Lt1] <- DOF type
  //
  // Consider the 4th entry of the dof_to_block_map, this represents the 
  // pressure DOF type. We want to consider this block type 9.


  //  /////////////////////////////////////////////////////////////////////////////
  //  // Start of artificial test data
  //  n_dof_types = 13;
  //  My_nmesh = 3;
  //  spatial_dim = 3;
  //  N_velocity_doftypes = My_nmesh * spatial_dim;
  //  My_ndof_types_in_mesh.resize(3,0);
  //  My_ndof_types_in_mesh[0] = 4;
  //  My_ndof_types_in_mesh[1] = 5;
  //  My_ndof_types_in_mesh[2] = 4;
  //  // End of artificial test data.

  Vector<unsigned> dof_to_block_map(n_dof_types,0);

  // Encapsulate temporary variables
  {
    // Logic: We observe that the dof_to_block_map has the following pattern:
    // [0, 1, 2, p, 3, 4, 5, Lp1, ... Lpx, 6, 7, 8, Lt1, ... Lty, ...]
    //
    // We can loop through the number of meshes,
    //  fill in the first $spatial_dimension number of entries with velocity vals
    //  then fill in the rest with Lagrange DOF types.

    // Index for the dof_to_block_map vector.
    unsigned temp_index = 0;

    // Value for the velocity DOF type.
    unsigned velocity_val = 0;

    // Value for the pressure/Lagrange DOF type.
    unsigned pressure_lgr_val = N_velocity_doftypes;

    // Loop through the number of meshes.
    for (unsigned mesh_i = 0; mesh_i < My_nmesh; mesh_i++)
    {
      // Fill in the velocity DOF types.
      for (unsigned dim_i = 0; dim_i < spatial_dim; dim_i++) 
      {
        //dof_to_block_map[temp_index] = dim_i + mesh_i*spatial_dim;// fancy way.
        dof_to_block_map[temp_index++] = velocity_val++;
      } // for

      // Loop through the pressure/Lagrange multiplier DOF types.
      unsigned ndof_type_in_mesh_i = My_ndof_types_in_mesh[mesh_i];
      for (unsigned doftype_i = spatial_dim; 
          doftype_i < ndof_type_in_mesh_i; doftype_i++) 
      {
        dof_to_block_map[temp_index++] = pressure_lgr_val++;
      } // for
    } // for
  } // Encapsulation

  // With the artificial test data, we get the following list:
  // 0 1 2 9 3  4  5  10  11  6  7  8  12
  //    std::cout << "Lgr::setup dof_to_block_map: " << std::endl; 
  //    for (unsigned tmp_i = 0; tmp_i < dof_to_block_map.size(); tmp_i++) 
  //    {
  //      std::cout << dof_to_block_map[tmp_i] << " "; 
  //    }
  //    std::cout << "\n" << std::endl; 
  //    pause("I'm back"); 


  // Call the block setup
  this->block_setup(dof_to_block_map);

  //    pause("Lgr::setup() done block_setup, about to print dof block dist nrow"); 

  //    for (unsigned block_i = 0; block_i < nblock_types(); block_i++) 
  //    {
  //      const unsigned nrow = dof_block_distribution_pt(block_i)->nrow();
  //      std::cout << "dof block " << block_i << ", nrow = " << nrow << std::endl; 
  //    }

  //    pause("done dof block dist nrow, about to print out each block.");


  //    std::cout << "Print 1, NOT replaced." << std::endl; 
  //    std::cout << "============================================" << std::endl; 
  //    std::cout << "============================================" << std::endl; 
  //    std::cout << "============================================" << std::endl; 
  //    std::cout << "============================================" << std::endl; 


  //  pause("Dumped!"); 


  //////////////////////////////////////////////////////////////////////////
  // Need to create the norms, used for Sigma, if required
  //////////////////////////////////////////////////////////////////////////

  // To do: This could be made more efficient by storing only the
  // augmented blocks. I will do this later.
  //
  // Extract the velocity block. Although we only require the infinity norm
  // a single direction of the velocity block, we extract them all since we
  // use the rest immediately (to perform the augmentation).
  DenseMatrix<CRDoubleMatrix*> v_aug_pt(N_velocity_doftypes,
      N_velocity_doftypes,0);

  for(unsigned row_i = 0; row_i < N_velocity_doftypes; row_i++)
  {
    for(unsigned col_i = 0; col_i < N_velocity_doftypes; col_i++)
    {
      v_aug_pt(row_i,col_i) = new CRDoubleMatrix;
      this->get_block(row_i,col_i,*v_aug_pt(row_i,col_i));
    } // for
  } // for



  //  unsigned tempnblocks = this->nblock_types();
  //  // RRR Checking the block dimensions.
  //  for (unsigned col_i = 0; col_i < tempnblocks; col_i++) 
  //  {
  //    CRDoubleMatrix tempblock;
  //    this->get_dof_level_block(0,col_i,tempblock);
  //    unsigned tempncol = tempblock.ncol();
  //    std::cout << "ncol: " << tempncol << std::endl; 
  //  }
  //std::cout << std::endl; 
  //
  //
  //  std::cout << "block ordering:" << std::endl; 
  //  
  //  // RRR Checking the block dimensions.
  //  for (unsigned col_i = 0; col_i < tempnblocks; col_i++) 
  //  {
  //    CRDoubleMatrix tempblock = this->get_block(0,col_i);
  //    unsigned tempncol = tempblock.ncol();
  //    std::cout << "ncol: " << tempncol << std::endl; 
  //    
  //  }
  //  pause("done!"); 


  if(Use_norm_f_for_scaling_sigma)
  {
    //    DenseMatrix<CRDoubleMatrix* > u_pt(My_nmesh,My_nmesh,0);
    //
    //    // Recall the blocking scheme:
    //    // 0 1 2 3  4  5  6  7  8
    //    // u v w u1 v1 w1 u2 v2 w2 ...
    //    // So we loop through the meshes and then get the first entry...
    //    // from v_aug_pt.
    //    for(unsigned row_i = 0; row_i < My_nmesh; row_i++)
    //     {
    //      for(unsigned col_i = 0; col_i < My_nmesh; col_i++)
    //       {
    //        u_pt(row_i,col_i) = v_aug_pt(row_i*spatial_dim,
    //                                     col_i*spatial_dim);
    //       } // for
    //     } //  for

    // Get the norm and set the scaling sigma.
    //    Scaling_sigma = -CRDoubleMatrixHelpers::inf_norm(u_pt)
    //     *Scaling_sigma_multiplier;

    Scaling_sigma = -CRDoubleMatrixHelpers::inf_norm(v_aug_pt);

    //    Scaling_sigma = CRDoubleMatrixHelpers::maximum_gershgorin_disk(v_aug_pt)
    //     *Scaling_sigma_multiplier;
  } // if(Use_default_f_scaling)

#ifdef PARANOID
  if(Scaling_sigma == 0.0)
  {
    std::ostringstream warning_stream;
    warning_stream << "WARNING: " << std::endl
      << "The scaling (Scaling_sigma) is " 
      << Scaling_sigma << std::endl
      << "Division by Scaling_sigma = 0 will implode the world."
      << std::endl;
    OomphLibWarning(warning_stream.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }
  if(Scaling_sigma > 0.0)
  {
    std::ostringstream warning_stream;
    warning_stream << "WARNING: " << std::endl
      << "The scaling (Scaling_sigma) is positive: " 
      << Scaling_sigma << std::endl
      << "Performance may be degraded."
      << std::endl;
    OomphLibWarning(warning_stream.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }
#endif


  //////////////////////////////////////////////////////////////////////////
  //                Now create the augmented fluid matrix.                //
  //////////////////////////////////////////////////////////////////////////


  // Loop through the Lagrange multipliers and do three things:
  // For each Lagrange block:
  //   Create an empty w_i, this will be constructed incrementally from the
  //   mass matrices (in the columns).
  //
  //   For each block column
  //     1) Store the location of mass matrix.
  //     2) Store the mass matrix.
  //     3) Create the partial w_i += diag(m)^2
  //
  //   1) Store the w_i for this Lagrange multiplier.
  //   Loop through the mass matrix location to determine the aug. row.
  //     Loop through the mm location to determine the aug. col.
  //       1) Create inv_wi from the w_i already created.
  //       1) Create augmentation m * inv(w_i) * m
  //       2) Add the augmentation as determined by aug. row and aug. col.

  // Storage for the W block.
  Vector<CRDoubleMatrix*> w_pt(N_lagrange_doftypes,0);

  // Note that we do not need to store all the inverse w_i since they
  // are only used once per Lagrange multiplier.
  for(unsigned l_i = 0; l_i < N_lagrange_doftypes; l_i++)
  {
    // Storage for the current Lagrange block mass matrices.
    Vector<CRDoubleMatrix*> mm_pt;
    Vector<CRDoubleMatrix*> mmt_pt;

    // Get the current Lagrange DOF type.
    // Recall that the ordering is:
    // | Fluid DOF types | Lagrange multiplier DOF types
    //   u v w uc vc wc p  L1, ....
    unsigned l_doftype = N_fluid_doftypes + l_i;

    // Store the mass matrix locations for the current Lagrange block.
    Vector<unsigned> mm_locations;

    // Store the number of mass matrices.
    unsigned n_mm = 0;

    // It is easier to create the diagonal matrix w_i after
    // I have extracted all the mass matrices for this Lagrange
    // block.

    // Go along the block columns for the current Lagrange block ROW.
    for(unsigned col_i = 0; col_i < N_velocity_doftypes; col_i++)
    {
      // Get the block matrix for this block column.
      CRDoubleMatrix* mm_temp_pt = new CRDoubleMatrix;
      this->get_block(l_doftype, col_i, *mm_temp_pt);

      if(mm_temp_pt->nnz() > 0)
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
        << "Perhaps the problem setup is incorrect." << std::endl;
      OomphLibWarning(warning_stream.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Get the transpose of the mass matrices.
    for(unsigned mm_i = 0; mm_i < n_mm; mm_i++)
    {
      // Get the block matrix for this block column.
      CRDoubleMatrix* mm_temp_pt = new CRDoubleMatrix;
      this->get_block(mm_locations[mm_i],l_doftype,*mm_temp_pt);

      if(mm_temp_pt->nnz() > 0)
      {
        mmt_pt.push_back(mm_temp_pt);
      }
      else
      {
        // There should be a non-zero mass matrix here, since we are getting
        // mass matrices from the L^T block, corresponding to L.
#ifdef PARANOID
        {
          std::ostringstream warning_stream;
          warning_stream << "WARNING:\n"
            << "The mass matrix block " << mm_locations[mm_i] 
            << " in L^T block " << l_i << " is zero.\n"
            << "Perhaps the problem setup is incorrect." << std::endl;
          OomphLibWarning(warning_stream.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }
    } // loop through the ROW of the Lagrange COLUMN.

    // The mass matrices for the current Lagrange block row is in mm_pt. 
    // The mass matrices for the current Lagrange block col is in mmt_pt. 
    // Now create the w_i and put it in w_pt.

    // Create both the w_i and inv_w_i matrices.
    w_pt[l_i] = new CRDoubleMatrix(this->Block_distribution_pt[l_doftype]);
    CRDoubleMatrix* inv_w_pt 
      = new CRDoubleMatrix(this->Block_distribution_pt[l_doftype]);

    // The w_i is either a
    // 1) diagonal or
    // 2) block diagonal
    // approximation of the actual W = LL^T block 
    // (LL^T has off-diagonal blocks).
    //
    // When w_i is inverted, we always take the diagonal approximation 
    // beforehand. In the case of (1), nothing needs to be done, 
    // but in (2) we need to extract the diagonal of the 
    // block-diagonal matrix.

    {
      // Get the number of local rows for this Lagrange block.
      unsigned long l_i_nrow_local 
        = this->Block_distribution_pt[l_doftype]->nrow_local();

      // The first row, for the column offset (required in parallel).
      unsigned l_i_first_row 
        = this->Block_distribution_pt[l_doftype]->first_row();

      // A vector to contain the results of mass matrices squared.
      Vector<double> w_i_diag_values(l_i_nrow_local,0);
      Vector<double> invw_i_diag_values(l_i_nrow_local,0);

      Vector<int> w_i_column_indices(l_i_nrow_local);
      Vector<int> w_i_row_start(l_i_nrow_local+1);

      // Component-wise, square and add all the diagonals.
      // Loop through the mass matrices.
      for(unsigned m_i = 0; m_i < n_mm; m_i++)
      {
        // Get the diagonal entries for the mass matrix.
        CRDoubleMatrix* temp_mm_sqrd_pt = new CRDoubleMatrix;
        temp_mm_sqrd_pt->build(mm_pt[m_i]->distribution_pt());
        mm_pt[m_i]->multiply(*mmt_pt[m_i],*temp_mm_sqrd_pt);

        Vector<double> m_diag = temp_mm_sqrd_pt->diagonal_entries();
        //        Vector<double> m_diag = mm_pt[m_i]->diagonal_entries();

        // Loop through the entries, add them.
        for(unsigned long row_i = 0; row_i < l_i_nrow_local; row_i++)
        {
          w_i_diag_values[row_i] += m_diag[row_i];
          // w_i_diag_values[row_i] += (m_diag[row_i]* m_diag[row_i]);
        }

        delete temp_mm_sqrd_pt; temp_mm_sqrd_pt = 0;
      }


      // NOTE: We do scale with sigma here!
      // Divide by Scaling_sigma and create the inverse of w.
      for(unsigned long row_i = 0; row_i < l_i_nrow_local; row_i++)
      {

        // w_i is a diagonal matrix, so take the inverse to
        // invert the matrix.
        invw_i_diag_values[row_i] = Scaling_sigma / w_i_diag_values[row_i];

        w_i_diag_values[row_i] /= Scaling_sigma;
        w_i_column_indices[row_i] = row_i + l_i_first_row;
        w_i_row_start[row_i] = row_i;
      }

      w_i_row_start[l_i_nrow_local] = l_i_nrow_local;

      // Theses are square matrices. So we use the l_i_nrow_global for the
      // number of columns.
      unsigned long l_i_nrow_global 
        = this->Block_distribution_pt[l_doftype]->nrow();
      w_pt[l_i]->build(l_i_nrow_global,
          w_i_diag_values,
          w_i_column_indices,
          w_i_row_start);
      inv_w_pt->build(l_i_nrow_global,
          invw_i_diag_values,
          w_i_column_indices,
          w_i_row_start);
    }


    ////////////////////////////////////////////////////////////////////////
    // Now we create the augmented matrix in v_aug_pt.
    // v_aug_pt is already re-ordered
    // Loop through the mm_locations
    for(unsigned ii = 0; ii < n_mm; ii++)
    {
      unsigned aug_i = mm_locations[ii];

      for(unsigned jj = 0; jj < n_mm; jj++)
      {
        unsigned aug_j = mm_locations[jj];

        // A temp pointer to store the intermediate results.
        CRDoubleMatrix* aug_pt = new CRDoubleMatrix;

        mmt_pt[ii]->multiply((*inv_w_pt),(*aug_pt));
        CRDoubleMatrix* another_mat_pt = new CRDoubleMatrix;

        aug_pt->multiply(*mm_pt[jj],(*another_mat_pt));

        //add_matrices(aug_pt,v_aug_pt(aug_i,aug_j));
        v_aug_pt(aug_i,aug_j)->add(*another_mat_pt,*v_aug_pt(aug_i,aug_j));

        delete aug_pt; aug_pt = 0;
        delete another_mat_pt; another_mat_pt = 0;
      } // loop jj
    } // loop ii

    delete inv_w_pt;
    for(unsigned m_i = 0; m_i < n_mm; m_i++)
    {
      delete mm_pt[m_i];
      delete mmt_pt[m_i];
    }
  } // loop through Lagrange multipliers.

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////


  //      std::cout << "\n" << std::endl; 
  //      std::cout << "\n" << std::endl; 
  //      std::cout << "Print 2, Setting replacement blocks." << std::endl; 
  //      std::cout << "============================================" << std::endl; 
  //      std::cout << "============================================" << std::endl; 
  //      std::cout << "============================================" << std::endl; 
  //      std::cout << "============================================" << std::endl; 

  // Recall that the DOF type ordering is:
  // The general ordering for the DOF types (for this program, in 3D) is
  //
  //  0 1 2 3   4 5 6 7  8  ..x   x+0 x+1 x+2 x+3 x+4
  // [u v w p] [u v w l1 l2 ...] [u   v   w   l1  l2 ...] ...
  //
  // Eg:
  // [0 1 2 3] [4  5  6   7   8 ] [9  10 11 12 ]
  // [u v w p] [up vp wp Lp1 Lp2] [ut vt wt Lt1]
  //
  // With these information we can construct the desired block structure:
  //   0 1 2   3  4  5    6  7  8    9   10  11  12
  // | u v w | up vp wp | ut vt wt | p | Lp1 Lp2 Lt1 |
  //
  // So we need to map:
  // 3 -> 4
  // 4 -> 5
  // 5 -> 6
  //
  // 6 -> 9
  // 7 -> 10
  // 8 -> 11
  //
  // Actually, it will be easier to go from:
  //
  // 1 -> 4
  // 2 -> 5
  // 3 -> 6
  //
  // 4 -> 9
  // 5 -> 10
  // 6 -> 11
  //
  // Since we just get the number of DOF types in the meshes, and we know
  // the number of meshes.

  //    oomph_info << "My_nmesh: " << My_nmesh << std::endl;
  //    for (unsigned mesh_i = 0; mesh_i < My_nmesh; mesh_i++) 
  //    {
  //      oomph_info << "Doftypes in mesh " 
  //                 << mesh_i << ": " 
  //                 << My_ndof_types_in_mesh[mesh_i] << std::endl; 
  //      
  //    }

  /// Replace all f blocks
  //  if(Replace_all_f_blocks)
  {
    // Create the inverse to dof_to_block_map
    Vector<unsigned> temp_block_to_dof_map(n_dof_types,0);
    for (unsigned dof_i = 0; dof_i < n_dof_types; dof_i++) 
    {
      temp_block_to_dof_map[dof_to_block_map[dof_i]] = dof_i;
    }

    // print it out to check.

    //    std::cout << "dof_to_block_map:" << std::endl;
    //    for (unsigned iii = 0; iii < dof_to_block_map.size(); iii++) 
    //    {
    //      std::cout << dof_to_block_map[iii] << std::endl; 
    //    }
    //
    //    std::cout << std::endl;
    //    std::cout << "block_to_dof_map: " << std::endl;
    //    for (unsigned iii = 0; iii < block_to_dof_map.size(); iii++) 
    //    {
    //      std::cout << block_to_dof_map[iii] << std::endl; 
    //    }


    // Now do the replacement of all blocks in v_aug_pt 
    for (unsigned block_row_i = 0; 
        block_row_i < N_velocity_doftypes; block_row_i++) 
    {
      unsigned temp_dof_row_i = temp_block_to_dof_map[block_row_i];
      for (unsigned block_col_i = 0; 
          block_col_i < N_velocity_doftypes; block_col_i++) 
      {
        unsigned temp_dof_col_i = temp_block_to_dof_map[block_col_i];
        this->set_replacement_dof_block(
            temp_dof_row_i,temp_dof_col_i,
            v_aug_pt(block_row_i,block_col_i));
      }
    }
  }
  //  else
  //  {
  //    unsigned blocked_row_i = spatial_dim;
  //    unsigned blocked_col_i = spatial_dim;
  //    unsigned doftype_row_running_total = My_ndof_types_in_mesh[0];
  //    unsigned doftype_col_running_total = My_ndof_types_in_mesh[0];
  //
  //    for (unsigned row_mesh_i = 1; row_mesh_i < My_nmesh; row_mesh_i++) 
  //    {
  //      // Get the indirection for the mesh
  //      for (unsigned row_dim_i = 0; row_dim_i < spatial_dim; row_dim_i++) 
  //      {
  //        const unsigned doftype_row_i = doftype_row_running_total 
  //          + row_dim_i;
  //
  //        // Now do the same for the columns
  //        for (unsigned col_mesh_i = 1; col_mesh_i < My_nmesh; col_mesh_i++) 
  //        {
  //          for (unsigned col_dim_i = 0; col_dim_i < spatial_dim; col_dim_i++)
  //          {
  //            const unsigned doftype_col_i = doftype_col_running_total
  //              + col_dim_i;
  //
  //            this->set_replacement_dof_block(
  //                doftype_row_i,doftype_col_i,
  //                v_aug_pt(blocked_row_i,blocked_col_i));
  //
  //            //            oomph_info << "(" << blocked_row_i << "," 
  //            //                              << blocked_col_i << ")" 
  //            //                       << " -> " 
  //            //                       << "(" << doftype_row_i << ","
  //            //                              << doftype_col_i << ")" << std::endl; 
  //
  //            blocked_col_i++;
  //          }
  //
  //          // Update the column dof type running total
  //          doftype_col_running_total += My_ndof_types_in_mesh[col_mesh_i];
  //        }
  //
  //        // reset the blocked col i and the running total for the next row
  //        blocked_col_i = spatial_dim;
  //        doftype_col_running_total = My_ndof_types_in_mesh[0];
  //
  //        // Update the blocked row_i
  //        blocked_row_i++;
  //      }
  //      // Update the row dof type running total
  //      doftype_row_running_total += My_ndof_types_in_mesh[row_mesh_i];
  //    }
  //  }


  // AT this point, we have created the augmented fluid block in v_aug_pt
  // and the w block in w_pt.
  //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  // Setup the fluid subsidiary preconditioner
  //
  // We solve the fluid block using different preconditioners.
  //
  // 1) For exact block preconditioning, we use the SuperLU solver.
  // For this, we have to form the fluid block by extracting the 
  // pressure matrices.
  //
  // 2) For Exact LSC block preconditioning we do not need to form the whole
  // fluid block since the pressure and velocity are solved separately.
  //
  if(Using_superlu_ns_preconditioner)
  {

    //      pause("Lgr::setup() done setting replacement blocks, about to check repl. blocks");
    //      MapMatrix<unsigned,CRDoubleMatrix*> re_block_pt = this->replacement_dof_block_pt();
    //
    //      for (unsigned dof_block_i = 0; dof_block_i < this->ndof_types(); dof_block_i++) 
    //      {
    //        for (unsigned dof_block_j = 0; dof_block_j < this->ndof_types(); dof_block_j++) 
    //        {
    //          CRDoubleMatrix* tmp_block_pt = re_block_pt.get(dof_block_i,dof_block_j);
    //          if(tmp_block_pt != 0)
    //          {
    //            const unsigned tmp_nrow = tmp_block_pt->nrow();
    //            const unsigned tmp_ncol = tmp_block_pt->ncol();
    //            std::cout << "block(" << dof_block_i <<","<<dof_block_j
    //                      <<") has been replaced, nrow: " 
    //                      << tmp_nrow << ", ncol: "<< tmp_ncol << std::endl;
    //          }
    //        }
    //      }
    //      pause("printed repl. block pt"); 
    //      
    //      
    //
    //
    //
    //
    //      std::cout << "Print 1, NOT replaced." << std::endl; 
    //      std::cout << "============================================" << std::endl; 
    //      std::cout << "============================================" << std::endl; 
    //      std::cout << "============================================" << std::endl; 
    //      std::cout << "============================================" << std::endl; 
    //      
    //      for (unsigned block_i = 0; block_i < nblock_types(); block_i++) 
    //      {
    //        for (unsigned block_j = 0; block_j < nblock_types(); block_j++) 
    //        {
    //          CRDoubleMatrix tmp_block = get_block(block_i,block_j);
    //          const unsigned tmp_block_nrow = tmp_block.nrow();
    //          const unsigned tmp_block_ncol = tmp_block.ncol();
    //          std::cout << "block(" << block_i  << "," << block_j << ")"
    //                    << " , nrow: " << tmp_block_nrow << ", col: " << tmp_block_ncol << std::endl; 
    //
    //          std::cout << "============================================" << std::endl; 
    //          std::cout << "\n" << std::endl; 
    //        }
    //      }
    //      
    //      pause("Got the nrow yo 222222222"); 
    //
    //      DenseMatrix<CRDoubleMatrix* > f_subblock_pt(N_fluid_doftypes,
    //                                                  N_fluid_doftypes,0);
    //      // put in v_aug_pt:
    //      for(unsigned v_i = 0; v_i < N_fluid_doftypes; v_i++)
    //      {
    //        for(unsigned v_j = 0; v_j < N_fluid_doftypes; v_j++)
    //        {
    //          f_subblock_pt(v_i,v_j) = new CRDoubleMatrix;
    //          this->get_block(v_i,v_j,*f_subblock_pt(v_i,v_j));
    //        }
    //      }


    //      DenseMatrix<CRDoubleMatrix* > f_subblock_pt(N_fluid_doftypes,
    //          N_fluid_doftypes,0);
    //      // put in v_aug_pt:
    //      for(unsigned v_i = 0; v_i < N_velocity_doftypes; v_i++)
    //      {
    //        for(unsigned v_j = 0; v_j < N_velocity_doftypes; v_j++)
    //        {
    //          f_subblock_pt(v_i,v_j) = v_aug_pt(v_i,v_j);
    //        }
    //      }
    //
    //      // Now get the B and B^T blocks. Note that we have left out the zero block.
    //
    //      // Fill in the pressure block B
    //      for(unsigned col_i = 0; col_i < N_velocity_doftypes; col_i++)
    //      {
    //        f_subblock_pt(N_velocity_doftypes,col_i) = new CRDoubleMatrix;
    //        this->get_block(N_velocity_doftypes,col_i,
    //            *f_subblock_pt(N_velocity_doftypes,col_i));
    //      }
    //
    //      // Fill in the pressure block B^T
    //      for(unsigned row_i = 0; row_i < N_velocity_doftypes; row_i++)
    //      {
    //        f_subblock_pt(row_i,N_velocity_doftypes) = new CRDoubleMatrix;
    //        this->get_block(row_i,N_velocity_doftypes,
    //            *f_subblock_pt(row_i,N_velocity_doftypes));
    //      }

    // Concatenate the sub matrices.
    //      CRDoubleMatrix* f_aug_pt = new CRDoubleMatrix;


    //      std::cout << "dof_block_dimension():" << std::endl;
    //      for (unsigned i = 0; i < 6; i++) 
    //      {
    //        std::cout << "i = " << i << ", val: " << this->dof_block_dimension(i) << std::endl;
    //        
    //      }
    //
    //      pause("bloodz"); 


    VectorMatrix<BlockSelector> f_aug_blocks(N_fluid_doftypes,
        N_fluid_doftypes);
    for (unsigned block_i = 0; block_i < N_fluid_doftypes; block_i++) 
    {
      for (unsigned block_j = 0; block_j < N_fluid_doftypes; block_j++) 
      {
        f_aug_blocks[block_i][block_j].select_block(block_i,block_j,true);
      }
    }

    ///////////////////////////////////////////////////////////////////////////////      
    ///////////////////////////////////////////////////////////////////////////////      
    ///////////////////////////////////////////////////////////////////////////////      

    CRDoubleMatrix f_aug_block 
      = this->get_concatenated_block(f_aug_blocks);

    //      pause("after get_concat_blocks"); 

    //      Vector<LinearAlgebraDistribution*> f_dist_pt(N_fluid_doftypes,0);
    //      for(unsigned f_i = 0; f_i < N_fluid_doftypes; f_i++)
    //      {
    //        f_dist_pt[f_i] = this->block_distribution_pt(f_i);
    //      }

    //      CRDoubleMatrixHelpers::concatenate_without_communication
    //        (f_dist_pt,f_subblock_pt,*f_aug_pt);

    // delete the sub F pointers, this would also delete the 
    // v_aug_pt
    //      for(unsigned row_i = 0; row_i < N_fluid_doftypes; row_i++)
    //      {
    //        for(unsigned col_i = 0; col_i < N_fluid_doftypes; col_i++)
    //        {
    //          delete f_subblock_pt(row_i, col_i);
    //          f_subblock_pt(row_i, col_i) = 0;
    //        }
    //      }

    if(Navier_stokes_preconditioner_pt == 0)
    {
      Navier_stokes_preconditioner_pt = new SuperLUPreconditioner;
    }

    Navier_stokes_preconditioner_pt->setup(&f_aug_block);

    //      Navier_stokes_preconditioner_pt->setup(f_aug_pt);

    //      delete f_aug_pt;
    //      f_aug_pt = 0;
  }
  else
  {
    //      std::cout << "\n" << std::endl; 
    //      std::cout << "\n" << std::endl; 
    //      std::cout << "\n" << std::endl; 

    //      pause("Lgr::setup() Got to the else block"); 

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    //      std::cout << "Lgr using a block subsidiary preconditioner" << std::endl; 


    //      // Get the rest of the f block.
    //      DenseMatrix<CRDoubleMatrix* > f_subblock_pt(N_fluid_doftypes,
    //          N_fluid_doftypes,0);
    //      // put in v_aug_pt:
    //      for(unsigned v_i = 0; v_i < N_velocity_doftypes; v_i++)
    //      {
    //        for(unsigned v_j = 0; v_j < N_velocity_doftypes; v_j++)
    //        {
    //          f_subblock_pt(v_i,v_j) = v_aug_pt(v_i,v_j);
    //        }
    //      }

    //      // Fill in the pressure block B plus 1
    //      for(unsigned col_i = 0; col_i < N_fluid_doftypes; col_i++)
    //      {
    //        f_subblock_pt(N_velocity_doftypes,col_i) = new CRDoubleMatrix;
    //
    //        this->get_block(N_velocity_doftypes,col_i,
    //            *f_subblock_pt(N_velocity_doftypes,col_i));
    //      }

    //      // Fill in the pressure block B^T
    //      for(unsigned row_i = 0; row_i < N_velocity_doftypes; row_i++)
    //      {
    //        f_subblock_pt(row_i,N_velocity_doftypes) = new CRDoubleMatrix;
    //
    //        this->get_block(row_i,N_velocity_doftypes,
    //            *f_subblock_pt(row_i,N_velocity_doftypes));
    //      }


    // Determine whether the NS preconditioner is a block preconditioner (and
    // therefore a subsidiary preconditioner)
#ifdef PARANOID
    BlockPreconditioner<CRDoubleMatrix>* navier_stokes_block_preconditioner_pt
      = dynamic_cast<BlockPreconditioner<CRDoubleMatrix>* >
      (Navier_stokes_preconditioner_pt);

    if(navier_stokes_block_preconditioner_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "Navier-stokes preconditioner is not a block\n"
        << "preconditioner." << std::endl;
      throw OomphLibError(
          error_message.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
#else
    BlockPreconditioner<CRDoubleMatrix>* navier_stokes_block_preconditioner_pt
      = static_cast<BlockPreconditioner<CRDoubleMatrix>* >
      (Navier_stokes_preconditioner_pt);
#endif


    // Tell the LSC preconditioner which DOF type should be treated as one
    // These DOF types are local to the LSC preconditioner.
    // For example, for the artificial test data:
    // 
    // ub vb wb p up vp wp Lp1 Lp2 ut vt wt Lt
    // 0  1  2  3 4  5  6  7   8   9  10 11  12 DOF type order.
    //
    // Assume that we gave this list to turn_into_subsidiary:
    // 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 MARK1 (see below)
    // ub vb wb up vp wp ut vt wt p
    // 0  1  2  4  5  6  9  10 11 3 turn_into_sub list.
    //
    // What the turn_into_sub list says is: your dof type
    // 0 is my dof type 0
    // 1 is my dof type 1
    // 2 is my dof type 2
    // 3 is my dof type 4 (we've skipped a pressure).
    // etc....
    // your dof type 9 is my dof type 3 (for the pressure)
    //
    // NOW: We have to tell the LSC preconditioner which dof types to
    // treat as one. The LSC preconditioner 
    // expects 2 or 3 velocity dof types, and one pressure DOF types.
    // We give it this list (Look at the list MARK1):
    // u [0, 3, 6]
    // v [1, 4, 7]
    // w [2, 5, 8]
    // p [9]
    //
    // Artificial test data:
    //    spatial_dim = 3;
    //    My_nmesh = 3;

    Vector<Vector<unsigned> > subsidiary_dof_type_coarsening_map;

    for (unsigned direction = 0; direction < spatial_dim; direction++)
    {
      Vector<unsigned> dir_doftypes_vec(My_nmesh,0);
      for (unsigned mesh_i = 0; mesh_i < My_nmesh; mesh_i++) 
      {
        dir_doftypes_vec[mesh_i] = spatial_dim*mesh_i+direction;
      }

      // Push it in!
      subsidiary_dof_type_coarsening_map.push_back(dir_doftypes_vec);
    }

    Vector<unsigned> ns_p_vec(1,0);
    // This is simply the number of velocity dof types,
    // But we work it out anyway so the artificial test data is used.
    ns_p_vec[0] = My_nmesh*spatial_dim; 

    subsidiary_dof_type_coarsening_map.push_back(ns_p_vec);

    //    std::cout << "Lgr: sub dof type coarsening map: " << std::endl; 
    //    for (unsigned i = 0; 
    //        i < subsidiary_dof_type_coarsening_map.size(); i++) 
    //    {
    //      for (unsigned j = 0; 
    //          j < subsidiary_dof_type_coarsening_map[i].size(); j++) 
    //      {
    //        std::cout << subsidiary_dof_type_coarsening_map[i][j] << " ";
    //      }
    //      std::cout << "\n"; 
    //    }
    //    pause("Outputted the dof type coarsen map."); 



    //    Doftype_coarsen_map_fine.resize(0,0);
    //    Vector<unsigned> tmp0;
    //    tmp0.push_back(0);
    //    tmp0.push_back(1);
    //    tmp0.push_back(2);
    //    tmp0.push_back(3);
    //    Vector<unsigned> tmp1;
    //    tmp1.push_back(4);
    //    tmp1.push_back(5);
    //    tmp1.push_back(6);
    //    tmp1.push_back(7);
    //    Vector<unsigned> tmp2;
    //    tmp2.push_back(8);
    //    tmp2.push_back(9);
    //    tmp2.push_back(10);
    //    tmp2.push_back(11);
    //    Vector<unsigned> tmp3;
    //    tmp3.push_back(12);
    //    tmp3.push_back(13);
    //    Vector<unsigned> tmp4;
    //    tmp4.push_back(14);
    //    tmp4.push_back(15);
    //    Doftype_coarsen_map_fine.push_back(tmp0);
    //    Doftype_coarsen_map_fine.push_back(tmp1);
    //    Doftype_coarsen_map_fine.push_back(tmp2);
    //    Doftype_coarsen_map_fine.push_back(tmp3);
    //    Doftype_coarsen_map_fine.push_back(tmp4);

    //    pause("Before turn into..."); 


    // Now we create the dof_number_in_master_map vector to pass to the
    // turn_into_subsidiary_block_preconditioner(...) function.
    //
    // This vector tells the subsidiary block preconditioner what it's DOF type
    // is in relation it's parent block preconditioner.
    //
    // In general, we want:
    //
    //   dof_number_in_master_map[subsidiary DOF type] = master DOF type
    //
    // 
    // Example: Using the example above, our problem has the natural 
    // DOF type ordering:
    //
    //  0 1 2 3   4  5  6   7   8    9  10 11 12   <- Natural DOF type ordering.
    // [u v w p] [up vp wp Lp1 Lp2] [ut vt wt Lt1] <- DOF type.
    //
    // We need to create a vector telling the subsidiary block preconditioner
    // which DOF types it will work with, in our case, it will work with the
    // fluid DOF types: u, v, w, p, up, vp, wp, ut, vt and wt.
    //
    // Furthermore, the ordering of the DOF type matters, most importantly,
    // THIS preconditioner needs to be aware of the natural DOF type required
    // by the subsidiary block preconditioner. The Navier-Stokes Schur Complement
    // block preconditioner works with $spatial_dimension + 1 DOF types, BUT, 
    // this block preconditioner has split the velocity DOF types into a finer 
    // grain DOF type structure required by the Navier-Stokes preconditioner. 
    // We will discuss coarsening DOF types later, for now we focus on the
    // dof_number_in_master_map vector.
    //
    // We start by discussing a simpler example, lets say this preconditioner 
    // does not have a finer DOF type splitting than required by the subsidiary
    // block preconditioner, so let the natural DOF type ordering of this
    // block preconditioner be:
    // 0 1 2 3 4  5
    // u v w p L1 L2
    //
    // This preconditioner must know the natural DOF type ordering required by
    // the subsidiary block preconditioner. For our implementation of the LSC
    // preconditioner, we require the velocity DOF types first, then the pressure
    // DOF type. In 3D, this would be:
    //
    // LSC block preconditioner DOF type ordering:
    // 0 1 2 3
    // u v w p
    //
    // We must map each of the subsidiary DOF type with a parent DOF type.
    //
    // So the dof_number_in_master_map would be:
    // [0, 1, 2, 3], which says:
    //
    // subsidiary DOF type 0 -> master DOF type 0
    // subsidiary DOF type 1 -> master DOF type 1
    // subsidiary DOF type 2 -> master DOF type 2
    // subsidiary DOF type 3 -> master DOF type 3
    //
    // This is a trivial mapping. 
    //
    //
    //
    // Another example, lets assume that the natural DOF type ordering of this 
    // preconditioner be:
    // 0  1  2 3 4 5
    // L1 L2 u v w p
    //
    // Then the dof_number_in_master_map would be:
    // [2, 3, 4, 5], which says:
    // 
    // subsidiary DOF type 0 -> master DOF type 2
    // subsidiary DOF type 1 -> master DOF type 3
    // subsidiary DOF type 2 -> master DOF type 4
    // subsidiary DOF type 3 -> master DOF type 5
    //
    //
    //
    // Another example, assume that this block preconditioner has the natural
    // DOF type ordering:
    // 0 1  2 3  4 5  6 7
    // v L1 p L2 u L3 w L4
    //
    // Then the dof_number_in_master_map would be
    // [4 0 6 2], which says: 
    //
    // subsidiary DOF type 0 -> master DOF type 4
    // subsidiary DOF type 1 -> master DOF type 0
    // subsidiary DOF type 2 -> master DOF type 6
    // subsidiary DOF type 3 -> master DOF type 2
    //
    //
    //
    //
    // BACK TO OUR PRECONDITIONER
    //
    // For your convenience, we repeat the natural DOF type ordering of this 
    // block preconditioner here:
    //
    //  0 1 2 3   4  5  6   7   8    9  10 11 12   <- Natural DOF type ordering.
    // [u v w p] [up vp wp Lp1 Lp2] [ut vt wt Lt1] <- DOF type.
    //
    // First, we ignore the fact that there are more DOF types than required by
    // the subsidiary block preconditioner. We give the subsidiary block 
    // preconditioner the following dof_number_in_master_map:
    // [0 1 2 4  5  6  9  10 11 3]
    //  u v w up vp wp ut vt wt p
    // 
    // Which basically says to the subsidiary block preconditioner:
    // "you will work with my DOF types 0, 1, 2, 4, etc.."
    //
    // We chose this DOF ordering as it is the same as our block ordering, and
    // so will be easier when we replace blocks (discussed later).
    // 
    
//    { // encapsulating temp vectors.

      // We create the mapping with the help of two vectors.
      // Noting that the first spatial dimension number of dof types in each mesh
      // corresponds to the velocity dof types, and the rest are either pressure
      // (in the case of the bulk mesh) or lagrange multiplier dof types, we
      // simply:
      // 1) Loop through the meshes
      // 2)   Loop through the spatial_dim, store the dof type in the v_vector.
      // 3)   Loop through the remaining dof types, store the dof type in the l
      //      vector.
      // 4) Concatenate the two vectors.

      //      Artificial test data
      //      n_dof_types = 13;
      //      My_nmesh = 3;
      //      spatial_dim = 3;
      //      My_ndof_types_in_mesh.resize(3,0);
      //      My_ndof_types_in_mesh[0] = 4;
      //      My_ndof_types_in_mesh[1] = 5;
      //      My_ndof_types_in_mesh[2] = 4;
      //      // Outputting My_ndof_types_in_mesh
      //      for (unsigned i = 0; i < My_ndof_types_in_mesh.size(); i++) 
      //      {
      //        std::cout << "My_ndof_inmesh: " 
      //                  << My_ndof_types_in_mesh[i] << std::endl; 
      //      }
      //      pause("Need someone."); 


      unsigned partial_sum_index = 0;

      Vector<unsigned> subsidiary_list_bcpl;
      for(unsigned mesh_i = 0; mesh_i < My_nmesh; mesh_i++)
      {
        // Store the velocity dof types.
        for(unsigned dim_i = 0; dim_i < spatial_dim; dim_i++)
        {
          subsidiary_list_bcpl.push_back(partial_sum_index + dim_i);
        } // for spatial_dim

        // Update the DOF index
        partial_sum_index += My_ndof_types_in_mesh[mesh_i];
      } // for My_nmesh

      // push back the pressure DOF type
      subsidiary_list_bcpl.push_back(spatial_dim);

//    } // end of encapsulating

    //    // Output for artificial test.
    //    std::cout << "subsidiary_list_bcpl:" << std::endl; 
    //    for (unsigned i = 0; i < subsidiary_list_bcpl.size(); i++) 
    //    {
    //      std::cout << subsidiary_list_bcpl[i] << std::endl;
    //    }
    //    // With the artificial test data, this should output:
    //    // 0, 1, 2, 4, 5, 6, 9, 10, 11, 3
    //    pause("Printed out subsidiary_list_bcpl"); 

    // The ns_dof_list will ensure that the NS preconditioner have the 
    // structure:
    // 0  1  2  3  4  5  6
    // ub vb up vp ut vt p
    navier_stokes_block_preconditioner_pt
      ->turn_into_subsidiary_block_preconditioner(
          this, subsidiary_list_bcpl, subsidiary_dof_type_coarsening_map);

    // Set the replacement blocks.
    //
    // Recall: block structure:
    //  0 1 2   3  4  5    6  7  8     9   10  11  12   <- block index
    // [u v w | up vp wp | ut vt wt ] [p | Lp1 Lp2 Lt1] <- DOF type

    //      for (unsigned row_i = spatial_dim; row_i < N_velocity_doftypes; row_i++) 
    //      {
    //        for (unsigned col_i = spatial_dim; col_i < N_velocity_doftypes; col_i++) 
    //        {
    //          navier_stokes_block_preconditioner_pt
    //            ->set_replacement_dof_block(row_i,col_i,
    //                f_subblock_pt(row_i,col_i));
    //        }
    //      }

    //    pause("About to call setup of NS prec."); 



    navier_stokes_block_preconditioner_pt
      ->setup(matrix_pt());

    //      for (unsigned i = 0; i < N_fluid_doftypes; i++)
    //      {
    //        for (unsigned j = 0; j < N_fluid_doftypes; j++) 
    //        {
    //          delete f_subblock_pt(i,j);
    //        }
    //      }
  } // else - NS prec is block preconditioner

  const unsigned v_aug_nrow = v_aug_pt.nrow();
  const unsigned v_aug_ncol = v_aug_pt.ncol();
  for (unsigned v_row = 0; v_row < v_aug_nrow; v_row++) 
  {
    for (unsigned v_col = 0; v_col < v_aug_ncol; v_col++) 
    {
      delete v_aug_pt(v_row,v_col);
      v_aug_pt(v_row,v_col) = 0;
    }
  }
  ///////////////////////////////////////////////////////////////////////////

  // Solver for the W block.
  W_preconditioner_pt.resize(N_lagrange_doftypes,0);
  for(unsigned l_i = 0; l_i < N_lagrange_doftypes; l_i++)
  {
    if(W_preconditioner_fct_pt == 0)
    {
      W_preconditioner_pt[l_i] = new SuperLUPreconditioner;
    }
    else
    {
      // We use the preconditioner provided.
      W_preconditioner_pt[l_i] = 
        (*W_preconditioner_fct_pt)();
    }

    // Is this a block preconditioner?
    BlockPreconditioner<CRDoubleMatrix>* w_block_preconditioner_pt = 
      dynamic_cast<BlockPreconditioner<CRDoubleMatrix>* >
      (W_preconditioner_pt[l_i]);

    if(w_block_preconditioner_pt == 0)
    {
      W_preconditioner_is_block_preconditioner = false;

      W_preconditioner_pt[l_i]->setup(w_pt[l_i]);
    }
    else
    {

      W_preconditioner_is_block_preconditioner = true;
      Vector<unsigned> l_mult_dof_map;
      l_mult_dof_map.push_back(N_fluid_doftypes + l_i);

      W_preconditioner_pt[l_i]->
        turn_into_subsidiary_block_preconditioner(this,l_mult_dof_map);
      W_preconditioner_pt[l_i]->setup(matrix_pt());
    }
  }

  // Delete w_pt(0,N_lagrange_doftypes)
  for (unsigned l_i = 0; l_i < N_lagrange_doftypes; l_i++) 
  {
    delete w_pt[l_i];
  }

  Preconditioner_has_been_setup = true;
} // end of LagrangeEnforcedFlowPreconditioner::setup

/// \short Function to set a new momentum matrix preconditioner 
/// (inexact solver)
void LagrangeEnforcedFlowPreconditioner::set_navier_stokes_preconditioner(
    Preconditioner* new_ns_preconditioner_pt)
{
  // Check if pointer is non-zero.
  if(new_ns_preconditioner_pt == 0)
  {
    std::ostringstream warning_stream;
    warning_stream << "WARNING: \n"
      << "The LSC preconditioner point is null.\n" 
      << "Using default (SuperLU) preconditioner.\n" 
      << std::endl;
    OomphLibWarning(warning_stream.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);

    Navier_stokes_preconditioner_pt = 0;
    Using_superlu_ns_preconditioner = true;
  }
  else
  {
    // If the default SuperLU preconditioner has been used
    // clean it up now...
    if (Using_superlu_ns_preconditioner 
        && Navier_stokes_preconditioner_pt != 0)
    {
      delete Navier_stokes_preconditioner_pt;
    }

    Navier_stokes_preconditioner_pt = new_ns_preconditioner_pt;
    Using_superlu_ns_preconditioner = false;
  }
}

//========================================================================
/// \short Clears the memory.
//========================================================================
void LagrangeEnforcedFlowPreconditioner::clean_up_memory()
{
  // clean the block preconditioner base class memory
  this->clear_block_preconditioner_base();

  // Delete the Navier-Stokes preconditioner pointer.
  //  for now, we do not delete the Navier-Stokes preconditioner pointer.
  // Since we did not create this preconditioner. We should...
  // Create function pointers for both the Lagrange and Navier-Stokes block.
  // At the moment, if SuperLU preconditioning is used for the NS block,
  // There would be a memory leak since we have not deleted this...
  //  delete Navier-Stokes_preconditioner_pt;
  //  Navier-Stokes_preconditioner_pt = 0;

  // Delete the W_preconditioner_pt
  unsigned n_lgr_mult_prec_size = W_preconditioner_pt.size();
  for (unsigned prec_i = 0; prec_i < n_lgr_mult_prec_size; prec_i++) 
  {
    delete W_preconditioner_pt[prec_i];
  }

  // Resize the vector to zero, to avoid any dangling pointers.
  W_preconditioner_pt.resize(0,0);
} // func LagrangeEnforcedFlowPreconditioner::clean_up_memory

}// namespace oomph
