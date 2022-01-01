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
#include "block_preconditioner.h"

namespace oomph
{
  /// Static boolean to allow block_matrix_test(...) to be run.
  /// Defaults to false.
  template<typename MATRIX>
  bool BlockPreconditioner<MATRIX>::Run_block_matrix_test = false;


  //============================================================================
  /// Determine the size of the matrix blocks and setup the
  /// lookup schemes relating the global degrees of freedom with
  /// their "blocks" and their indices (row/column numbers) in those
  /// blocks.
  /// The distributions of the preconditioner and the blocks are
  /// automatically specified (and assumed to be uniform) at this
  /// stage.
  /// This method should be used if any block contains more than one
  /// type of DOF. The argument vector dof_to_block_map should be of length
  /// ndof. Each element should contain an integer indicating the block number
  /// corresponding to that type of DOF.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::block_setup(
    const Vector<unsigned>& dof_to_block_map_in)
  {
#ifdef PARANOID
    // Subsidiary preconditioners don't really need the meshes
    if (this->is_master_block_preconditioner())
    {
      std::ostringstream err_msg;
      unsigned n = nmesh();
      if (n == 0)
      {
        err_msg << "No meshes have been set for this block preconditioner!\n"
                << "Set one with set_nmesh(...), set_mesh(...)" << std::endl;
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        for (unsigned m = 0; m < n; m++)
        {
          if (Mesh_pt[m] == 0)
          {
            err_msg << "The mesh pointer to mesh " << m << " is null!\n"
                    << "Set a non-null one with set_mesh(...)" << std::endl;
            throw OomphLibError(
              err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
          }
        }
      }
    }
#endif

    // Create a copy of the vector input so that we can modify it below
    Vector<unsigned> dof_to_block_map = dof_to_block_map_in;

    if (is_subsidiary_block_preconditioner())
    {
#ifdef PARANOID
      // Get the size of the Doftype_in_master_preconditioner_coarse.
      unsigned para_doftype_in_master_preconditioner_coarse_size =
        Doftype_in_master_preconditioner_coarse.size();

      // Check that the Doftype_in_master_preconditioner_coarse vector is not
      // empty. This must be set (via the function
      // turn_into_subsidiary_block_preconditioner) if this is a
      // subsidiary block preconditioner.
      if (para_doftype_in_master_preconditioner_coarse_size == 0)
      {
        std::ostringstream err_msg;
        err_msg << "The mapping from the dof types of the master "
                << "block preconditioner \n"
                << "to the subsidiary block preconditioner is empty.\n"
                << "Doftype_in_master_preconditioner_coarse.size() == 0 \n"
                << "has turn_into_subsidiary_block_preconditioner(...)\n"
                << "been called with the correct parameters?\n"
                << std::endl;
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }


      // PARANOID checks for Doftype_coarsen_map_coarse
      // This is also set in the function
      // turn_into_subsidiary_block_preconditioner(...).
      //
      // The Doftype_coarsen_map_coarse vector must satisfy two conditions
      // for it to be valid.
      //
      // 1) The dof type numbers in the dof_coarsen_map vector must be
      //    unique. For example, it does not make sense to have the vector
      //    [[0,1][1,2]] because the first inner vector says
      //    "treat dof types 0 and 1 as dof type 0" and the second inner vector
      //    says "treat dof type 1 and 2 as dof type 1", but dof type 1 is
      //    already being treated as dof type 0.
      //
      // 2) Every SUBSIDIARY dof type must be mapped to a dof type in the
      //    Doftype_coarsen_map_coarse vector.
      //    For example, if there are 5 dof types (passed down from the master
      //    block preconditioner), and this block subsidiary block
      //    preconditioner only deals with 3 dof types, then all 5 dof types
      //    must be mapped to a dof type in the subsidiary preconditioner. For
      //    example if the dof_map is [1,2,3,4,5], then the subsidiary block
      //    preconditioner knows that 5 dof types have been passed down. But if
      //    it only works with three dof types, we MUST have three inner vectors
      //    in the doftype_coarsen_map vector (which corresponds to dof types 0,
      //    1 and 2), the union of the dof types in the three inner vectors must
      //    contain dof types 0, 1, 2, 3 and 4 exactly once. It cannot contain,
      //    say, 0, 1, 5, 7, 9, even though it passes the uniqueness check. We
      //    ensure this by two conditions:
      //
      //    2.1) The Doftype_coarsen_map_coarse vector must contain the same
      //         number of dof types as the dof_map vector.
      //         In other words, recall that Doftype_coarsen_map_coarse is a
      //         2D vector, this must contain the same number of vectors as
      //         there are elements in the dof_to_block_map_in vector.
      //
      //    2.2) The maximum element in the doftype_coarsen_map_coarse vector
      //         is the length of the dof_map vector minus 1.

      // A set is deal for checking the above three conditions, we shall insert
      // all the elements in the doftype_coarsen_map_coarse vector into this
      // set.
      std::set<unsigned> doftype_map_set;

      // Condition (1): Check for uniqueness by inserting all the values of
      // Doftype_coarsen_map_coarse into a set.
      unsigned para_doftype_coarsen_map_coarse_size =
        Doftype_coarsen_map_coarse.size();

      // Loop through the outer vector of Doftype_coarsen_map_coarse
      // then loop through the inner vectors and attempt to insert each
      // element of Doftype_coarsen_map_coarse into doftype_map_set.
      //
      // The inner for loop will throw an error if we cannot insert the
      // element, this means that it is already inserted and thus not unique.
      for (unsigned i = 0; i < para_doftype_coarsen_map_coarse_size; i++)
      {
        // Loop through the inner vector
        unsigned para_doftype_coarsen_map_coarse_i_size =
          Doftype_coarsen_map_coarse[i].size();
        for (unsigned j = 0; j < para_doftype_coarsen_map_coarse_i_size; j++)
        {
          // Attempt to insert all the values of the inner vector into a set.
          std::pair<std::set<unsigned>::iterator, bool> doftype_map_ret =
            doftype_map_set.insert(Doftype_coarsen_map_coarse[i][j]);

          if (!doftype_map_ret.second)
          {
            std::ostringstream err_msg;
            err_msg << "Error: the doftype number "
                    << Doftype_coarsen_map_coarse[i][j]
                    << " is already inserted." << std::endl;
            throw OomphLibError(
              err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

      // Condition (2.1): Check that the doftype_map_set describes as many
      // values as doftype_in_master_preconditioner_coarse. I.e. if dof_map
      // contains 5 dof types, then the doftype_coarsen_map_coarse vector must
      // also contain 5 dof types.
      if (para_doftype_in_master_preconditioner_coarse_size !=
          doftype_map_set.size())
      {
        std::ostringstream err_msg;
        err_msg << "The size of doftype_in_master_preconditioner_coarse "
                << "must be the same as the total\n"
                << "number of values in the doftype_coarsen_map_coarse vector."
                << std::endl;
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Condition (2.2): Check that the maximum element in the
      // doftype_coarsen_map_coarse vector is the length of the
      // doftype_in_master_preconditioner_coarse minus 1.
      unsigned para_doftype_in_master_preconditioner_coarse_size_minus_one =
        para_doftype_in_master_preconditioner_coarse_size - 1;
      if (para_doftype_in_master_preconditioner_coarse_size_minus_one !=
          *doftype_map_set.rbegin())
      {
        std::ostringstream err_msg;
        err_msg << "The maximum dof type number in the "
                << "doftype_coarsen_map vector must be "
                << para_doftype_in_master_preconditioner_coarse_size_minus_one
                << std::endl;
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Set the mapping from the master preconditioner DOF types to the
      // subsidiary preconditioner DOF types.
      //
      // IMPORTANT: Since DOF types may be coarsened in the master block
      // preconditioner, this may no longer reflect the actual underlying dof
      // types. We must get the actual underlying dof types for the
      // block_setup(...) function to work properly so all the look up schemes
      // for this (subsidiary) block preconditioner is correct and works
      // properly, this is for backwards compatibility purposes and to make sure
      // Richard Muddle's still works at this (subsidiary) level, although it
      // may not be used.
      //
      // If we do not want to make it backwards compatible, we may as well
      // kill the block_setup(...) for subsidiary block preconditioners -
      // but other thing may break. Do it at your own risk (take time to
      // fully understand the whole block preconditioning framework code).

      // Create the corresponding Doftype_in_master_preconditioner_fine and
      // Doftype_coarsen_map_fine vectors.

      // First resize the vectors.
      Doftype_in_master_preconditioner_fine.resize(0);
      Doftype_coarsen_map_fine.resize(0);

      // The Doftype_in_master_preconditioner_fine vector is easy.  We know that
      // the Doftype_coarsen_map_fine in the master preconditioner must be
      // constructed already. So we simply loop through the values in
      // doftype_in_master_preconditioner_coarse, then get the most fine grain
      // dof types from the master preconditioner's Doftype_coarsen_map_fine
      // vector.
      //
      // For example, if the master preconditioner has the vector:
      // Doftype_coarsen_map_fine = [0,1,2,3][4,5,6,7][8,9,10,11][12,13][14,15]
      //
      // and passes the two vectors
      // doftype_in_master_preconditioner_coarse = [1,2,3]
      // doftype_coarsen_map_coarse = [[0][1,2]]
      //
      // Then we want
      // Doftype_in_master_preconditioner_fine = [4,5,6,7,8,9,10,11,12,13]
      //
      // We achieve this by looking up the corresponding fine dof types in the
      // masters' Doftype_coarsen_map_fine vector which corresponds to the
      // values in Doftype_in_master_preconditioner_coarse.
      //
      // That is, the values in Doftype_in_master_preconditioner_coarse gives us
      // the index of sub vector we want in the master's
      // Doftype_coarsen_map_fine vector.

#ifdef PARANOID
      // Check that the master block preconditioner's Doftype_coarsen_map_fine
      // is set up. Under the current implementation, this would always be set
      // up properly, but we check it just in case!
      if (master_block_preconditioner_pt()->doftype_coarsen_map_fine().size() ==
          0)
      {
        std::ostringstream err_msg;
        err_msg << "The master block preconditioner's "
                << "Doftype_coarsen_map_fine is not\n"
                << "set up properly.\n"
                << "\n"
                << "This vector is constructed in the function "
                << "block_setup(...).\n"
                << std::endl;
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      unsigned doftype_in_master_preconditioner_coarse_size =
        Doftype_in_master_preconditioner_coarse.size();
      for (unsigned i = 0; i < doftype_in_master_preconditioner_coarse_size;
           i++)
      {
        // The index of the sub vector we want.
        unsigned subvec_index = Doftype_in_master_preconditioner_coarse[i];

        // Get the corresponding most fine grain sub vector from the master
        // block preconditioner
        Vector<unsigned> tmp_master_dof_subvec =
          Master_block_preconditioner_pt->get_fine_grain_dof_types_in(
            subvec_index);

        Doftype_in_master_preconditioner_fine.insert(
          Doftype_in_master_preconditioner_fine.end(),
          tmp_master_dof_subvec.begin(),
          tmp_master_dof_subvec.end());
      }

      // The Doftype_coarsen_map_fine vector is a bit more tricky.
      // The Doftype_coarsen_map_coarse vector describes which coarse dof types
      // of THIS preconditioner are grouped together. We have to translate this
      // into the most fine grain dof types.
      //
      // For example, if
      // Doftype_coarsen_map_coarse = [[0][1,2]]
      // Doftype_in_master_preconditioner_coarse = [1,2,3]
      //
      // and the MASTER preconditioner has:
      // Doftype_coarsen_map_fine= [[0,1,2,3][4,5,6,7][8,9,10,11][12,13][14,15]]
      //
      // Then [[0][1,2]] tell us that the most fine grain DOF types 1 of the
      // master preconditioner most be grouped together, and the most fine
      // grained dof types 2 and 3 of the master preconditioner must be grouped
      // together.
      //
      // This gives the vector [[4,5,6,7] [8,9,10,11,12,13]], translating this
      // into the local DOF types of this preconditioner we have
      // Doftype_coarsen_map_fine = [[0,1,2,3][4,5,6,7,8,9]]. This corresponds
      // with the Doftype_in_master_preconditioner_fine vector we created above:
      // Doftype_in_master_preconditioner_fine = [4,5,6,7,8,9,10,11,12,13]
      //
      // Together, the master block preconditioner says to THIS subsidiary block
      // preconditioner "work on my DOF types [4,5,6,7,8,9,10,11,12,13], but
      // group your DOF type [0,1,2,3] together as DOF type 0 and [4,5,6,7,8,9]
      // together together as DOF type 1".
      //
      // Think of it like this: For each DOF type in Doftype_coarsen_map_coarse
      // we look at how many values this corresponds to in the master
      // preconditioner. In this case, Doftype_coarsen_map_coarse:
      //
      // 1 - corresponds to fine DOF types 0,1,2,3 in this preconditioner,
      // and 4,5,6,7 in the master preconditioner;
      //
      // 2 - corresponds to fine DOF types 4,5,6,7 in this preconditioner,
      // and 8,9,10,11 in the master preconditioner;
      //
      // 3 - corresponds to fine DOF types 8,9 in this preconditioner,
      // and 12,13 in the master preconditioner.
      //
      // Thus Doftype_coarsen_map_fine = [[0,1,2,3][4,5,6,7,8,9]]
      //
      /// /////////////////////////////////////////////////////////////////////
      //
      // How to do this: First we create a 2D vector which has the corresponds
      // to the fine dof types in the master preconditioner but starting from
      // 0. For example, take the above example (repeated below):
      //   Passed to this prec by the master prec:
      //   Doftype_coarsen_map_coarse = [[0][1,2]]
      //   Doftype_in_master_preconditioner_coarse = [1,2,3]
      //
      // and the MASTER preconditioner has:
      // Doftype_coarsen_map_fine= [[0,1,2,3][4,5,6,7][8,9,10,11][12,13][14,15]]
      //
      // Step 1:
      // Then, the temp 2D vector we want to create is:
      // master_fine_doftype_translated = [[0 1 2 3], [4,5,6,7], [8,9]]
      // This comes from using Doftype_in_master_preconditioner_coarse
      // then get the number of fine dof types in the master.
      //
      // Step 2:
      // Then:
      //   Loop through the vector Doftype_coarsen_map_coarse,
      //     Loop over the inner vectors in Doftype_coarsen_map_coarse
      //       Each element in the inner vector corresponds to a vector in
      //       master_fine_doftype_translated. We push in the vectors of
      //       master_fine_doftype_translated intp Doftype_coarsen_map_fine
      //

      Vector<Vector<unsigned>> master_fine_doftype_translated;
      unsigned dof_type_index = 0;
      for (unsigned i = 0; i < doftype_in_master_preconditioner_coarse_size;
           i++)
      {
        // How many fine DOF types are in the master's
        // Doftype_in_master_preconditioner_coarse[i]?
        unsigned coarse_dof = Doftype_in_master_preconditioner_coarse[i];

        unsigned n_master_fine_doftypes =
          Master_block_preconditioner_pt->nfine_grain_dof_types_in(coarse_dof);

        Vector<unsigned> tmp_sub_vec;
        for (unsigned j = 0; j < n_master_fine_doftypes; j++)
        {
          tmp_sub_vec.push_back(dof_type_index);
          dof_type_index++;
        }
        master_fine_doftype_translated.push_back(tmp_sub_vec);
      }


      // master_fine_doftype_translated now contains vectors with values are
      // from 0, 1, 2, ..,
      //
      // Now read out the values of master_fine_doftype_translated and place
      // them in order according to Doftype_coarsen_map_coarse.
      unsigned doftype_coarsen_map_coarse_size =
        Doftype_coarsen_map_coarse.size();
      for (unsigned i = 0; i < doftype_coarsen_map_coarse_size; i++)
      {
        Vector<unsigned> tmp_vec;
        unsigned doftype_coarsen_map_coarse_i_size =
          Doftype_coarsen_map_coarse[i].size();
        for (unsigned j = 0; j < doftype_coarsen_map_coarse_i_size; j++)
        {
          unsigned subvec_i = Doftype_coarsen_map_coarse[i][j];

          tmp_vec.insert(tmp_vec.end(),
                         master_fine_doftype_translated[subvec_i].begin(),
                         master_fine_doftype_translated[subvec_i].end());
        }

        Doftype_coarsen_map_fine.push_back(tmp_vec);
      }

      // Get the number of block types (and DOF types) in this preconditioner
      // from the length of the dof_map vector.
      Internal_ndof_types = Doftype_in_master_preconditioner_fine.size();

      // Nblock_types is later updated in block_setup(...)
      Internal_nblock_types = Internal_ndof_types;

      // Compute number of rows in this (sub) preconditioner using data from
      // the master.
      Nrow = 0;
      for (unsigned b = 0; b < Internal_ndof_types; b++)
      {
        Nrow += this->internal_dof_block_dimension(b);
      }

#ifdef PARANOID
      if (Nrow == 0)
      {
        std::ostringstream error_message;
        error_message
          << "Nrow=0 in subsidiary preconditioner. This seems fishy and\n"
          << "suggests that block_setup() was not called for the \n"
          << "master block preconditioner yet.";
        throw OomphLibWarning(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    // If this is a master block preconditioner, then set the
    // Doftype_coarsen_map_fine and Doftype_coarsen_map_coarse to the
    // identity. Recall that the Doftype_coarsen_map_fine maps the dof types
    // that this preconditioner requires with the most fine grain dof types (the
    // internal dof types) and the Doftype_coarsen_map_coarse maps the dof
    // types that this preconditioner requires with the dof types which this
    // preconditioner is given from a master preconditioner (these dof types may
    // or may not be coarsened). In the case of the master preconditioner, these
    // are the same (since dof types are not coarsened), furthermore the
    // identity mapping is provided to say that dof type 0 maps to dof type 0,
    // dof type 1 maps to dof type 1,
    // dof type 2 maps to dof type 2,
    // etc...
    //
    // If this is not a master block preconditioner, then the vectors
    // Doftype_coarsen_map_fine and Doftype_coarsen_map_coarse is handled
    // by the turn_into_subsidiary_block_preconditioner(...) function.
    if (is_master_block_preconditioner())
    {
      // How many dof types does this preconditioner work with?
      unsigned n_external_dof_types = dof_to_block_map.size();

      // Note: at the master level, the n_external_dof_types should be the same
      // as the internal_ndof_types(), since the dof_to_block_map MUST describe
      // the mapping between every dof type (not yet coarsened - so it is the
      // same number as the internal dof types) to the block types. But we
      // distinguish them for clarity. We also check that this is the case.
#ifdef PARANOID
      unsigned n_internal_dof_types = internal_ndof_types();

      if (n_internal_dof_types != n_external_dof_types)
      {
        std::ostringstream err_msg;
        err_msg
          << "The internal ndof types and the length of the dof_to_block_map\n"
          << "vector is not the same. Since this is the master block "
          << "preconditioner,\n"
          << "you must describe which block each DOF type belongs to,\n"
          << "no more, no less."
          << "internal_ndof_types = " << n_internal_dof_types << "\n"
          << "dof_to_block_map.size() = " << n_external_dof_types << "\n";
        throw OomphLibWarning(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Clear and reserve space.
      Doftype_coarsen_map_fine.clear();
      Doftype_coarsen_map_coarse.clear();
      Doftype_coarsen_map_fine.reserve(n_external_dof_types);
      Doftype_coarsen_map_coarse.reserve(n_external_dof_types);

      // Now push back the identity mapping.
      for (unsigned i = 0; i < n_external_dof_types; i++)
      {
        // Create a vector and push it in.
        Vector<unsigned> tmp_vec(1, i);
        Doftype_coarsen_map_fine.push_back(tmp_vec);
        Doftype_coarsen_map_coarse.push_back(tmp_vec);
      }
    }
    else
    // Else this is a subsidiary block preconditioner.
    {
      // Both the Doftype_coarsen_map_fine and Doftype_coarsen_map_coarse
      // vectors must be already be handled by the
      // turn_into_subsidiary_block_preconditioner(...) function. We check this.
#ifdef PARANOID
      if ((Doftype_coarsen_map_fine.size() == 0) ||
          (Doftype_coarsen_map_coarse.size() == 0))
      {
        std::ostringstream err_msg;
        err_msg << "Either the Doftype_coarsen_map_fine or the \n"
                << "Doftype_coarsen_map_coarse vectors is of size 0.\n"
                << "Did you remember to call the function "
                << "turn_into_subsidiary_block_preconditioner(...)?";
        throw OomphLibWarning(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }


    // Now we create the vector Block_to_dof_map_coarse.
    // Recall that the vector describe which dof types are in which block with
    // the relationship:
    //
    // Block_to_dof_map_coarse[block_number] = Vector[dof types];
    //
    // Note that this is not the internal (underlying) dof type.
    // Nor is this in relation to the parent block preconditioner's dof type.
    // The number of elements in it is the same as dof_to_block_map vector.
    //
    // Since the dof type coarsening feature is added later, we encapsulate this
    // bit of the code so it does not affect things below.
    {
      // Check that the dof_to_block map "makes sense" for the
      // Doftype_coarsen_map_coarse.
      // The Doftype_coarsen_map_coarse describes which doftypes should be
      // considered as a single doftype in THIS preconditioner.
      //
      // For example, if this preconditioner is the LSC block preconditioner
      // applied to a 3D problem, it expects 4 doftypes:
      // 3 velocity, [u, v, w] and 1 pressure [p],
      // giving us the dof type ordering
      // [u v w p].
      //
      // The LSC preconditioner groups the velocity and pressure doftypes
      // separately, thus the dof_to_block_map will be:
      // [0 0 0 1]
      //
      // Then the Doftype_coarsen_map_coarse MUST have length 4, to identify
      // which of the OTHER (possibly coarsened) dof types should be grouped
      // together to be considered as a single dof types of THIS preconditioner.
      //
      // For example, if the preconditioner above this one has the dof type
      // ordering:
      // 0  1  2  3  4  5  6  7  8  9
      // ub vb wb up vp wp ut vt wt p
      // Then we want to tell THIS preconditioner which dof types belongs to
      // u, v, w and p, by providing the optional argument
      // Doftype_coarsen_map_coarse to the
      // turn_into_subsidiary_block_preconditioner(...) function
      // [[0 3 6] [1 4 7] [2 5 8] [9]]
      //
      // So, it is important that the length of dof_to_block_map is the same as
      // the length of Doftype_coarsen_map_coarse. We check this.
      unsigned dof_to_block_map_size = dof_to_block_map.size();

#ifdef PARANOID
      if (dof_to_block_map_size != Doftype_coarsen_map_coarse.size())
      {
        std::ostringstream err_msg;
        err_msg
          << "The size of dof_to_block_map and Doftype_coarsen_map_coarse is "
             "not "
          << "the same.\n"
          << "dof_to_block_map.size() = " << dof_to_block_map_size << "\n"
          << "Doftype_coarsen_map_coarse.size() = "
          << Doftype_coarsen_map_coarse.size() << ".\n"
          << "One of the two list is incorrect, please look at the comments\n"
          << "in the source code for more details.";
        throw OomphLibWarning(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Create the Block_to_dof_map_coarse from
      // the dof_to_block_map and Doftype_coarsen_map_coarse.

      // find the maximum block number
      unsigned max_block_number =
        *std::max_element(dof_to_block_map.begin(), dof_to_block_map.end());

      // Now we do the following:
      // Lets say the Doftype_coarsen_map_coarse is:
      // [0 3 6]
      // [1 4 7]
      // [2 5 8]
      // [9]
      //
      // (this is the same as the above example)
      //
      // and the dof_to_block_map is [0 0 0 1].
      //
      // Then we need to form the Block_to_dof_map_coarse:
      // [0 3 6 1 4 7 2 5 8]
      // [9]

      // Clear anything in the Block_to_dof_map_coarse
      Block_to_dof_map_coarse.clear();

      const unsigned tmp_nblock = max_block_number + 1;

      Block_to_dof_map_coarse.resize(tmp_nblock);

      for (unsigned i = 0; i < dof_to_block_map_size; i++)
      {
        Block_to_dof_map_coarse[dof_to_block_map[i]].push_back(i);
      }

      Block_to_dof_map_fine.clear();
      Block_to_dof_map_fine.resize(tmp_nblock);
      for (unsigned block_i = 0; block_i < tmp_nblock; block_i++)
      {
        // get the dof types in this block.
        const unsigned ndof_in_block = Block_to_dof_map_coarse[block_i].size();
        for (unsigned dof_i = 0; dof_i < ndof_in_block; dof_i++)
        {
          const unsigned coarsened_dof_i =
            Block_to_dof_map_coarse[block_i][dof_i];

          // Insert the most fine grain dofs which this dof_i corresponds to
          // into block_i
          Vector<unsigned> dof_i_dofs =
            Doftype_coarsen_map_fine[coarsened_dof_i];

          Block_to_dof_map_fine[block_i].insert(
            Block_to_dof_map_fine[block_i].end(),
            dof_i_dofs.begin(),
            dof_i_dofs.end());
        }
      }

      // Now set the dof_to_block_map to the identify.
      // NOTE: We are now using the internal n dof types. This is because the
      // dof type coarsening feature was built on top of the existing block
      // preconditioning framework which does not handle coarsening of dof
      // types. Hence, under the hood, it still works with the most fine grain
      // dof types and does not do any coarsening.

      // Locally cache the internal ndof types (using access function because
      // the Internal_ndof_type variable may not be set up yet if this is a
      // master preconditioner).
      unsigned tmp_internal_ndof_types = internal_ndof_types();

      dof_to_block_map.resize(tmp_internal_ndof_types, 0);

      for (unsigned i = 0; i < tmp_internal_ndof_types; i++)
      {
        dof_to_block_map[i] = i;
      }
    } // end of Block_to_dof_map_coarse encapsulation

#ifdef PARANOID

    // Check that the meshes are ok. This only needs to be done in the master
    // because subsidiary preconditioners don't do anything with the meshes
    // here.
    if (is_master_block_preconditioner())
    {
      // This is declared as local_nmesh because there are other variables
      // called nmesh floating about... but this will not exist if PARANOID is
      // switched on.
      unsigned local_nmesh = nmesh();

      // Check that some mesh pointers have been assigned.
      if (local_nmesh == 0)
      {
        std::ostringstream error_msg;
        error_msg << "Cannot setup blocks because no meshes have been set.";
        throw OomphLibError(
          error_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Each mesh must contain elements with the same number of dof.
      // A stricter check is to ensure that the mesh contains only one type of
      // elements. This is relaxed in same cases.
      for (unsigned mesh_i = 0; mesh_i < local_nmesh; mesh_i++)
      {
        // The number of elements in the current mesh.
        unsigned n_element = mesh_pt(mesh_i)->nelement();

        // When the bulk mesh is distributed, there may not be any elements
        // in the surface mesh(es).
        if (n_element > 0)
        {
          // The string of the first element in the current mesh.
          std::string first_element_string =
            typeid(*(mesh_pt(mesh_i)->element_pt(0))).name();

          // If there are multiple element types in the current mesh,
          // we can at least make sure that they contain the same types of DOFs.
          if (bool(Allow_multiple_element_type_in_mesh[mesh_i]))
          {
            // The ndof types of the first element.
            unsigned first_element_ndof_type =
              mesh_pt(mesh_i)->element_pt(0)->ndof_types();

            // Loop through the meshes and compare the number of types of DOFs.
            for (unsigned el_i = 1; el_i < n_element; el_i++)
            {
              // The ndof type of the current element.
              unsigned current_element_ndof_type =
                mesh_pt(mesh_i)->element_pt(el_i)->ndof_types();

              // The string of the current element.
              std::string current_element_string =
                typeid(*(mesh_pt(mesh_i)->element_pt(el_i))).name();

              // Compare against the first element.
              if (current_element_ndof_type != first_element_ndof_type)
              {
                std::ostringstream error_message;
                error_message
                  << "Elements in the same mesh MUST have the same number of "
                     "types "
                  << "of DOFs.\n"
                  << "The element in mesh " << mesh_i << ", at position "
                  << el_i << " is: \n"
                  << current_element_string << ", \n"
                  << "with ndof types: " << current_element_ndof_type << ".\n"
                  << "The first element in the same mesh is: \n"
                  << first_element_string << ", \n"
                  << "with ndof types: " << first_element_ndof_type << ".\n";
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
            }
          }
          else
          // There should be only one type of elements in the current mesh.
          // Check that this is the case!
          {
            // Loop through the elements in the current mesh.
            for (unsigned el_i = 1; el_i < n_element; el_i++)
            {
              // The string of the current element.
              std::string current_element_string =
                typeid(*(mesh_pt(mesh_i)->element_pt(el_i))).name();

              // Compare against the first element.
              if (current_element_string.compare(first_element_string) != 0)
              {
                std::ostringstream error_message;
                error_message
                  << "By default, a mesh containing block preconditionable "
                  << "elements must contain only one type of element.\n"
                  << "The element in mesh " << mesh_i << ", at position "
                  << el_i << " is: \n"
                  << current_element_string << "\n"
                  << "The first element in the same mesh is: \n"
                  << first_element_string << "\n"
                  << "If this is correct, consider calling the set_mesh(...) "
                     "with\n"
                  << "the optional argument set true to allow multiple "
                     "element\n"
                  << "types in the same mesh.\n"
                  << "Note: A minimal requirement is that the elements in the "
                     "same\n"
                  << "mesh MUST have the same number of DOF types.";
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
            }
          }
        }
      }
    }

#endif
    // clear the memory
    this->clear_block_preconditioner_base();

    // get my_rank and nproc
#ifdef OOMPH_HAS_MPI
    unsigned my_rank = comm_pt()->my_rank();
    unsigned nproc = comm_pt()->nproc();
#endif


    /// //////////////////////////////////////////////////////////////////////////
    // start of master block preconditioner only operations
    /// //////////////////////////////////////////////////////////////////////////
#ifdef OOMPH_HAS_MPI
    unsigned* nreq_sparse = new unsigned[nproc]();
    unsigned* nreq_sparse_for_proc = new unsigned[nproc]();
    unsigned** index_in_dof_block_sparse_send = new unsigned*[nproc]();
    unsigned** dof_number_sparse_send = new unsigned*[nproc]();
    Vector<MPI_Request> send_requests_sparse;
    Vector<MPI_Request> recv_requests_sparse;
#endif

    // If this preconditioner is the master preconditioner then we need
    // to assemble the vectors : Dof_number
    //                           Index_in_dof_block
    if (is_master_block_preconditioner())
    {
      // Get the number of dof types in each mesh.
      Ndof_types_in_mesh.assign(nmesh(), 0);
      for (unsigned i = 0; i < nmesh(); i++)
      {
        Ndof_types_in_mesh[i] = mesh_pt(i)->ndof_types();
      }
      // Setup the distribution of this preconditioner, assumed to be the same
      // as the matrix if the matrix is distributable.
      if (dynamic_cast<DistributableLinearAlgebraObject*>(matrix_pt()))
      {
        this->build_distribution(
          dynamic_cast<DistributableLinearAlgebraObject*>(matrix_pt())
            ->distribution_pt());
      }
      else
      {
        LinearAlgebraDistribution dist(comm_pt(), matrix_pt()->nrow(), false);
        this->build_distribution(dist);
      }
      Nrow = matrix_pt()->nrow();

      // Boolean to indicate whether the matrix is actually distributed,
      // ie distributed and on more than one processor.
      bool matrix_distributed =
        (this->distribution_pt()->distributed() &&
         this->distribution_pt()->communicator_pt()->nproc() > 1);


      // Matrix must be a CR matrix.
      CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt());

      if (cr_matrix_pt == 0)
      {
        std::ostringstream error_message;
        error_message << "Block setup for distributed matrices only works "
                      << "for CRDoubleMatrices";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }


      // Get distribution.
      unsigned first_row = this->distribution_pt()->first_row();
      unsigned nrow_local = this->distribution_pt()->nrow_local();
      unsigned last_row = first_row + nrow_local - 1;

#ifdef OOMPH_HAS_MPI
      // storage for the rows required by each processor in the dense
      // block lookup storage scheme
      // dense_required_rows(p,0) is the minimum global index required by proc p
      //                 ...(p,1) is the maximum global index required by proc p
      DenseMatrix<unsigned> dense_required_rows(nproc, 2);
      for (unsigned p = 0; p < nproc; p++)
      {
        dense_required_rows(p, 0) = this->distribution_pt()->first_row(p);
        dense_required_rows(p, 1) = this->distribution_pt()->first_row(p) +
                                    this->distribution_pt()->nrow_local(p) - 1;
      }

      // determine the global rows That are not in the range first_row to
      // first_row+nrow_local for which we should store the
      // Dof_index and Index_in_dof_block for
      // then send the lists to other processors
      std::set<unsigned> sparse_global_rows_for_block_lookup;
      if (matrix_distributed)
      {
        unsigned nnz = cr_matrix_pt->nnz();
        int* column_index = cr_matrix_pt->column_index();
        for (unsigned i = 0; i < nnz; i++)
        {
          unsigned ci = column_index[i];
          if (ci < first_row || ci > last_row)
          {
            sparse_global_rows_for_block_lookup.insert(ci);
          }
        }
      }

      int nsparse = sparse_global_rows_for_block_lookup.size();

      Global_index_sparse.resize(0);
      std::copy(sparse_global_rows_for_block_lookup.begin(),
                sparse_global_rows_for_block_lookup.end(),
                std::back_inserter(Global_index_sparse));

      Index_in_dof_block_sparse.resize(nsparse);
      Dof_number_sparse.resize(nsparse);
      sparse_global_rows_for_block_lookup.clear();

      Vector<MPI_Request> recv_requests_sparse_nreq;
      if (matrix_distributed)
      {
        MPI_Aint base_displacement_sparse;
        MPI_Get_address(nreq_sparse, &base_displacement_sparse);

        int zero = 0;
        for (unsigned p = 0; p < nproc; p++)
        {
          // determine the global eqn numbers required by this processor
          // that can be classified by processor p
          int begin = 0;
          for (int i = 0; i < nsparse; ++i)
          {
            if (Global_index_sparse[i] < dense_required_rows(p, 0))
            {
              ++begin;
            }
            else
            {
              if (Global_index_sparse[i] <= dense_required_rows(p, 1))
              {
                ++nreq_sparse[p];
              }
              else
              {
                break;
              }
            }
          }

          // if this processor has rows to be classified by proc p
          if (nreq_sparse[p] > 0)
          {
            // send the number of global eqn numbers
            MPI_Request req1;
            MPI_Isend(&nreq_sparse[p],
                      1,
                      MPI_UNSIGNED,
                      p,
                      31,
                      comm_pt()->mpi_comm(),
                      &req1);
            send_requests_sparse.push_back(req1);

            // send the global eqn numbers
            MPI_Request req2;
            MPI_Isend(&Global_index_sparse[begin],
                      nreq_sparse[p],
                      MPI_UNSIGNED,
                      p,
                      32,
                      comm_pt()->mpi_comm(),
                      &req2);
            send_requests_sparse.push_back(req2);

            // post the recvs for the data that will be returned

            // the datatypes, displacements, lengths for the two datatypes
            MPI_Datatype types[2];
            MPI_Aint displacements[2];
            int lengths[2];

            // index in dof block
            MPI_Type_contiguous(nreq_sparse[p], MPI_UNSIGNED, &types[0]);
            MPI_Type_commit(&types[0]);
            MPI_Get_address(&Index_in_dof_block_sparse[begin],
                            &displacements[0]);
            displacements[0] -= base_displacement_sparse;
            lengths[0] = 1;

            // dof number
            MPI_Type_contiguous(nreq_sparse[p], MPI_UNSIGNED, &types[1]);
            MPI_Type_commit(&types[1]);
            MPI_Get_address(&Dof_number_sparse[begin], &displacements[1]);
            displacements[1] -= base_displacement_sparse;
            lengths[1] = 1;

            // build the final type
            MPI_Datatype recv_type;
            MPI_Type_create_struct(
              2, lengths, displacements, types, &recv_type);
            MPI_Type_commit(&recv_type);
            MPI_Type_free(&types[0]);
            MPI_Type_free(&types[1]);

            // and recv
            MPI_Request req;
            MPI_Irecv(
              nreq_sparse, 1, recv_type, p, 33, comm_pt()->mpi_comm(), &req);
            recv_requests_sparse.push_back(req);
            MPI_Type_free(&recv_type);
          }

          // if no communication required, confirm this
          if (nreq_sparse[p] == 0)
          {
            MPI_Request req1;
            MPI_Isend(
              &zero, 1, MPI_UNSIGNED, p, 31, comm_pt()->mpi_comm(), &req1);
            send_requests_sparse.push_back(req1);
          }

          //
          MPI_Request req;
          MPI_Irecv(&nreq_sparse_for_proc[p],
                    1,
                    MPI_UNSIGNED,
                    p,
                    31,
                    comm_pt()->mpi_comm(),
                    &req);
          recv_requests_sparse_nreq.push_back(req);
        }
      }
#endif

      // resize the storage
      Dof_number_dense.resize(nrow_local);
      Index_in_dof_block_dense.resize(nrow_local);

      // zero the number of dof types
      Internal_ndof_types = 0;

#ifdef PARANOID
      // Vector to keep track of previously assigned block numbers
      // to check consistency between multiple assignments.
      Vector<int> previously_assigned_block_number(nrow_local,
                                                   Data::Is_unclassified);
#endif

      // determine whether the problem is distribution
      bool problem_distributed = false;

      // the problem method distributed() is only accessible with MPI
#ifdef OOMPH_HAS_MPI
      problem_distributed = any_mesh_distributed();
#endif

      // if the problem is not distributed
      if (!problem_distributed)
      {
        // Offset for the block type in the overall system.
        // Different meshes contain different block-preconditionable
        // elements -- their blocks are added one after the other.
        unsigned dof_offset = 0;

        // Loop over all meshes.
        for (unsigned m = 0; m < nmesh(); m++)
        {
          // Number of elements in this mesh.
          unsigned n_element = mesh_pt(m)->nelement();

          // Find the number of block types that the elements in this mesh
          // are in charge of.
          unsigned ndof_in_element = ndof_types_in_mesh(m);
          Internal_ndof_types += ndof_in_element;

          for (unsigned e = 0; e < n_element; e++)
          {
            // List containing pairs of global equation number and
            // dof number for each global dof in an element.
            std::list<std::pair<unsigned long, unsigned>> dof_lookup_list;

            // Get list of blocks associated with the element's global unknowns.
            mesh_pt(m)->element_pt(e)->get_dof_numbers_for_unknowns(
              dof_lookup_list);

            // Loop over all entries in the list
            // and store the block number.
            typedef std::list<std::pair<unsigned long, unsigned>>::iterator IT;
            for (IT it = dof_lookup_list.begin(); it != dof_lookup_list.end();
                 it++)
            {
              unsigned long global_dof = it->first;
              if (global_dof >= unsigned(first_row) &&
                  global_dof <= unsigned(last_row))
              {
                unsigned dof_number = (it->second) + dof_offset;
                Dof_number_dense[global_dof - first_row] = dof_number;

#ifdef PARANOID
                // Check consistency of block numbers if assigned multiple times
                if (previously_assigned_block_number[global_dof - first_row] <
                    0)
                {
                  previously_assigned_block_number[global_dof - first_row] =
                    dof_number;
                }
#endif
              }
            }
          }

          // About to do the next mesh which contains block preconditionable
          // elements of a different type; all the dofs that these elements are
          // "in charge of" differ from the ones considered so far.
          // Bump up the block counter to make sure we're not overwriting
          // anything here
          dof_offset += ndof_in_element;
        }

#ifdef PARANOID
        // check that every global equation number has been allocated a dof type
        for (unsigned i = 0; i < nrow_local; i++)
        {
          if (previously_assigned_block_number[i] < 0)
          {
            std::ostringstream error_message;
            error_message << "Not all degrees of freedom have had DOF type "
                          << "numbers allocated. Dof number " << i
                          << " is unallocated.";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
#endif
      }
      // else the problem is distributed
      else
      {
#ifdef OOMPH_HAS_MPI
        // Offset for the block type in the overall system.
        // Different meshes contain different block-preconditionable
        // elements -- their blocks are added one after the other...
        unsigned dof_offset = 0;

        // the set of global degrees of freedom and their corresponding dof
        // number on this processor
        std::map<unsigned long, unsigned> my_dof_map;

        // Loop over all meshes
        for (unsigned m = 0; m < nmesh(); m++)
        {
          // Number of elements in this mesh
          unsigned n_element = this->mesh_pt(m)->nelement();

          // Find the number of block types that the elements in this mesh
          // are in charge of
          unsigned ndof_in_element = ndof_types_in_mesh(m);
          Internal_ndof_types += ndof_in_element;

          // Loop over all elements
          for (unsigned e = 0; e < n_element; e++)
          {
            // if the element is not a halo element
            if (!this->mesh_pt(m)->element_pt(e)->is_halo())
            {
              // List containing pairs of global equation number and
              // dof number for each global dof in an element
              std::list<std::pair<unsigned long, unsigned>> dof_lookup_list;

              // Get list of blocks associated with the element's global
              // unknowns
              this->mesh_pt(m)->element_pt(e)->get_dof_numbers_for_unknowns(
                dof_lookup_list);

              // update the block numbers and put it in the map.
              typedef std::list<std::pair<unsigned long, unsigned>>::iterator
                IT;
              for (IT it = dof_lookup_list.begin(); it != dof_lookup_list.end();
                   it++)
              {
                it->second = (it->second) + dof_offset;
                my_dof_map[it->first] = it->second;
              }
            }
          }

          // About to do the next mesh which contains block preconditionable
          // elements of a different type; all the dofs that these elements are
          // "in charge of" differ from the ones considered so far.
          // Bump up the block counter to make sure we're not overwriting
          // anything here
          dof_offset += ndof_in_element;
        }

        // next copy the map of my dofs to two vectors to send
        unsigned my_ndof = my_dof_map.size();
        unsigned long* my_global_dofs = new unsigned long[my_ndof];
        unsigned* my_dof_numbers = new unsigned[my_ndof];
        typedef std::map<unsigned long, unsigned>::iterator IT;
        unsigned pt = 0;
        for (IT it = my_dof_map.begin(); it != my_dof_map.end(); it++)
        {
          my_global_dofs[pt] = it->first;
          my_dof_numbers[pt] = it->second;
          pt++;
        }

        // and then clear the map
        my_dof_map.clear();

        // count up how many DOFs need to be sent to each processor
        int* first_dof_to_send = new int[nproc];
        int* ndof_to_send = new int[nproc];
        unsigned ptr = 0;
        for (unsigned p = 0; p < nproc; p++)
        {
          first_dof_to_send[p] = 0;
          ndof_to_send[p] = 0;
          while (ptr < my_ndof &&
                 my_global_dofs[ptr] < dense_required_rows(p, 0))
          {
            ptr++;
          }
          first_dof_to_send[p] = ptr;
          while (ptr < my_ndof &&
                 my_global_dofs[ptr] <= dense_required_rows(p, 1))
          {
            ndof_to_send[p]++;
            ptr++;
          }
        }

        // next communicate to each processor how many dofs it expects to recv
        int* ndof_to_recv = new int[nproc];
        MPI_Alltoall(ndof_to_send,
                     1,
                     MPI_INT,
                     ndof_to_recv,
                     1,
                     MPI_INT,
                     comm_pt()->mpi_comm());

        // the base displacements for the sends
        MPI_Aint base_displacement;
        MPI_Get_address(my_global_dofs, &base_displacement);

#ifdef PARANOID
        // storage for paranoid check to ensure that every row as been
        // imported
        std::vector<bool> dof_recv(nrow_local, false);
#endif

        // next send and recv
        Vector<MPI_Request> send_requests;
        Vector<MPI_Request> recv_requests;
        Vector<unsigned long*> global_dofs_recv(nproc, 0);
        Vector<unsigned*> dof_numbers_recv(nproc, 0);
        Vector<unsigned> proc;
        for (unsigned p = 0; p < nproc; p++)
        {
          if (p != my_rank)
          {
            // send
            if (ndof_to_send[p] > 0)
            {
              // the datatypes, displacements, lengths for the two datatypes
              MPI_Datatype types[2];
              MPI_Aint displacements[2];
              int lengths[2];

              // my global dofs
              MPI_Type_contiguous(
                ndof_to_send[p], MPI_UNSIGNED_LONG, &types[0]);
              MPI_Type_commit(&types[0]);
              MPI_Get_address(my_global_dofs + first_dof_to_send[p],
                              &displacements[0]);
              displacements[0] -= base_displacement;
              lengths[0] = 1;

              // my dof numbers
              MPI_Type_contiguous(ndof_to_send[p], MPI_UNSIGNED, &types[1]);
              MPI_Type_commit(&types[1]);
              MPI_Get_address(my_dof_numbers + first_dof_to_send[p],
                              &displacements[1]);
              displacements[1] -= base_displacement;
              lengths[1] = 1;

              // build the final type
              MPI_Datatype send_type;
              MPI_Type_create_struct(
                2, lengths, displacements, types, &send_type);
              MPI_Type_commit(&send_type);
              MPI_Type_free(&types[0]);
              MPI_Type_free(&types[1]);

              // and send
              MPI_Request req;
              MPI_Isend(my_global_dofs,
                        1,
                        send_type,
                        p,
                        2,
                        comm_pt()->mpi_comm(),
                        &req);
              send_requests.push_back(req);
              MPI_Type_free(&send_type);
            }

            // and recv
            if (ndof_to_recv[p] > 0)
            {
              // resize the storage
              global_dofs_recv[p] = new unsigned long[ndof_to_recv[p]];
              dof_numbers_recv[p] = new unsigned[ndof_to_recv[p]];
              proc.push_back(p);

              // the datatypes, displacements, lengths for the two datatypes
              MPI_Datatype types[2];
              MPI_Aint displacements[2];
              int lengths[2];

              // my global dofs
              MPI_Type_contiguous(
                ndof_to_recv[p], MPI_UNSIGNED_LONG, &types[0]);
              MPI_Type_commit(&types[0]);
              MPI_Get_address(global_dofs_recv[p], &displacements[0]);
              displacements[0] -= base_displacement;
              lengths[0] = 1;

              // my dof numbers
              MPI_Type_contiguous(ndof_to_recv[p], MPI_UNSIGNED, &types[1]);
              MPI_Type_commit(&types[1]);
              MPI_Get_address(dof_numbers_recv[p], &displacements[1]);
              displacements[1] -= base_displacement;
              lengths[1] = 1;

              // build the final type
              MPI_Datatype recv_type;
              MPI_Type_create_struct(
                2, lengths, displacements, types, &recv_type);
              MPI_Type_commit(&recv_type);
              MPI_Type_free(&types[0]);
              MPI_Type_free(&types[1]);

              // and recv
              MPI_Request req;
              MPI_Irecv(my_global_dofs,
                        1,
                        recv_type,
                        p,
                        2,
                        comm_pt()->mpi_comm(),
                        &req);
              recv_requests.push_back(req);
              MPI_Type_free(&recv_type);
            }
          }
          // send to self
          else
          {
            unsigned u = first_dof_to_send[p] + ndof_to_recv[p];
            for (unsigned i = first_dof_to_send[p]; i < u; i++)
            {
#ifdef PARANOID
              // indicate that this dof has ben recv
              dof_recv[my_global_dofs[i] - first_row] = true;
#endif
              Dof_number_dense[my_global_dofs[i] - first_row] =
                my_dof_numbers[i];
            }
          }
        }

        // recv and store the data
        unsigned c_recv = recv_requests.size();
        while (c_recv > 0)
        {
          // wait for any communication to finish
          int req_number;
          MPI_Waitany(
            c_recv, &recv_requests[0], &req_number, MPI_STATUS_IGNORE);
          recv_requests.erase(recv_requests.begin() + req_number);
          c_recv--;

          // determine the source processor
          unsigned p = proc[req_number];
          proc.erase(proc.begin() + req_number);

          // import the data
          for (int i = 0; i < ndof_to_recv[p]; i++)
          {
#ifdef PARANOID
            // indicate that this dof has ben recv
            dof_recv[global_dofs_recv[p][i] - first_row] = true;
#endif
            Dof_number_dense[global_dofs_recv[p][i] - first_row] =
              dof_numbers_recv[p][i];
          }

          // delete the data
          delete[] global_dofs_recv[p];
          delete[] dof_numbers_recv[p];
        }

        // finally wait for the send requests to complete as we are leaving
        // an MPI block of code
        unsigned csr = send_requests.size();
        if (csr)
        {
          MPI_Waitall(csr, &send_requests[0], MPI_STATUS_IGNORE);
        }

        // clean up
        delete[] ndof_to_send;
        delete[] first_dof_to_send;
        delete[] ndof_to_recv;
        delete[] my_global_dofs;
        delete[] my_dof_numbers;
#ifdef PARANOID
        unsigned all_recv = true;
        for (unsigned i = 0; i < nrow_local; i++)
        {
          if (!dof_recv[i])
          {
            all_recv = false;
          }
        }
        if (!all_recv)
        {
          std::ostringstream error_message;
          error_message << "Not all the DOF numbers required were received";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
#else
        std::ostringstream error_message;
        error_message
          << "The problem appears to be distributed, MPI is required";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
#endif
      }
#ifdef OOMPH_HAS_MPI
      Vector<unsigned*> sparse_rows_for_proc(nproc, 0);
      Vector<MPI_Request> sparse_rows_for_proc_requests;
      if (matrix_distributed)
      {
        // wait for number of sparse rows each processor requires
        // post recvs for that data
        if (recv_requests_sparse_nreq.size() > 0)
        {
          MPI_Waitall(recv_requests_sparse_nreq.size(),
                      &recv_requests_sparse_nreq[0],
                      MPI_STATUS_IGNORE);
        }
        for (unsigned p = 0; p < nproc; ++p)
        {
          if (nreq_sparse_for_proc[p] > 0)
          {
            MPI_Request req;
            sparse_rows_for_proc[p] = new unsigned[nreq_sparse_for_proc[p]];
            MPI_Irecv(sparse_rows_for_proc[p],
                      nreq_sparse_for_proc[p],
                      MPI_UNSIGNED,
                      p,
                      32,
                      comm_pt()->mpi_comm(),
                      &req);
            sparse_rows_for_proc_requests.push_back(req);
          }
        }
      }
#endif


      // for every global degree of freedom required by this processor we now
      // have the corresponding dof number

      // clear the Ndof_in_dof_block storage
      Dof_dimension.assign(Internal_ndof_types, 0);

      // first consider a non distributed matrix
      if (!matrix_distributed)
      {
        // set the Index_in_dof_block
        unsigned nrow = this->distribution_pt()->nrow();
        Index_in_dof_block_dense.resize(nrow);
        Index_in_dof_block_dense.initialise(0);
        for (unsigned i = 0; i < nrow; i++)
        {
          Index_in_dof_block_dense[i] = Dof_dimension[Dof_number_dense[i]];
          Dof_dimension[Dof_number_dense[i]]++;
        }
      }

      // next a distributed matrix
      else
      {
#ifdef OOMPH_HAS_MPI


        // first compute how many instances of each dof are on this
        // processor
        unsigned* my_nrows_in_dof_block = new unsigned[Internal_ndof_types];
        for (unsigned i = 0; i < Internal_ndof_types; i++)
        {
          my_nrows_in_dof_block[i] = 0;
        }
        for (unsigned i = 0; i < nrow_local; i++)
        {
          my_nrows_in_dof_block[Dof_number_dense[i]]++;
        }

        // next share the data
        unsigned* nrow_in_dof_block_recv =
          new unsigned[Internal_ndof_types * nproc];
        MPI_Allgather(my_nrows_in_dof_block,
                      Internal_ndof_types,
                      MPI_UNSIGNED,
                      nrow_in_dof_block_recv,
                      Internal_ndof_types,
                      MPI_UNSIGNED,
                      comm_pt()->mpi_comm());
        delete[] my_nrows_in_dof_block;

        // compute my first dof index and Nrows_in_dof_block
        Vector<unsigned> my_first_dof_index(Internal_ndof_types, 0);
        for (unsigned i = 0; i < Internal_ndof_types; i++)
        {
          for (unsigned p = 0; p < my_rank; p++)
          {
            my_first_dof_index[i] +=
              nrow_in_dof_block_recv[p * Internal_ndof_types + i];
          }
          Dof_dimension[i] = my_first_dof_index[i];
          for (unsigned p = my_rank; p < nproc; p++)
          {
            Dof_dimension[i] +=
              nrow_in_dof_block_recv[p * Internal_ndof_types + i];
          }
        }
        delete[] nrow_in_dof_block_recv;

        // next compute Index in dof block
        Index_in_dof_block_dense.resize(nrow_local);
        Index_in_dof_block_dense.initialise(0);
        Vector<unsigned> dof_counter(Internal_ndof_types, 0);
        for (unsigned i = 0; i < nrow_local; i++)
        {
          Index_in_dof_block_dense[i] =
            my_first_dof_index[Dof_number_dense[i]] +
            dof_counter[Dof_number_dense[i]];
          dof_counter[Dof_number_dense[i]]++;
        }

        // the base displacements for the sends
        if (sparse_rows_for_proc_requests.size() > 0)
        {
          MPI_Waitall(sparse_rows_for_proc_requests.size(),
                      &sparse_rows_for_proc_requests[0],
                      MPI_STATUS_IGNORE);
        }
        MPI_Aint base_displacement;
        MPI_Get_address(dof_number_sparse_send, &base_displacement);
        unsigned first_row = this->distribution_pt()->first_row();
        for (unsigned p = 0; p < nproc; ++p)
        {
          if (nreq_sparse_for_proc[p] > 0)
          {
            // construct the data
            index_in_dof_block_sparse_send[p] =
              new unsigned[nreq_sparse_for_proc[p]];
            dof_number_sparse_send[p] = new unsigned[nreq_sparse_for_proc[p]];
            for (unsigned i = 0; i < nreq_sparse_for_proc[p]; ++i)
            {
              unsigned r = sparse_rows_for_proc[p][i];
              r -= first_row;
              index_in_dof_block_sparse_send[p][i] =
                Index_in_dof_block_dense[r];
              dof_number_sparse_send[p][i] = Dof_number_dense[r];
            }
            delete[] sparse_rows_for_proc[p];

            // send the data
            // the datatypes, displacements, lengths for the two datatypes
            MPI_Datatype types[2];
            MPI_Aint displacements[2];
            int lengths[2];

            // index in dof block
            MPI_Type_contiguous(
              nreq_sparse_for_proc[p], MPI_UNSIGNED, &types[0]);
            MPI_Type_commit(&types[0]);
            MPI_Get_address(index_in_dof_block_sparse_send[p],
                            &displacements[0]);
            displacements[0] -= base_displacement;
            lengths[0] = 1;

            // dof number
            MPI_Type_contiguous(
              nreq_sparse_for_proc[p], MPI_UNSIGNED, &types[1]);
            MPI_Type_commit(&types[1]);
            MPI_Get_address(dof_number_sparse_send[p], &displacements[1]);
            displacements[1] -= base_displacement;
            lengths[1] = 1;

            // build the final type
            MPI_Datatype send_type;
            MPI_Type_create_struct(
              2, lengths, displacements, types, &send_type);
            MPI_Type_commit(&send_type);
            MPI_Type_free(&types[0]);
            MPI_Type_free(&types[1]);

            // and recv
            MPI_Request req;
            MPI_Isend(dof_number_sparse_send,
                      1,
                      send_type,
                      p,
                      33,
                      comm_pt()->mpi_comm(),
                      &req);
            send_requests_sparse.push_back(req);
            MPI_Type_free(&send_type);
          }
          else
          {
            index_in_dof_block_sparse_send[p] = 0;
            dof_number_sparse_send[p] = 0;
          }
        }
#endif
      }
    }

    /// //////////////////////////////////////////////////////////////////////////
    // end of master block preconditioner only operations
    /// //////////////////////////////////////////////////////////////////////////

    // compute the number of rows in each block

#ifdef PARANOID
    // check the vector is the correct length
    if (dof_to_block_map.size() != Internal_ndof_types)
    {
      std::ostringstream error_message;
      error_message << "The dof_to_block_map vector (size="
                    << dof_to_block_map.size()
                    << ") must be of size Internal_ndof_types="
                    << Internal_ndof_types;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // find the maximum block number RAYAY use std::max_element
    unsigned max_block_number = 0;
    for (unsigned i = 0; i < Internal_ndof_types; i++)
    {
      if (dof_to_block_map[i] > max_block_number)
      {
        max_block_number = dof_to_block_map[i];
      }
    }

    // resize the storage the the block to dof map
    Block_number_to_dof_number_lookup.clear();
    Block_number_to_dof_number_lookup.resize(max_block_number + 1);
    Ndof_in_block.clear();
    Ndof_in_block.resize(max_block_number + 1);

    // resize storage
    Dof_number_to_block_number_lookup.resize(Internal_ndof_types);

    // build the storage for the two maps (block to dof) and (dof to block)
    for (unsigned i = 0; i < Internal_ndof_types; i++)
    {
      Dof_number_to_block_number_lookup[i] = dof_to_block_map[i];
      Block_number_to_dof_number_lookup[dof_to_block_map[i]].push_back(i);
      Ndof_in_block[dof_to_block_map[i]]++;
    }

#ifdef PARANOID
    // paranoid check that every block number has at least one DOF associated
    // with it
    for (unsigned i = 0; i < max_block_number + 1; i++)
    {
      if (Block_number_to_dof_number_lookup[i].size() == 0)
      {
        std::ostringstream error_message;
        error_message << "block number " << i
                      << " does not have any DOFs associated with it";
        throw OomphLibWarning(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Update the number of blocks types.
    Internal_nblock_types = max_block_number + 1;

    // Distributed or not, depends on if we have more than one processor.
    bool distributed = this->master_distribution_pt()->distributed();

    // Create the new block distributions.
    Internal_block_distribution_pt.resize(Internal_nblock_types);
    for (unsigned i = 0; i < Internal_nblock_types; i++)
    {
      unsigned block_dim = 0;
      for (unsigned j = 0; j < Ndof_in_block[i]; j++)
      {
        block_dim +=
          internal_dof_block_dimension(Block_number_to_dof_number_lookup[i][j]);
      }
      Internal_block_distribution_pt[i] =
        new LinearAlgebraDistribution(comm_pt(), block_dim, distributed);
    }

    // Work out the distribution of the dof-level blocks.
    // Since several dof types may be coarsened into a single dof type.
    // We get the dof-level block distributions from the parent preconditioner.

    // How many dof types are there?
    if (is_subsidiary_block_preconditioner())
    {
      // Delete any pre-existing distributions.
      const unsigned dof_block_distribution_size =
        Dof_block_distribution_pt.size();
      for (unsigned dof_i = 0; dof_i < dof_block_distribution_size; dof_i++)
      {
        delete Dof_block_distribution_pt[dof_i];
      }
      const unsigned ndofs = this->ndof_types();
      Dof_block_distribution_pt.resize(ndofs, 0);

      // For each dof type, work out how many parent preconditioner dof types
      // are in it.
      for (unsigned dof_i = 0; dof_i < ndofs; dof_i++)
      {
        // For each external dof, we get the dofs coarsened into it (from the
        // parent preconditioner level, not the most fine grain level).
        const unsigned ncoarsened_dofs_in_dof_i =
          Doftype_coarsen_map_coarse[dof_i].size();
        Vector<LinearAlgebraDistribution*> tmp_dist_pt(ncoarsened_dofs_in_dof_i,
                                                       0);
        for (unsigned parent_dof_i = 0; parent_dof_i < ncoarsened_dofs_in_dof_i;
             parent_dof_i++)
        {
          tmp_dist_pt[parent_dof_i] =
            master_block_preconditioner_pt()->dof_block_distribution_pt(
              Doftype_in_master_preconditioner_coarse
                [Doftype_coarsen_map_coarse[dof_i][parent_dof_i]]);
        }

        Dof_block_distribution_pt[dof_i] = new LinearAlgebraDistribution;


        LinearAlgebraDistributionHelpers::concatenate(
          tmp_dist_pt, *Dof_block_distribution_pt[dof_i]);
      }
    }

    // Create Block_distribution_pt
    {
      // Delete any existing distributions in Block_distribution_pt.
      // (This should already be deleted in clear_block_preconditioner_base(...)
      // but we are just being extra safe!).
      unsigned n_existing_block_dist = Block_distribution_pt.size();
      for (unsigned dist_i = 0; dist_i < n_existing_block_dist; dist_i++)
      {
        delete Block_distribution_pt[dist_i];
      }

      Block_distribution_pt.clear();

      // Work out the distributions of the concatenated blocks.
      unsigned super_block_size = Block_to_dof_map_coarse.size();
      Block_distribution_pt.resize(super_block_size, 0);
      for (unsigned super_block_i = 0; super_block_i < super_block_size;
           super_block_i++)
      {
        unsigned sub_block_size = Block_to_dof_map_coarse[super_block_i].size();
        Vector<LinearAlgebraDistribution*> tmp_dist_pt(sub_block_size, 0);

        for (unsigned sub_block_i = 0; sub_block_i < sub_block_size;
             sub_block_i++)
        {
          tmp_dist_pt[sub_block_i] = dof_block_distribution_pt(
            Block_to_dof_map_coarse[super_block_i][sub_block_i]);
        }

        Block_distribution_pt[super_block_i] = new LinearAlgebraDistribution;

        LinearAlgebraDistributionHelpers::concatenate(
          tmp_dist_pt, *Block_distribution_pt[super_block_i]);
      }

    } // Creating Block_distribution_pt.


    // Create the distribution of the preconditioner matrix,
    // if this preconditioner is a subsidiary preconditioner then it stored
    // at Distribution_pt;
    // if this preconditioner is a master preconditioner then it is stored
    // at Internal_preconditioner_matrix_distribution_pt.
    LinearAlgebraDistribution dist;
    LinearAlgebraDistributionHelpers::concatenate(
      Internal_block_distribution_pt, dist);

    // Build the distribution.
    if (is_subsidiary_block_preconditioner())
    {
      this->build_distribution(dist);
    }
    else
    {
      Internal_preconditioner_matrix_distribution_pt =
        new LinearAlgebraDistribution(dist);
    }

    Preconditioner_matrix_distribution_pt = new LinearAlgebraDistribution;
    LinearAlgebraDistributionHelpers::concatenate(
      Block_distribution_pt, *Preconditioner_matrix_distribution_pt);

    // Clear all distributions in Auxiliary_block_distribution_pt, except for
    // the one which corresponds to the preconditioner matrix distribution. This
    // is already deleted by clear_block_preconditioner_base(...)

    // Create the key which corresponds to
    // preconditioner_matrix_distribution_pt.
    {
      const unsigned nblocks = Block_distribution_pt.size();
      Vector<unsigned> preconditioner_matrix_key(nblocks, 0);
      for (unsigned i = 0; i < nblocks; i++)
      {
        preconditioner_matrix_key[i] = i;
      }

      // Now iterate through Auxiliary_block_distribution_pt and delete
      // everything except for the value which corresponds to
      // preconditioner_matrix_key.
      std::map<Vector<unsigned>, LinearAlgebraDistribution*>::iterator iter =
        Auxiliary_block_distribution_pt.begin();
      while (iter != Auxiliary_block_distribution_pt.end())
      {
        if (iter->first != preconditioner_matrix_key)
        {
          delete iter->second;
          iter++;
        }
        else
        {
          ++iter;
        }
      }

      // Clear it just to be safe!
      Auxiliary_block_distribution_pt.clear();

      // Insert the preconditioner matrix distribution.
      insert_auxiliary_block_distribution(
        preconditioner_matrix_key, Preconditioner_matrix_distribution_pt);
    } // End of Auxiliary_block_distribution_pt encapsulation.

    // Clearing up after comm to assemble sparse lookup schemes.
#ifdef OOMPH_HAS_MPI
    if (send_requests_sparse.size() > 0)
    {
      MPI_Waitall(send_requests_sparse.size(),
                  &send_requests_sparse[0],
                  MPI_STATUS_IGNORE);
    }
    if (recv_requests_sparse.size() > 0)
    {
      MPI_Waitall(recv_requests_sparse.size(),
                  &recv_requests_sparse[0],
                  MPI_STATUS_IGNORE);
    }
    for (unsigned p = 0; p < nproc; p++)
    {
      delete[] index_in_dof_block_sparse_send[p];
      delete[] dof_number_sparse_send[p];
    }
    delete[] index_in_dof_block_sparse_send;
    delete[] dof_number_sparse_send;
    delete[] nreq_sparse;
    delete[] nreq_sparse_for_proc;
#endif

    // Next we assemble the lookup schemes for the rows
    // if the matrix is not distributed then we assemble Global_index
    // if the matrix is distributed then Rows_to_send_..., Rows_to_recv_... etc.
    if (!distributed)
    {
      // Resize the storage.
      Global_index.resize(Internal_nblock_types);
      for (unsigned b = 0; b < Internal_nblock_types; b++)
      {
        Global_index[b].resize(Internal_block_distribution_pt[b]->nrow());
      }

      // Compute:
      unsigned nrow = this->master_nrow();
      for (unsigned i = 0; i < nrow; i++)
      {
        // the dof type number;
        int dof_number = this->internal_dof_number(i);
        if (dof_number >= 0)
        {
          // the block number;
          unsigned block_number = Dof_number_to_block_number_lookup[dof_number];

          // the index in the block.
          unsigned index_in_block = 0;
          unsigned ptr = 0;
          while (int(Block_number_to_dof_number_lookup[block_number][ptr]) !=
                 dof_number)
          {
            index_in_block += internal_dof_block_dimension(
              Block_number_to_dof_number_lookup[block_number][ptr]);
            ptr++;
          }
          index_in_block += internal_index_in_dof(i);
          Global_index[block_number][index_in_block] = i;
        }
      }
    }
    // otherwise the matrix is distributed
    else
    {
#ifdef OOMPH_HAS_MPI

      // the pointer to the master distribution
      const LinearAlgebraDistribution* master_distribution_pt =
        this->master_distribution_pt();

      // resize the nrows... storage
      Nrows_to_send_for_get_block.resize(Internal_nblock_types, nproc);
      Nrows_to_send_for_get_block.initialise(0);
      Nrows_to_send_for_get_ordered.resize(nproc);
      Nrows_to_send_for_get_ordered.initialise(0);

      // loop over my rows
      unsigned nrow_local = master_distribution_pt->nrow_local();
      unsigned first_row = master_distribution_pt->first_row();
      for (unsigned i = 0; i < nrow_local; i++)
      {
        // the block number
        int b = this->internal_block_number(first_row + i);

        // check that the DOF i is associated with this preconditioner
        if (b >= 0)
        {
          // the block index
          unsigned j = this->internal_index_in_block(first_row + i);

          // the processor this row will be sent to
          unsigned block_p = 0;
          while (!(Internal_block_distribution_pt[b]->first_row(block_p) <= j &&
                   (Internal_block_distribution_pt[b]->first_row(block_p) +
                      Internal_block_distribution_pt[b]->nrow_local(block_p) >
                    j)))
          {
            block_p++;
          }

          // and increment the counter
          Nrows_to_send_for_get_block(b, block_p)++;
          Nrows_to_send_for_get_ordered[block_p]++;
        }
      }

      // resize the storage for Nrows_to_recv
      Nrows_to_recv_for_get_block.resize(Internal_nblock_types, nproc);
      Nrows_to_recv_for_get_block.initialise(0);
      Nrows_to_recv_for_get_ordered.resize(nproc);
      Nrows_to_recv_for_get_ordered.initialise(0);

      // next we send the number of rows that will be sent by this processor
      Vector<unsigned*> nrows_to_send(nproc, 0);
      Vector<unsigned*> nrows_to_recv(nproc, 0);
      Vector<MPI_Request> send_requests_nrow;
      Vector<MPI_Request> recv_requests_nrow;
      Vector<unsigned> proc;
      for (unsigned p = 0; p < nproc; p++)
      {
        if (p != my_rank)
        {
          // send
          proc.push_back(p);
          nrows_to_send[p] = new unsigned[Internal_nblock_types];
          for (unsigned b = 0; b < Internal_nblock_types; b++)
          {
            nrows_to_send[p][b] = Nrows_to_send_for_get_block(b, p);
          }
          MPI_Request s_req;
          MPI_Isend(nrows_to_send[p],
                    Internal_nblock_types,
                    MPI_UNSIGNED,
                    p,
                    3,
                    comm_pt()->mpi_comm(),
                    &s_req);
          send_requests_nrow.push_back(s_req);

          // recv
          nrows_to_recv[p] = new unsigned[Internal_nblock_types];
          MPI_Request r_req;
          MPI_Irecv(nrows_to_recv[p],
                    Internal_nblock_types,
                    MPI_UNSIGNED,
                    p,
                    3,
                    comm_pt()->mpi_comm(),
                    &r_req);
          recv_requests_nrow.push_back(r_req);
        }
        // send to self
        else
        {
          for (unsigned b = 0; b < Internal_nblock_types; b++)
          {
            Nrows_to_recv_for_get_block(b, p) =
              Nrows_to_send_for_get_block(b, p);
          }
          Nrows_to_recv_for_get_ordered[p] = Nrows_to_send_for_get_ordered[p];
        }
      }

      // create some temporary storage for the global row indices that will
      // be received from another processor.
      DenseMatrix<int*> block_rows_to_send(Internal_nblock_types, nproc, 0);
      Vector<int*> ordered_rows_to_send(nproc, 0);

      // resize the rows... storage
      Rows_to_send_for_get_block.resize(Internal_nblock_types, nproc);
      Rows_to_send_for_get_block.initialise(0);
      Rows_to_send_for_get_ordered.resize(nproc);
      Rows_to_send_for_get_ordered.initialise(0);
      Rows_to_recv_for_get_block.resize(Internal_nblock_types, nproc);
      Rows_to_recv_for_get_block.initialise(0);

      // resize the storage
      for (unsigned p = 0; p < nproc; p++)
      {
        for (unsigned b = 0; b < Internal_nblock_types; b++)
        {
          Rows_to_send_for_get_block(b, p) =
            new int[Nrows_to_send_for_get_block(b, p)];
          if (p != my_rank)
          {
            block_rows_to_send(b, p) =
              new int[Nrows_to_send_for_get_block(b, p)];
          }
          else
          {
            Rows_to_recv_for_get_block(b, p) =
              new int[Nrows_to_send_for_get_block(b, p)];
          }
        }
        Rows_to_send_for_get_ordered[p] =
          new int[Nrows_to_send_for_get_ordered[p]];
      }


      // loop over my rows to allocate the nrows
      DenseMatrix<unsigned> ptr_block(Internal_nblock_types, nproc, 0);
      for (unsigned i = 0; i < nrow_local; i++)
      {
        // the block number
        int b = this->internal_block_number(first_row + i);

        // check that the DOF i is associated with this preconditioner
        if (b >= 0)
        {
          // the block index
          unsigned j = this->internal_index_in_block(first_row + i);

          // the processor this row will be sent to
          unsigned block_p = 0;
          while (!(Internal_block_distribution_pt[b]->first_row(block_p) <= j &&
                   (Internal_block_distribution_pt[b]->first_row(block_p) +
                      Internal_block_distribution_pt[b]->nrow_local(block_p) >
                    j)))
          {
            block_p++;
          }

          // and store the row
          Rows_to_send_for_get_block(b, block_p)[ptr_block(b, block_p)] = i;
          if (block_p != my_rank)
          {
            block_rows_to_send(b, block_p)[ptr_block(b, block_p)] =
              j - Internal_block_distribution_pt[b]->first_row(block_p);
          }
          else
          {
            Rows_to_recv_for_get_block(b, block_p)[ptr_block(b, block_p)] =
              j - Internal_block_distribution_pt[b]->first_row(block_p);
          }
          ptr_block(b, block_p)++;
        }
      }

      // next block ordered
      for (unsigned p = 0; p < nproc; ++p)
      {
        int pt = 0;
        for (unsigned b = 0; b < Internal_nblock_types; ++b)
        {
          for (unsigned i = 0; i < Nrows_to_send_for_get_block(b, p); ++i)
          {
            Rows_to_send_for_get_ordered[p][pt] =
              Rows_to_send_for_get_block(b, p)[i];
            pt++;
          }
        }
      }

      // next process the nrow recvs as they complete

      // recv and store the data
      unsigned c = recv_requests_nrow.size();
      while (c > 0)
      {
        // wait for any communication to finish
        int req_number;
        MPI_Waitany(c, &recv_requests_nrow[0], &req_number, MPI_STATUS_IGNORE);
        recv_requests_nrow.erase(recv_requests_nrow.begin() + req_number);
        c--;

        // determine the source processor
        unsigned p = proc[req_number];
        proc.erase(proc.begin() + req_number);

        // copy the data to its final storage
        Nrows_to_recv_for_get_ordered[p] = 0;
        for (unsigned b = 0; b < Internal_nblock_types; b++)
        {
          Nrows_to_recv_for_get_block(b, p) = nrows_to_recv[p][b];
          Nrows_to_recv_for_get_ordered[p] += nrows_to_recv[p][b];
        }

        // and clear
        delete[] nrows_to_recv[p];
      }

      // resize the storage for the incoming rows data
      Rows_to_recv_for_get_ordered.resize(nproc, 0);
      for (unsigned p = 0; p < nproc; p++)
      {
        if (p != my_rank)
        {
          for (unsigned b = 0; b < Internal_nblock_types; b++)
          {
            Rows_to_recv_for_get_block(b, p) =
              new int[Nrows_to_recv_for_get_block(b, p)];
          }
        }
      }

      // compute the number of sends and recv from this processor
      // to each other processor
      Vector<unsigned> nsend_for_rows(nproc, 0);
      Vector<unsigned> nrecv_for_rows(nproc, 0);
      for (unsigned p = 0; p < nproc; p++)
      {
        if (p != my_rank)
        {
          for (unsigned b = 0; b < Internal_nblock_types; b++)
          {
            if (Nrows_to_send_for_get_block(b, p) > 0)
            {
              nsend_for_rows[p]++;
            }
            if (Nrows_to_recv_for_get_block(b, p) > 0)
            {
              nrecv_for_rows[p]++;
            }
          }
        }
      }

      // finally post the sends and recvs
      MPI_Aint base_displacement;
      MPI_Get_address(matrix_pt(), &base_displacement);
      Vector<MPI_Request> req_rows;
      for (unsigned p = 0; p < nproc; p++)
      {
        if (p != my_rank)
        {
          // send
          if (nsend_for_rows[p] > 0)
          {
            MPI_Datatype send_types[nsend_for_rows[p]];
            MPI_Aint send_displacements[nsend_for_rows[p]];
            int send_sz[nsend_for_rows[p]];
            unsigned send_ptr = 0;
            for (unsigned b = 0; b < Internal_nblock_types; b++)
            {
              if (Nrows_to_send_for_get_block(b, p) > 0)
              {
                MPI_Type_contiguous(Nrows_to_send_for_get_block(b, p),
                                    MPI_INT,
                                    &send_types[send_ptr]);
                MPI_Type_commit(&send_types[send_ptr]);
                MPI_Get_address(block_rows_to_send(b, p),
                                &send_displacements[send_ptr]);
                send_displacements[send_ptr] -= base_displacement;
                send_sz[send_ptr] = 1;
                send_ptr++;
              }
            }
            MPI_Datatype final_send_type;
            MPI_Type_create_struct(nsend_for_rows[p],
                                   send_sz,
                                   send_displacements,
                                   send_types,
                                   &final_send_type);
            MPI_Type_commit(&final_send_type);
            for (unsigned i = 0; i < nsend_for_rows[p]; i++)
            {
              MPI_Type_free(&send_types[i]);
            }
            MPI_Request send_req;
            MPI_Isend(matrix_pt(),
                      1,
                      final_send_type,
                      p,
                      4,
                      comm_pt()->mpi_comm(),
                      &send_req);
            req_rows.push_back(send_req);
            MPI_Type_free(&final_send_type);
          }

          // recv
          if (nrecv_for_rows[p] > 0)
          {
            MPI_Datatype recv_types[nrecv_for_rows[p]];
            MPI_Aint recv_displacements[nrecv_for_rows[p]];
            int recv_sz[nrecv_for_rows[p]];
            unsigned recv_ptr = 0;
            for (unsigned b = 0; b < Internal_nblock_types; b++)
            {
              if (Nrows_to_recv_for_get_block(b, p) > 0)
              {
                MPI_Type_contiguous(Nrows_to_recv_for_get_block(b, p),
                                    MPI_INT,
                                    &recv_types[recv_ptr]);
                MPI_Type_commit(&recv_types[recv_ptr]);
                MPI_Get_address(Rows_to_recv_for_get_block(b, p),
                                &recv_displacements[recv_ptr]);
                recv_displacements[recv_ptr] -= base_displacement;
                recv_sz[recv_ptr] = 1;
                recv_ptr++;
              }
            }
            MPI_Datatype final_recv_type;
            MPI_Type_create_struct(nrecv_for_rows[p],
                                   recv_sz,
                                   recv_displacements,
                                   recv_types,
                                   &final_recv_type);
            MPI_Type_commit(&final_recv_type);
            for (unsigned i = 0; i < nrecv_for_rows[p]; i++)
            {
              MPI_Type_free(&recv_types[i]);
            }
            MPI_Request recv_req;
            MPI_Irecv(matrix_pt(),
                      1,
                      final_recv_type,
                      p,
                      4,
                      comm_pt()->mpi_comm(),
                      &recv_req);
            req_rows.push_back(recv_req);
            MPI_Type_free(&final_recv_type);
          }
        }
      }

      // cleaning up Waitalls


      // wait for the recv requests so we can compute
      // Nrows_to_recv_for_get_ordered
      unsigned n_req_rows = req_rows.size();
      if (n_req_rows)
      {
        MPI_Waitall(n_req_rows, &req_rows[0], MPI_STATUS_IGNORE);
      }

      // resize the storage
      Rows_to_recv_for_get_ordered.resize(nproc);
      Rows_to_recv_for_get_ordered.initialise(0);

      // construct block offset
      Vector<int> vec_offset(Internal_nblock_types, 0);
      for (unsigned b = 1; b < Internal_nblock_types; ++b)
      {
        vec_offset[b] = vec_offset[b - 1] +
                        Internal_block_distribution_pt[b - 1]->nrow_local();
      }

      //
      for (unsigned p = 0; p < nproc; p++)
      {
        int pt = 0;
        Rows_to_recv_for_get_ordered[p] =
          new int[Nrows_to_recv_for_get_ordered[p]];
        for (unsigned b = 0; b < Internal_nblock_types; b++)
        {
          for (unsigned i = 0; i < Nrows_to_recv_for_get_block(b, p); i++)
          {
            Rows_to_recv_for_get_ordered[p][pt] =
              Rows_to_recv_for_get_block(b, p)[i] + vec_offset[b];
            pt++;
          }
        }
      }

      // clean up
      for (unsigned p = 0; p < nproc; p++)
      {
        if (p != my_rank)
        {
          for (unsigned b = 0; b < Internal_nblock_types; b++)
          {
            delete[] block_rows_to_send(b, p);
          }
          if (Nrows_to_send_for_get_ordered[p] > 0)
          {
            delete[] ordered_rows_to_send[p];
          }
        }
      }

      // and the send reqs
      unsigned n_req_send_nrow = send_requests_nrow.size();
      if (n_req_send_nrow)
      {
        MPI_Waitall(n_req_send_nrow, &send_requests_nrow[0], MPI_STATUS_IGNORE);
      }
      for (unsigned p = 0; p < nproc; p++)
      {
        delete[] nrows_to_send[p];
      }
#endif
    }

    // If we asked for output of blocks to a file then do it.
    if (block_output_on()) output_blocks_to_files(Output_base_filename);
  }

  //============================================================================
  //??ds
  /// Function to turn this preconditioner into a
  /// subsidiary preconditioner that operates within a bigger
  /// "master block preconditioner (e.g. a Navier-Stokes 2x2 block
  /// preconditioner dealing with the fluid sub-blocks within a
  /// 3x3 FSI preconditioner. Once this is done the master block
  /// preconditioner deals with the block setup etc.
  /// The vector block_map must specify the dof number in the
  /// master preconditioner that corresponds to a block number in this
  /// preconditioner. ??ds horribly misleading comment!
  /// The length of the vector is used to determine the number of
  /// blocks in this preconditioner therefore it must be correctly sized.
  /// This calls the other turn_into_subsidiary_block_preconditioner(...)
  /// function providing an empty doftype_to_doftype_map vector.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::turn_into_subsidiary_block_preconditioner(
    BlockPreconditioner<MATRIX>* master_block_prec_pt,
    const Vector<unsigned>& doftype_in_master_preconditioner_coarse)
  {
    // Create the identity dof_coarsen_map
    Vector<Vector<unsigned>> doftype_coarsen_map_coarse;
    unsigned doftype_in_master_preconditioner_coarse_size =
      doftype_in_master_preconditioner_coarse.size();

    for (unsigned dof_i = 0;
         dof_i < doftype_in_master_preconditioner_coarse_size;
         dof_i++)
    {
      // Create a vector of size 1 and value i,
      // then push it into the dof_coarsen_map vector.
      Vector<unsigned> tmp_vec(1, dof_i);
      doftype_coarsen_map_coarse.push_back(tmp_vec);
    }

    // Call the other turn_into_subsidiary_block_preconditioner function.
    turn_into_subsidiary_block_preconditioner(
      master_block_prec_pt,
      doftype_in_master_preconditioner_coarse,
      doftype_coarsen_map_coarse);
  }


  //============================================================================
  /// Function to turn this block preconditioner into a
  /// subsidiary block preconditioner that operates within a bigger
  /// master block preconditioner (e.g. a Navier-Stokes 2x2 block
  /// preconditioner dealing with the fluid sub-blocks within a
  /// 3x3 FSI preconditioner. Once this is done the master block
  /// preconditioner deals with the block setup etc.
  ///
  /// The vector doftype_map must specify the dof type in the
  /// master preconditioner that corresponds to a dof type in this block
  /// preconditioner.
  ///
  /// In general, we want:
  /// doftype_map[doftype in subsidiary prec] = doftype in master prec.
  ///
  /// It tells this block preconditioner which dof types of the master
  /// block preconditioner it is working with.
  ///
  /// The length of the vector is used to determine the number of
  /// dof types in THIS block preconditioner therefore it must be correctly
  /// sized.
  ///
  /// For example, let the master block preconditioner have 5 dof types in total
  /// and a 1-4 dof type splitting where the block (0,0) corresponds to
  /// dof type 0 and the block (1,1) corresponds to dof types 1, 2, 3 and 4
  /// (i.e. it would have given to block_setup the vector [0,1,1,1,1]).
  /// Furthermore, it solves (1,1) block with subsidiary block preconditioner.
  /// Then the doftype_map passed to this function of the subsidiary block
  /// preconditioner would be [1, 2, 3, 4].
  ///
  /// Dof type coarsening (following on from the example above):
  /// Let the subsidiary block preconditioner (THIS block preconditioner)
  /// only works with two DOF types, then the master block preconditioner must
  /// "coarsen" the dof types by providing the optional argument
  /// doftype_coarsen_map vector.
  ///
  /// The doftype_coarsen_map vector (in this case) might be [[0,1], [2,3]]
  /// telling the subsidiary block preconditioner that the SUBSIDIARY dof types
  /// 0 and 1 should be treated as dof type 0 and the subsidiary dof types 2
  /// and 3 should be treated as subsidiary dof type 1.
  ///
  /// If no doftype_coarsen_map vector is provided, then the identity is
  /// used automatically (see the turn_into_subsidiary_block_preconditioner(...)
  /// function with only two arguments). In the above case, the identity
  /// doftype_coarsen_map vector for the subsidiary block preconditioner
  /// would be the 2D vector [[0], [1], [2], [3]] which means
  /// dof type 0 is treated as dof type 0,
  /// dof type 1 is treated as dof type 1,
  /// dof type 2 is treated as dof type 2, and
  /// dof type 3 is treated as dof type 3.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::turn_into_subsidiary_block_preconditioner(
    BlockPreconditioner<MATRIX>* master_block_prec_pt,
    const Vector<unsigned>& doftype_in_master_preconditioner_coarse,
    const Vector<Vector<unsigned>>& doftype_coarsen_map_coarse)
  {
    // Set the master block preconditioner pointer
    Master_block_preconditioner_pt = master_block_prec_pt;

    // Set the Doftype_coarsen_map_coarse.
    Doftype_coarsen_map_coarse = doftype_coarsen_map_coarse;

    Doftype_in_master_preconditioner_coarse =
      doftype_in_master_preconditioner_coarse;
  } // end of turn_into_subsidiary_block_preconditioner(...)


  //============================================================================
  /// Determine the size of the matrix blocks and setup the
  /// lookup schemes relating the global degrees of freedom with
  /// their "blocks" and their indices (row/column numbers) in those
  /// blocks.
  /// The distributions of the preconditioner and the blocks are
  /// automatically specified (and assumed to be uniform) at this
  /// stage.
  /// This method should be used if each DOF type corresponds to a
  /// unique block type.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::block_setup()
  {
#ifdef PARANOID

    // Subsidiary preconditioners don't really need the meshes
    if (this->is_master_block_preconditioner())
    {
      std::ostringstream err_msg;
      unsigned n = nmesh();
      if (n == 0)
      {
        err_msg << "No meshes have been set for this block preconditioner!\n"
                << "Set one with set_nmesh(...), set_mesh(...)" << std::endl;
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        for (unsigned m = 0; m < n; m++)
        {
          if (Mesh_pt[m] == 0)
          {
            err_msg << "The mesh pointer to mesh " << m << " is null!\n"
                    << "Set a non-null one with set_mesh(...)" << std::endl;
            throw OomphLibError(
              err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
          }
        }
      }
    }
#endif

    // Get the number of dof types.
    unsigned internal_n_dof_types = ndof_types();

    // Build the dof to block map - assume that each type of dof corresponds
    // to a different type of block.
    Vector<unsigned> dof_to_block_lookup(internal_n_dof_types);
    for (unsigned i = 0; i < internal_n_dof_types; i++)
    {
      dof_to_block_lookup[i] = i;
    }

    // call the block setup method
    this->block_setup(dof_to_block_lookup);
  }


  //============================================================================
  /// Get the block matrices required for the block preconditioner. Takes a
  /// pointer to a matrix of bools that indicate if a specified sub-block is
  /// required for the preconditioning operation. Computes the required block
  /// matrices, and stores pointers to them in the matrix block_matrix_pt. If an
  /// entry in block_matrix_pt is equal to NULL that sub-block has not been
  /// requested and is therefore not available.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::get_blocks(
    DenseMatrix<bool>& required_blocks,
    DenseMatrix<MATRIX*>& block_matrix_pt) const
  {
    // Cache number of block types
    const unsigned n_block_types = nblock_types();

#ifdef PARANOID
    // If required blocks matrix pointer is not the correct size then abort.
    if ((required_blocks.nrow() != n_block_types) ||
        (required_blocks.ncol() != n_block_types))
    {
      std::ostringstream error_message;
      error_message << "The size of the matrix of bools required_blocks "
                    << "(which indicates which blocks are required) is not the "
                    << "right size, required_blocks is "
                    << required_blocks.ncol() << " x " << required_blocks.nrow()
                    << ", whereas it should "
                    << "be " << n_block_types << " x " << n_block_types;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // If block matrix pointer is not the correct size then abort.
    if ((block_matrix_pt.nrow() != n_block_types) ||
        (block_matrix_pt.ncol() != n_block_types))
    {
      std::ostringstream error_message;
      error_message << "The size of the block matrix pt is not the "
                    << "right size, block_matrix_pt is "
                    << block_matrix_pt.ncol() << " x " << block_matrix_pt.nrow()
                    << ", whereas it should "
                    << "be " << n_block_types << " x " << n_block_types;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

#endif

    // Loop over the blocks
    for (unsigned i = 0; i < n_block_types; i++)
    {
      for (unsigned j = 0; j < n_block_types; j++)
      {
        // If block(i,j) is required then create a matrix and fill it in.
        if (required_blocks(i, j))
        {
          //??ds might want to remove this use of new as well?
          block_matrix_pt(i, j) = new MATRIX;
          get_block(i, j, *block_matrix_pt(i, j));
        }

        // Otherwise set pointer to null.
        else
        {
          block_matrix_pt(i, j) = 0;
        }
      }
    }
  }

  //============================================================================
  /// Takes the naturally ordered vector and extracts the blocks
  /// indicated by the block number (the values) in the Vector
  /// block_vec_number all at once, then concatenates them without
  /// communication. Here, the values in block_vec_number is the block number
  /// in the current preconditioner.
  /// This is a non-const function because distributions may be created
  /// and stored in Auxiliary_block_distribution_pt for future use.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::get_concatenated_block_vector(
    const Vector<unsigned>& block_vec_number,
    const DoubleVector& v,
    DoubleVector& w)
  {
#ifdef PARANOID

    // Check if v is built.
    if (!v.built())
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // v must have the same distribution as the upper-most master block
    // preconditioner, since the upper-most master block preconditioner
    // should have the same distribution as the matrix pointed to
    // by matrix_pt().
    if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the global vector v must match the "
              << " specified master_distribution_pt(). \n"
              << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check to see if there are more blocks defined in the block_vec_number
    // vector than the number of block types. This is not allowed.
    const unsigned para_nblock_types = nblock_types();
    const unsigned para_block_vec_number_size = block_vec_number.size();
    if (para_block_vec_number_size > para_nblock_types)
    {
      std::ostringstream err_msg;
      err_msg << "You have requested " << para_block_vec_number_size
              << " number of blocks, (block_vec_number.size() is "
              << para_block_vec_number_size << ").\n"
              << "But there are only " << para_nblock_types
              << " nblock_types.\n"
              << "Please make sure that block_vec_number is correctly sized.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check if any block numbers defined in block_vec_number is equal to or
    // greater than the number of block types.
    // E.g. if there are 5 block types, we can only have block numbers:
    //  0, 1, 2, 3 and 4.
    for (unsigned i = 0; i < para_block_vec_number_size; i++)
    {
      const unsigned para_required_block = block_vec_number[i];
      if (para_required_block >= para_nblock_types)
      {
        std::ostringstream err_msg;
        err_msg << "block_vec_number[" << i << "] is " << para_required_block
                << ".\n"
                << "But there are only " << para_nblock_types
                << " nblock_types.\n";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }

    // Check that no block number is inserted twice.
    std::set<unsigned> para_set;
    for (unsigned b = 0; b < para_block_vec_number_size; b++)
    {
      std::pair<std::set<unsigned>::iterator, bool> para_set_ret;
      para_set_ret = para_set.insert(block_vec_number[b]);

      if (!para_set_ret.second)
      {
        std::ostringstream err_msg;
        err_msg << "Error: the block number " << block_vec_number[b]
                << " appears twice.\n";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Number of blocks to get.
    const unsigned n_block = block_vec_number.size();

    // Each block is made of dof types. We get the most fine grain dof types.
    // Most fine grain in the sense that these are the dof types that belongs
    // in this block before any coarsening of dof types has taken place.
    // The ordering of the dof types matters, this is handled properly when
    // creating the Block_to_dof_map_fine vector and must be respected here.
    // I.e. we cannot arbitrarily insert dof types (even if they are correct)
    // in the vector most_fine_grain_dof.
    Vector<unsigned> most_fine_grain_dof;
    for (unsigned b = 0; b < n_block; b++)
    {
      const unsigned mapped_b = block_vec_number[b];
      most_fine_grain_dof.insert(most_fine_grain_dof.end(),
                                 Block_to_dof_map_fine[mapped_b].begin(),
                                 Block_to_dof_map_fine[mapped_b].end());
    }

    // Get all the dof level vectors in one go.
    Vector<DoubleVector> dof_block_vector;
    internal_get_block_vectors(most_fine_grain_dof, v, dof_block_vector);

    // Next we need to build the output DoubleVector w with the correct
    // distribution: the concatenation of the distributions of all the
    // dof-level vectors. This is the same as the concatenation of the
    // distributions of the blocks within this preconditioner.
    //
    // So we first check if it exists already, if not, we create it and
    // store it for future use. We store it because concatenation of
    // distributions requires communication, so concatenation of
    // distributions on-the-fly should be avoided.
    std::map<Vector<unsigned>, LinearAlgebraDistribution*>::const_iterator iter;

    // Attempt to get an iterator pointing to the pair with the value
    // block_vec_number.
    iter = Auxiliary_block_distribution_pt.find(block_vec_number);

    if (iter != Auxiliary_block_distribution_pt.end())
    // If it exists, build w with the distribution pointed to
    // by pair::second.
    {
      w.build(iter->second);
    }
    else
    // Else, we need to create the distribution and store it in
    // Auxiliary_block_distribution_pt.
    {
      Vector<LinearAlgebraDistribution*> tmp_vec_dist_pt(n_block, 0);
      for (unsigned b = 0; b < n_block; b++)
      {
        tmp_vec_dist_pt[b] = Block_distribution_pt[block_vec_number[b]];
      }

      // Note that the distribution is created with new but not deleted here.
      // This is handled in the clean up functions.
      LinearAlgebraDistribution* tmp_dist_pt = new LinearAlgebraDistribution;
      LinearAlgebraDistributionHelpers::concatenate(tmp_vec_dist_pt,
                                                    *tmp_dist_pt);

      // Store the pair of Vector<unsigned> and LinearAlgebraDistribution*
      insert_auxiliary_block_distribution(block_vec_number, tmp_dist_pt);

      // Build w.
      w.build(tmp_dist_pt);
    }

    // Now concatenate all the dof level vectors into the vector w.
    DoubleVectorHelpers::concatenate_without_communication(dof_block_vector, w);

  } // get_concatenated_block_vector(...)

  //============================================================================
  /// Takes concatenated block ordered vector, b, and copies its
  // entries to the appropriate entries in the naturally ordered vector, v.
  // Here the values in block_vec_number indicates which blocks the vector
  // b is a concatenation of. The block number are those in the current
  // preconditioner. If the preconditioner is a subsidiary block
  // preconditioner the other entries in v that are not associated with it
  // are left alone.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::return_concatenated_block_vector(
    const Vector<unsigned>& block_vec_number,
    const DoubleVector& w,
    DoubleVector& v) const
  {
#ifdef PARANOID

    // Check if v is built.
    if (!v.built())
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // v must have the same distribution as the upper-most master block
    // preconditioner, since the upper-most master block preconditioner
    // should have the same distribution as the matrix pointed to
    // by matrix_pt().
    if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the global vector v must match the "
              << " specified master_distribution_pt(). \n"
              << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check to see if there are more blocks defined in the block_vec_number
    // vector than the number of block types. This is not allowed.
    const unsigned para_block_vec_number_size = block_vec_number.size();
    const unsigned para_n_block = nblock_types();
    if (para_block_vec_number_size > para_n_block)
    {
      std::ostringstream err_msg;
      err_msg << "Trying to return " << para_block_vec_number_size
              << " block vectors.\n"
              << "But there are only " << para_n_block << " block types.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check if any block numbers defined in block_vec_number is equal to or
    // greater than the number of block types.
    // E.g. if there are 5 block types, we can only have block numbers:
    //  0, 1, 2, 3 and 4.
    for (unsigned b = 0; b < para_block_vec_number_size; b++)
    {
      const unsigned para_required_block = block_vec_number[b];
      if (para_required_block > para_n_block)
      {
        std::ostringstream err_msg;
        err_msg << "block_vec_number[" << b << "] is " << para_required_block
                << ".\n"
                << "But there are only " << para_n_block << " block types.\n";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }

    // Check that no block number is inserted twice.
    std::set<unsigned> para_set;
    for (unsigned b = 0; b < para_block_vec_number_size; b++)
    {
      std::pair<std::set<unsigned>::iterator, bool> para_set_ret;
      para_set_ret = para_set.insert(block_vec_number[b]);

      if (!para_set_ret.second)
      {
        std::ostringstream err_msg;
        err_msg << "Error: the block number " << block_vec_number[b]
                << " appears twice.\n";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }

    // Check that w is built.
    if (!w.built())
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the block vector w must be setup.";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check that the distributions defined by block_vec_number is correct for
    // the distribution from w.
    // Recall that w is the concatenation of the block vectors defined by
    // the values in block_vec_number. We check that this is the case.
    Vector<LinearAlgebraDistribution*> para_vec_dist_pt(
      para_block_vec_number_size, 0);

    for (unsigned b = 0; b < para_block_vec_number_size; b++)
    {
      para_vec_dist_pt[b] = Block_distribution_pt[block_vec_number[b]];
    }

    LinearAlgebraDistribution para_tmp_dist;

    LinearAlgebraDistributionHelpers::concatenate(para_vec_dist_pt,
                                                  para_tmp_dist);

    if (*w.distribution_pt() != para_tmp_dist)
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the block vector w does not match \n"
              << "the concatenation of the block distributions defined in \n"
              << "block_vec_number.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Number of blocks to return.
    const unsigned n_block = block_vec_number.size();

    // Each block is made of dof types. We get the most fine grain dof types.
    // Most fine grain in the sense that these are the dof types that belongs
    // in this block before any coarsening of dof types has taken place.
    // The ordering of the dof types matters, this is handled properly when
    // creating the Block_to_dof_map_fine vector and must be respected here.
    // I.e. we cannot arbitrarily insert dof types (even if they are correct)
    // in the vector most_fine_grain_dof.
    Vector<unsigned> most_fine_grain_dof;
    for (unsigned b = 0; b < n_block; b++)
    {
      const unsigned mapped_b = block_vec_number[b];
      most_fine_grain_dof.insert(most_fine_grain_dof.end(),
                                 Block_to_dof_map_fine[mapped_b].begin(),
                                 Block_to_dof_map_fine[mapped_b].end());
    }

    // The number of most fine grain dof types associated with the blocks
    // defined by block_vec_number.
    const unsigned ndof = most_fine_grain_dof.size();

    // Build each dof level vector with the correct distribution.
    Vector<DoubleVector> dof_vector(ndof);
    for (unsigned d = 0; d < ndof; d++)
    {
      dof_vector[d].build(
        internal_block_distribution_pt(most_fine_grain_dof[d]));
    }

    // Perform the splitting of w into the most fine grain dof level vectors.
    DoubleVectorHelpers::split_without_communication(w, dof_vector);

    // Return all the dof level vectors in one go.
    internal_return_block_vectors(most_fine_grain_dof, dof_vector, v);
  } // return_concatenated_block_vector(...)

  //============================================================================
  /// Takes the naturally ordered vector and rearranges it into a
  /// vector of sub vectors corresponding to the blocks, so s[b][i] contains
  /// the i-th entry in the vector associated with block b.
  /// Note: If the preconditioner is a subsidiary preconditioner then only the
  /// sub-vectors associated with the blocks of the subsidiary preconditioner
  /// will be included. Hence the length of v is master_nrow() whereas the
  /// total length of the s vectors is the sum of the lengths of the
  /// individual block vectors defined in block_vec_number.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::get_block_vectors(
    const Vector<unsigned>& block_vec_number,
    const DoubleVector& v,
    Vector<DoubleVector>& s) const
  {
#ifdef PARANOID

    // Check if v is built.
    if (!v.built())
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // v must have the same distribution as the upper-most master block
    // preconditioner, since the upper-most master block preconditioner
    // should have the same distribution as the matrix pointed to
    // by matrix_pt().
    if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the global vector v must match the "
              << " specified master_distribution_pt(). \n"
              << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check to see if there are more blocks defined in the block_vec_number
    // vector than the number of block types. This is not allowed.
    const unsigned para_nblock_types = nblock_types();
    const unsigned para_block_vec_number_size = block_vec_number.size();
    if (para_block_vec_number_size > para_nblock_types)
    {
      std::ostringstream err_msg;
      err_msg << "You have requested " << para_block_vec_number_size
              << " number of blocks, (block_vec_number.size() is "
              << para_block_vec_number_size << ").\n"
              << "But there are only " << para_nblock_types
              << " nblock_types.\n"
              << "Please make sure that block_vec_number is correctly sized.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check if any block numbers defined in block_vec_number is equal to or
    // greater than the number of block types.
    // E.g. if there are 5 block types, we can only have block numbers:
    //  0, 1, 2, 3 and 4.
    for (unsigned i = 0; i < para_block_vec_number_size; i++)
    {
      const unsigned para_required_block = block_vec_number[i];
      if (para_required_block > para_nblock_types)
      {
        std::ostringstream err_msg;
        err_msg << "block_vec_number[" << i << "] is " << para_required_block
                << ".\n"
                << "But there are only " << para_nblock_types
                << " nblock_types.\n";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
    // Check that no block number is inserted twice.
    std::set<unsigned> para_set;
    for (unsigned b = 0; b < para_block_vec_number_size; b++)
    {
      std::pair<std::set<unsigned>::iterator, bool> para_set_ret;
      para_set_ret = para_set.insert(block_vec_number[b]);

      if (!para_set_ret.second)
      {
        std::ostringstream err_msg;
        err_msg << "Error: the block number " << block_vec_number[b]
                << " appears twice.\n";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Number of blocks to get.
    const unsigned n_block = block_vec_number.size();
    s.resize(n_block);

    // Each block is made of dof types. We get the most fine grain dof types.
    // Most fine grain in the sense that these are the dof types that belongs
    // in this block before any coarsening of dof types has taken place.
    // The ordering of the dof types matters, this is handled properly when
    // creating the Block_to_dof_map_fine vector and must be respected here.
    // I.e. we cannot arbitrarily insert dof types (even if they are correct)
    // in the vector most_fine_grain_dof.
    Vector<unsigned> most_fine_grain_dof;
    for (unsigned b = 0; b < n_block; b++)
    {
      const unsigned mapped_b = block_vec_number[b];

      most_fine_grain_dof.insert(most_fine_grain_dof.end(),
                                 Block_to_dof_map_fine[mapped_b].begin(),
                                 Block_to_dof_map_fine[mapped_b].end());
    }

    // Get all the dof level vectors in one go.
    Vector<DoubleVector> dof_vector;
    internal_get_block_vectors(most_fine_grain_dof, v, dof_vector);

    // For each block vector requested,
    // build the block s[b],
    // concatenate the corresponding dof vector

    // Since all the dof vectors are in dof_vector,
    // we need to loop through this.
    // The offset helps us loop through this.
    unsigned offset = 0;

    for (unsigned b = 0; b < n_block; b++)
    {
      // The actual block number required.
      const unsigned mapped_b = block_vec_number[b];

      // How many most fine grain dofs are in this block?
      const unsigned n_dof = Block_to_dof_map_fine[mapped_b].size();

      if (n_dof == 1)
      // No need to concatenate, just copy the DoubleVector.
      {
        s[b] = dof_vector[offset];
      }
      else
      // Concatenate the relevant dof vectors into s[b].
      {
        s[b].build(Block_distribution_pt[mapped_b], 0);
        Vector<DoubleVector*> tmp_vec_pt(n_dof, 0);
        for (unsigned vec_i = 0; vec_i < n_dof; vec_i++)
        {
          tmp_vec_pt[vec_i] = &dof_vector[offset + vec_i];
        }

        DoubleVectorHelpers::concatenate_without_communication(tmp_vec_pt,
                                                               s[b]);
      }

      // Update the offset.
      offset += n_dof;
    }
  } // get_block_vectors(...)


  //============================================================================
  /// Takes the naturally ordered vector and rearranges it into a
  /// vector of sub vectors corresponding to the blocks, so s[b][i] contains
  /// the i-th entry in the vector associated with block b.
  /// Note: If the preconditioner is a subsidiary preconditioner then only the
  /// sub-vectors associated with the blocks of the subsidiary preconditioner
  /// will be included. Hence the length of v is master_nrow() whereas the
  /// total length of the s vectors is Nrow.
  /// This is simply a wrapper around the other get_block_vectors(...) function
  /// where the block_vec_number Vector is the identity, i.e.
  /// block_vec_number is [0, 1, ..., nblock_types - 1].
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::get_block_vectors(
    const DoubleVector& v, Vector<DoubleVector>& s) const
  {
    // Get the number of blocks in this block preconditioner.
    const unsigned n_block = nblock_types();

    // Create the identity vector.
    Vector<unsigned> required_block(n_block, 0);
    for (unsigned i = 0; i < n_block; i++)
    {
      required_block[i] = i;
    }

    // Call the other function which does the work.
    get_block_vectors(required_block, v, s);
  }

  //============================================================================
  /// Takes the naturally ordered vector and
  /// rearranges it into a vector of sub vectors corresponding to the blocks,
  /// so s[b][i] contains the i-th entry in the vector associated with block b.
  /// The block_vec_number indicates which blocks we want.
  /// These blocks and vectors are those corresponding to the internal blocks.
  /// Note: If the preconditioner is a subsidiary preconditioner then only the
  /// sub-vectors associated with the blocks of the subsidiary preconditioner
  /// will be included. Hence the length of v is master_nrow() whereas the
  /// total length of the s vectors is the sum of the Nrow of the sub vectors.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::internal_get_block_vectors(
    const Vector<unsigned>& block_vec_number,
    const DoubleVector& v,
    Vector<DoubleVector>& s) const
  {
#ifdef PARANOID
    if (!v.built())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must match the "
                    << " specified master_distribution_pt(). \n"
                    << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Number of block types
    // const unsigned nblock = this->internal_nblock_types();
    const unsigned nblock = block_vec_number.size();

    // if + only one processor
    //    + more than one processor but matrix_pt is not distributed
    // then use the serial get_block method
    if (this->distribution_pt()->communicator_pt()->nproc() == 1 ||
        !this->distribution_pt()->distributed())
    {
      // Vector of vectors for each section of residual vector
      s.resize(nblock);

      // pointer to the data in v
      const double* v_pt = v.values_pt();

      // setup the block vector and then insert the data
      for (unsigned b = 0; b < nblock; b++)
      {
        const unsigned required_block = block_vec_number[b];
        s[b].build(Internal_block_distribution_pt[required_block], 0.0);
        double* s_pt = s[b].values_pt();
        unsigned nrow = s[b].nrow();
        for (unsigned i = 0; i < nrow; i++)
        {
          s_pt[i] = v_pt[this->Global_index[required_block][i]];
        }
      }
    }
    // otherwise use mpi
    else
    {
#ifdef OOMPH_HAS_MPI
      // my rank
      unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

      // the number of processors
      unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

      // build the vectors
      s.resize(nblock);
      for (unsigned b = 0; b < nblock; b++)
      {
        const unsigned required_block = block_vec_number[b];
        s[b].build(Internal_block_distribution_pt[required_block], 0.0);
      }

      // determine the maximum number of rows to be sent or recv
      // and determine the number of blocks each processor will send and recv
      // communication for
      Vector<int> nblock_send(nproc, 0);
      Vector<int> nblock_recv(nproc, 0);
      unsigned max_n_send_or_recv = 0;
      for (unsigned p = 0; p < nproc; p++)
      {
        for (unsigned b = 0; b < nblock; b++)
        {
          const unsigned required_block = block_vec_number[b];
          max_n_send_or_recv = std::max(
            max_n_send_or_recv, Nrows_to_send_for_get_block(required_block, p));
          max_n_send_or_recv = std::max(
            max_n_send_or_recv, Nrows_to_recv_for_get_block(required_block, p));
          if (Nrows_to_send_for_get_block(required_block, p) > 0)
          {
            nblock_send[p]++;
          }
          if (Nrows_to_recv_for_get_block(required_block, p) > 0)
          {
            nblock_recv[p]++;
          }
        }
      }

      // create a vectors of 1s the size of the nblock for the mpi indexed
      // data types
      int* block_lengths = new int[max_n_send_or_recv];
      for (unsigned i = 0; i < max_n_send_or_recv; i++)
      {
        block_lengths[i] = 1;
      }

      // perform the sends and receives
      Vector<MPI_Request> requests;
      for (unsigned p = 0; p < nproc; p++)
      {
        // send and recv with other processors
        if (p != my_rank)
        {
          // send
          if (nblock_send[p] > 0)
          {
            // create the datatypes vector
            MPI_Datatype block_send_types[nblock_send[p]];

            // create the datatypes
            unsigned ptr = 0;
            for (unsigned b = 0; b < nblock; b++)
            {
              const unsigned required_block = block_vec_number[b];

              if (Nrows_to_send_for_get_block(required_block, p) > 0)
              {
                MPI_Type_indexed(Nrows_to_send_for_get_block(required_block, p),
                                 block_lengths,
                                 Rows_to_send_for_get_block(required_block, p),
                                 MPI_DOUBLE,
                                 &block_send_types[ptr]);
                MPI_Type_commit(&block_send_types[ptr]);
                ptr++;
              }
            }

            // compute the displacements and lengths
            MPI_Aint displacements[nblock_send[p]];
            int lengths[nblock_send[p]];
            for (int i = 0; i < nblock_send[p]; i++)
            {
              lengths[i] = 1;
              displacements[i] = 0;
            }

            // build the final datatype
            MPI_Datatype type_send;
            MPI_Type_create_struct(nblock_send[p],
                                   lengths,
                                   displacements,
                                   block_send_types,
                                   &type_send);
            MPI_Type_commit(&type_send);

            // send
            MPI_Request send_req;
            MPI_Isend(const_cast<double*>(v.values_pt()),
                      1,
                      type_send,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &send_req);
            MPI_Type_free(&type_send);
            for (int i = 0; i < nblock_send[p]; i++)
            {
              MPI_Type_free(&block_send_types[i]);
            }
            requests.push_back(send_req);
          }

          // recv
          if (nblock_recv[p] > 0)
          {
            // create the datatypes vector
            MPI_Datatype block_recv_types[nblock_recv[p]];

            // and the displacements
            MPI_Aint displacements[nblock_recv[p]];

            // and the lengths
            int lengths[nblock_recv[p]];

            // all displacements are computed relative to s[0] values
            MPI_Aint displacements_base;
            MPI_Get_address(s[0].values_pt(), &displacements_base);

            // now build
            unsigned ptr = 0;
            for (unsigned b = 0; b < nblock; b++)
            {
              const unsigned required_block = block_vec_number[b];

              if (Nrows_to_recv_for_get_block(required_block, p) > 0)
              {
                MPI_Type_indexed(Nrows_to_recv_for_get_block(required_block, p),
                                 block_lengths,
                                 Rows_to_recv_for_get_block(required_block, p),
                                 MPI_DOUBLE,
                                 &block_recv_types[ptr]);
                MPI_Type_commit(&block_recv_types[ptr]);
                MPI_Get_address(s[b].values_pt(), &displacements[ptr]);
                displacements[ptr] -= displacements_base;
                lengths[ptr] = 1;
                ptr++;
              }
            }

            // build the final data type
            MPI_Datatype type_recv;
            MPI_Type_create_struct(nblock_recv[p],
                                   lengths,
                                   displacements,
                                   block_recv_types,
                                   &type_recv);
            MPI_Type_commit(&type_recv);

            // recv
            MPI_Request recv_req;
            MPI_Irecv(s[0].values_pt(),
                      1,
                      type_recv,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &recv_req);
            MPI_Type_free(&type_recv);
            for (int i = 0; i < nblock_recv[p]; i++)
            {
              MPI_Type_free(&block_recv_types[i]);
            }
            requests.push_back(recv_req);
          }
        }

        // communicate with self
        else
        {
          const double* v_values_pt = v.values_pt();
          for (unsigned b = 0; b < nblock; b++)
          {
            const unsigned required_block = block_vec_number[b];

            double* w_values_pt = s[b].values_pt();
            for (unsigned i = 0;
                 i < Nrows_to_send_for_get_block(required_block, p);
                 i++)
            {
              w_values_pt[Rows_to_recv_for_get_block(required_block, p)[i]] =
                v_values_pt[Rows_to_send_for_get_block(required_block, p)[i]];
            }
          }
        }
      }

      // and then just wait
      unsigned c = requests.size();
      Vector<MPI_Status> stat(c);
      if (c)
      {
        MPI_Waitall(c, &requests[0], &stat[0]);
      }
      delete[] block_lengths;

#else
      // throw error
      std::ostringstream error_message;
      error_message << "The preconditioner is distributed and on more than one "
                    << "processor. MPI is required.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    }
  }

  //============================================================================
  /// A helper function, takes the naturally ordered vector and
  /// rearranges it into a vector of sub vectors corresponding to the blocks,
  /// so s[b][i] contains the i-th entry in the vector associated with block b.
  /// The block_vec_number indicates which blocks we want.
  /// These blocks and vectors are those corresponding to the internal blocks.
  /// Note: If the preconditioner is a subsidiary preconditioner then only the
  /// sub-vectors associated with the blocks of the subsidiary preconditioner
  /// will be included. Hence the length of v is master_nrow() whereas the
  /// total length of the s vectors is the sum of the Nrow of the sub vectors.
  /// This is simply a wrapper around the other internal_get_block_vectors(...)
  /// function with the identity block_vec_number vector.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::internal_get_block_vectors(
    const DoubleVector& v, Vector<DoubleVector>& s) const
  {
    // Number of block types
    const unsigned nblock = this->internal_nblock_types();
    Vector<unsigned> block_vec_number(nblock, 0);
    for (unsigned b = 0; b < nblock; b++)
    {
      block_vec_number[b] = b;
    }

    internal_get_block_vectors(block_vec_number, v, s);
  }

  //============================================================================
  /// Takes the vector of block vectors, s, and copies its entries into
  /// the naturally ordered vector, v. If this is a subsidiary block
  /// preconditioner only those entries in v that are associated with its
  /// blocks are affected. The block_vec_number indicates which block the
  /// vectors in s came from. The block number corresponds to the block
  /// numbers in this preconditioner.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::return_block_vectors(
    const Vector<unsigned>& block_vec_number,
    const Vector<DoubleVector>& s,
    DoubleVector& v) const
  {
#ifdef PARANOID

    // Check if v is built.
    if (!v.built())
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // v must have the same distribution as the upper-most master block
    // preconditioner, since the upper-most master block preconditioner
    // should have the same distribution as the matrix pointed to
    // by matrix_pt().
    if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the global vector v must match the "
              << " specified master_distribution_pt(). \n"
              << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check if the number of vectors in s is the same as the number of block
    // numbers described in block_vec_number.
    const unsigned para_block_vec_number_size = block_vec_number.size();
    const unsigned para_s_size = s.size();
    if (para_block_vec_number_size != para_s_size)
    {
      std::ostringstream err_msg;
      err_msg << "block_vec_number.size() is " << para_block_vec_number_size
              << "\n."
              << "s.size() is " << para_s_size << ".\n"
              << "But they must be the same size!\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check to see if there are more blocks defined in the block_vec_number
    // vector than the number of block types. This is not allowed.
    const unsigned para_n_block = nblock_types();
    if (para_block_vec_number_size > para_n_block)
    {
      std::ostringstream err_msg;
      err_msg << "Trying to return " << para_block_vec_number_size
              << " block vectors.\n"
              << "But there are only " << para_n_block << " block types.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check if any block numbers defined in block_vec_number is equal to or
    // greater than the number of block types.
    // E.g. if there are 5 block types, we can only have block numbers:
    //  0, 1, 2, 3 and 4.
    for (unsigned b = 0; b < para_block_vec_number_size; b++)
    {
      const unsigned para_required_block = block_vec_number[b];
      if (para_required_block > para_n_block)
      {
        std::ostringstream err_msg;
        err_msg << "block_vec_number[" << b << "] is " << para_required_block
                << ".\n"
                << "But there are only " << para_n_block << " block types.\n";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }

    // Check that no block number is inserted twice.
    std::set<unsigned> para_set;
    for (unsigned b = 0; b < para_block_vec_number_size; b++)
    {
      std::pair<std::set<unsigned>::iterator, bool> para_set_ret;
      para_set_ret = para_set.insert(block_vec_number[b]);

      if (!para_set_ret.second)
      {
        std::ostringstream err_msg;
        err_msg << "Error: the block number " << block_vec_number[b]
                << " appears twice.\n";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }

    // Check to see that all the vectors in s are built
    // (since we are trying to return them).
    for (unsigned b = 0; b < para_block_vec_number_size; b++)
    {
      if (!s[b].built())
      {
        std::ostringstream err_msg;
        err_msg << "The distribution of the block vector s[" << b
                << "] must be setup.\n";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }

    // Since these are built, we check that the distributions are correct.
    // This are incorrect if the block numbers in block_vec_number and
    // the vectors in s does not match.
    for (unsigned b = 0; b < para_block_vec_number_size; b++)
    {
      if (*(s[b].distribution_pt()) !=
          *(Block_distribution_pt[block_vec_number[b]]))
      {
        std::ostringstream error_message;
        error_message
          << "The distribution of the block vector " << b << " must match the"
          << " specified distribution at "
          << "Block_distribution_pt[" << block_vec_number[b] << "].\n"
          << "The distribution of the Block_distribution_pt is determined by\n"
          << "the vector block_vec_number. Perhaps it is incorrect?\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Number of blocks to get.
    const unsigned n_block = block_vec_number.size();

    // Each block is made of dof types. We get the most fine grain dof types.
    // Most fine grain in the sense that these are the dof types that belongs
    // in this block before any coarsening of dof types has taken place.
    // The ordering of the dof types matters, this is handled properly when
    // creating the Block_to_dof_map_fine vector and must be respected here.
    // I.e. we cannot arbitrarily insert dof types (even if they are correct)
    // in the vector most_fine_grain_dof.
    Vector<unsigned> most_fine_grain_dof;
    for (unsigned b = 0; b < n_block; b++)
    {
      const unsigned mapped_b = block_vec_number[b];

      most_fine_grain_dof.insert(most_fine_grain_dof.end(),
                                 Block_to_dof_map_fine[mapped_b].begin(),
                                 Block_to_dof_map_fine[mapped_b].end());
    }

    // Split all the blocks into it's most fine grain dof vector.
    Vector<DoubleVector> dof_vector(most_fine_grain_dof.size());

    unsigned offset = 0;

    // Perform the splitting for each block.
    for (unsigned b = 0; b < n_block; b++)
    {
      // The actual block number.
      const unsigned mapped_b = block_vec_number[b];

      // How many most fine grain dof types are associated with this block?
      const unsigned ndof = Block_to_dof_map_fine[mapped_b].size();

      if (ndof == 1)
      // No need to split, just copy.
      {
        dof_vector[offset] = s[b];
      }
      else
      // Need to split s[b] into it's most fine grain dof vectors
      {
        // To store pointers to the dof vectors associated with this block.
        Vector<DoubleVector*> tmp_dof_vector_pt(ndof, 0);

        for (unsigned d = 0; d < ndof; d++)
        {
          const unsigned offset_plus_d = offset + d;

          // build the dof vector.
          dof_vector[offset_plus_d].build(
            Internal_block_distribution_pt[most_fine_grain_dof[offset_plus_d]]);

          // Store the pointer.
          tmp_dof_vector_pt[d] = &dof_vector[offset_plus_d];
        }

        // Split without communication.
        DoubleVectorHelpers::split_without_communication(s[b],
                                                         tmp_dof_vector_pt);
      }

      // Update the offset!
      offset += ndof;
    }

    // Return the block vectors all in one go.
    internal_return_block_vectors(most_fine_grain_dof, dof_vector, v);
  } // return_block_vectors(...)


  //============================================================================
  /// Takes the vector of block vectors, s, and copies its entries into
  /// the naturally ordered vector, v. If this is a subsidiary block
  /// preconditioner only those entries in v that are associated with its
  /// blocks are affected. The block_vec_number indicates which block the
  /// vectors in s came from. The block number corresponds to the block
  /// numbers in this preconditioner.
  /// This is simply a wrapper around the other return_block_vectors(...)
  /// function where the block_vec_number Vector is the identity, i.e.
  /// block_vec_number is [0, 1, ..., nblock_types - 1].
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::return_block_vectors(
    const Vector<DoubleVector>& s, DoubleVector& v) const
  {
    // The number of block types in this preconditioner.
    const unsigned n_block = nblock_types();

    // Create the identity vector.
    Vector<unsigned> required_block(n_block, 0);
    for (unsigned i = 0; i < n_block; i++)
    {
      required_block[i] = i;
    }

    // Call the other return_block_vectors function which does the work.
    return_block_vectors(required_block, s, v);
  } // return_block_vectors(...)

  //============================================================================
  /// Takes the naturally ordered vector and
  /// rearranges it into a vector of sub vectors corresponding to the blocks,
  /// so s[b][i] contains the i-th entry in the vector associated with block b.
  /// The block_vec_number indicates which blocks we want.
  /// These blocks and vectors are those corresponding to the internal blocks.
  /// Note: If the preconditioner is a subsidiary preconditioner then only the
  /// sub-vectors associated with the blocks of the subsidiary preconditioner
  /// will be included. Hence the length of v is master_nrow() whereas the
  /// total length of the s vectors is the sum of the Nrow of the sub vectors.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::internal_return_block_vectors(
    const Vector<unsigned>& block_vec_number,
    const Vector<DoubleVector>& s,
    DoubleVector& v) const
  {
    // the number of blocks
    const unsigned nblock = block_vec_number.size();

#ifdef PARANOID
    if (!v.built())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must match the "
                    << " specified master_distribution_pt(). \n"
                    << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    for (unsigned b = 0; b < nblock; b++)
    {
      if (!s[b].built())
      {
        std::ostringstream error_message;
        error_message << "The distribution of the block vector " << b
                      << " must be setup.";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      const unsigned required_block = block_vec_number[b];
      if (*(s[b].distribution_pt()) !=
          *(Internal_block_distribution_pt[required_block]))
      {
        std::ostringstream error_message;
        error_message
          << "The distribution of the block vector " << b << " must match the"
          << " specified distribution at Internal_block_distribution_pt[" << b
          << "]";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // if + only one processor
    //    + more than one processor but matrix_pt is not distributed
    // then use the serial get_block method
    if (this->distribution_pt()->communicator_pt()->nproc() == 1 ||
        !this->distribution_pt()->distributed())
    {
      double* v_pt = v.values_pt();
      for (unsigned b = 0; b < nblock; b++)
      {
        const unsigned required_block = block_vec_number[b];

        const double* s_pt = s[b].values_pt();
        unsigned nrow = this->internal_block_dimension(required_block);
        for (unsigned i = 0; i < nrow; i++)
        {
          v_pt[this->Global_index[required_block][i]] = s_pt[i];
        }
      }
    }
    // otherwise use mpi
    else
    {
#ifdef OOMPH_HAS_MPI

      // my rank
      unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

      // the number of processors
      unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

      // determine the maximum number of rows to be sent or recv
      // and determine the number of blocks each processor will send and recv
      // communication for
      Vector<int> nblock_send(nproc, 0);
      Vector<int> nblock_recv(nproc, 0);
      unsigned max_n_send_or_recv = 0;
      for (unsigned p = 0; p < nproc; p++)
      {
        for (unsigned b = 0; b < nblock; b++)
        {
          const unsigned required_block = block_vec_number[b];

          max_n_send_or_recv = std::max(
            max_n_send_or_recv, Nrows_to_send_for_get_block(required_block, p));
          max_n_send_or_recv = std::max(
            max_n_send_or_recv, Nrows_to_recv_for_get_block(required_block, p));
          if (Nrows_to_send_for_get_block(required_block, p) > 0)
          {
            nblock_recv[p]++;
          }
          if (Nrows_to_recv_for_get_block(required_block, p) > 0)
          {
            nblock_send[p]++;
          }
        }
      }

      // create a vectors of 1s the size of the nblock for the mpi indexed
      // data types
      int* block_lengths = new int[max_n_send_or_recv];
      for (unsigned i = 0; i < max_n_send_or_recv; i++)
      {
        block_lengths[i] = 1;
      }

      // perform the sends and receives
      Vector<MPI_Request> requests;
      for (unsigned p = 0; p < nproc; p++)
      {
        // send and recv with other processors
        if (p != my_rank)
        {
          // recv
          if (nblock_recv[p] > 0)
          {
            // create the datatypes vector
            MPI_Datatype block_recv_types[nblock_recv[p]];

            // create the datatypes
            unsigned ptr = 0;
            for (unsigned b = 0; b < nblock; b++)
            {
              const unsigned required_block = block_vec_number[b];

              if (Nrows_to_send_for_get_block(required_block, p) > 0)
              {
                MPI_Type_indexed(Nrows_to_send_for_get_block(required_block, p),
                                 block_lengths,
                                 Rows_to_send_for_get_block(required_block, p),
                                 MPI_DOUBLE,
                                 &block_recv_types[ptr]);
                MPI_Type_commit(&block_recv_types[ptr]);
                ptr++;
              }
            }

            // compute the displacements and lengths
            MPI_Aint displacements[nblock_recv[p]];
            int lengths[nblock_recv[p]];
            for (int i = 0; i < nblock_recv[p]; i++)
            {
              lengths[i] = 1;
              displacements[i] = 0;
            }

            // build the final datatype
            MPI_Datatype type_recv;
            MPI_Type_create_struct(nblock_recv[p],
                                   lengths,
                                   displacements,
                                   block_recv_types,
                                   &type_recv);
            MPI_Type_commit(&type_recv);

            // recv
            MPI_Request recv_req;
            MPI_Irecv(v.values_pt(),
                      1,
                      type_recv,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &recv_req);
            MPI_Type_free(&type_recv);
            for (int i = 0; i < nblock_recv[p]; i++)
            {
              MPI_Type_free(&block_recv_types[i]);
            }
            requests.push_back(recv_req);
          }

          // send
          if (nblock_send[p] > 0)
          {
            // create the datatypes vector
            MPI_Datatype block_send_types[nblock_send[p]];

            // and the displacements
            MPI_Aint displacements[nblock_send[p]];

            // and the lengths
            int lengths[nblock_send[p]];

            // all displacements are computed relative to s[0] values
            MPI_Aint displacements_base;
            MPI_Get_address(const_cast<double*>(s[0].values_pt()),
                            &displacements_base);

            // now build
            unsigned ptr = 0;
            for (unsigned b = 0; b < nblock; b++)
            {
              const unsigned required_block = block_vec_number[b];

              if (Nrows_to_recv_for_get_block(required_block, p) > 0)
              {
                MPI_Type_indexed(Nrows_to_recv_for_get_block(required_block, p),
                                 block_lengths,
                                 Rows_to_recv_for_get_block(required_block, p),
                                 MPI_DOUBLE,
                                 &block_send_types[ptr]);
                MPI_Type_commit(&block_send_types[ptr]);
                MPI_Get_address(const_cast<double*>(s[b].values_pt()),
                                &displacements[ptr]);
                displacements[ptr] -= displacements_base;
                lengths[ptr] = 1;
                ptr++;
              }
            }

            // build the final data type
            MPI_Datatype type_send;
            MPI_Type_create_struct(nblock_send[p],
                                   lengths,
                                   displacements,
                                   block_send_types,
                                   &type_send);
            MPI_Type_commit(&type_send);

            // send
            MPI_Request send_req;
            MPI_Isend(const_cast<double*>(s[0].values_pt()),
                      1,
                      type_send,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &send_req);
            MPI_Type_free(&type_send);
            for (int i = 0; i < nblock_send[p]; i++)
            {
              MPI_Type_free(&block_send_types[i]);
            }
            requests.push_back(send_req);
          }
        }

        // communicate wih self
        else
        {
          double* v_values_pt = v.values_pt();
          for (unsigned b = 0; b < nblock; b++)
          {
            const unsigned required_block = block_vec_number[b];

            const double* w_values_pt = s[b].values_pt();
            for (unsigned i = 0;
                 i < Nrows_to_send_for_get_block(required_block, p);
                 i++)
            {
              v_values_pt[Rows_to_send_for_get_block(required_block, p)[i]] =
                w_values_pt[Rows_to_recv_for_get_block(required_block, p)[i]];
            }
          }
        }
      }

      // and then just wait
      unsigned c = requests.size();
      Vector<MPI_Status> stat(c);
      if (c)
      {
        MPI_Waitall(c, &requests[0], &stat[0]);
      }
      delete[] block_lengths;

#else
      // throw error
      std::ostringstream error_message;
      error_message << "The preconditioner is distributed and on more than one "
                    << "processor. MPI is required.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    }
  }

  //============================================================================
  /// A helper function, takes the naturally ordered vector and
  /// rearranges it into a vector of sub vectors corresponding to the blocks,
  /// so s[b][i] contains the i-th entry in the vector associated with block b.
  /// The block_vec_number indicates which blocks we want.
  /// These blocks and vectors are those corresponding to the internal blocks.
  /// Note: If the preconditioner is a subsidiary preconditioner then only the
  /// sub-vectors associated with the blocks of the subsidiary preconditioner
  /// will be included. Hence the length of v is master_nrow() whereas the
  /// total length of the s vectors is the sum of the Nrow of the sub vectors.
  /// This is simply a wrapper around the other internal_get_block_vectors(...)
  /// function with the identity block_vec_number vector.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::internal_return_block_vectors(
    const Vector<DoubleVector>& s, DoubleVector& v) const
  {
    // the number of blocks
    const unsigned nblock = this->internal_nblock_types();
    Vector<unsigned> block_vec_number(nblock, 0);
    for (unsigned b = 0; b < nblock; b++)
    {
      block_vec_number[b] = b;
    }

    internal_return_block_vectors(block_vec_number, s, v);
  }

  //============================================================================
  /// A helper function, takes the naturally ordered vector, v,
  /// and extracts the n-th block vector, b.
  /// Here n is the block number in the current preconditioner.
  /// NOTE: The ordering of the vector b is the same as the
  /// ordering of the block matrix from internal_get_block(...).
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::internal_get_block_vector(
    const unsigned& b, const DoubleVector& v, DoubleVector& w) const
  {
#ifdef PARANOID
    // the number of blocks
    const unsigned n_blocks = this->internal_nblock_types();

    // paranoid check that block i is in this block preconditioner
    if (b >= n_blocks)
    {
      std::ostringstream error_message;
      error_message
        << "Requested block  vector " << b
        << ", however this preconditioner has internal_nblock_types() "
        << "= " << internal_nblock_types() << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (!v.built())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must match the "
                    << " specified master_distribution_pt(). \n"
                    << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // rebuild the block vector
    w.build(Internal_block_distribution_pt[b], 0.0);

    // if + only one processor
    //    + more than one processor but matrix_pt is not distributed
    // then use the serial get_block method
    if (this->distribution_pt()->communicator_pt()->nproc() == 1 ||
        !this->distribution_pt()->distributed())
    {
      double* w_pt = w.values_pt();
      const double* v_pt = v.values_pt();
      unsigned n_row = w.nrow();
      for (unsigned i = 0; i < n_row; i++)
      {
        w_pt[i] = v_pt[this->Global_index[b][i]];
      }
    }
    // otherwise use mpi
    else
    {
#ifdef OOMPH_HAS_MPI

      // my rank
      unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

      // the number of processors
      unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

      // determine the maximum number of rows to be sent or recv
      unsigned max_n_send_or_recv = 0;
      for (unsigned p = 0; p < nproc; p++)
      {
        max_n_send_or_recv =
          std::max(max_n_send_or_recv, Nrows_to_send_for_get_block(b, p));
        max_n_send_or_recv =
          std::max(max_n_send_or_recv, Nrows_to_recv_for_get_block(b, p));
      }

      // create a vectors of 1s (the size of the nblock for the mpi indexed
      // data types
      int* block_lengths = new int[max_n_send_or_recv];
      for (unsigned i = 0; i < max_n_send_or_recv; i++)
      {
        block_lengths[i] = 1;
      }

      // perform the sends and receives
      Vector<MPI_Request> requests;
      for (unsigned p = 0; p < nproc; p++)
      {
        // send and recv with other processors
        if (p != my_rank)
        {
          if (Nrows_to_send_for_get_block(b, p) > 0)
          {
            // create the send datatype
            MPI_Datatype type_send;
            MPI_Type_indexed(Nrows_to_send_for_get_block(b, p),
                             block_lengths,
                             Rows_to_send_for_get_block(b, p),
                             MPI_DOUBLE,
                             &type_send);
            MPI_Type_commit(&type_send);

            // send
            MPI_Request send_req;
            MPI_Isend(const_cast<double*>(v.values_pt()),
                      1,
                      type_send,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &send_req);
            MPI_Type_free(&type_send);
            requests.push_back(send_req);
          }

          if (Nrows_to_recv_for_get_block(b, p) > 0)
          {
            // create the recv datatype
            MPI_Datatype type_recv;
            MPI_Type_indexed(Nrows_to_recv_for_get_block(b, p),
                             block_lengths,
                             Rows_to_recv_for_get_block(b, p),
                             MPI_DOUBLE,
                             &type_recv);
            MPI_Type_commit(&type_recv);

            // recv
            MPI_Request recv_req;
            MPI_Irecv(w.values_pt(),
                      1,
                      type_recv,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &recv_req);
            MPI_Type_free(&type_recv);
            requests.push_back(recv_req);
          }
        }

        // communicate with self
        else
        {
          double* w_values_pt = w.values_pt();
          const double* v_values_pt = v.values_pt();
          for (unsigned i = 0; i < Nrows_to_send_for_get_block(b, p); i++)
          {
            w_values_pt[Rows_to_recv_for_get_block(b, p)[i]] =
              v_values_pt[Rows_to_send_for_get_block(b, p)[i]];
          }
        }
      }

      // and then just wait
      unsigned c = requests.size();
      Vector<MPI_Status> stat(c);
      if (c)
      {
        MPI_Waitall(c, &requests[0], &stat[0]);
      }
      delete[] block_lengths;

#else
      // throw error
      std::ostringstream error_message;
      error_message << "The preconditioner is distributed and on more than one "
                    << "processor. MPI is required.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    }
  }

  //============================================================================
  /// Takes the naturally ordered vector, v and returns the n-th
  /// block vector, b. Here n is the block number in the current
  /// preconditioner.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::get_block_vector(const unsigned& b,
                                                     const DoubleVector& v,
                                                     DoubleVector& w) const
  {
#ifdef PARANOID
    // the number of blocks
    const unsigned para_n_blocks = nblock_types();

    // paranoid check that block i is in this block preconditioner
    if (b >= para_n_blocks)
    {
      std::ostringstream err_msg;
      err_msg << "Requested block vector " << b
              << ", however this preconditioner has only " << para_n_blocks
              << " block types"
              << ".\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    if (!v.built())
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*(v.distribution_pt()) != *(this->master_distribution_pt()))
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the global vector v must match the "
              << " specified master_distribution_pt(). \n"
              << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Recall that, the relationship between the external blocks and the
    // external dof types, as seen by the preconditioner writer is stored in the
    // mapping Block_to_dof_map_coarse.
    //
    // However, each dof type could have been coarsened! The relationship
    // between the dof types of this preconditioner and the parent
    // preconditioner is stored in the mapping Doftype_coarsen_map_coarse. The
    // dof numbers in this map is relative to this preconditioner.
    //
    // Finally, the relationship between the dof types of this preconditioner
    // and the most fine grain dof types is stored in the mapping
    // Doftype_coarsen_map_fine. Again, the dof numbers in this map is relative
    // to this preconditioner.
    //
    // Furthermore, we note that concatenation of vectors without communication
    //  is associative, but not commutative. I.e.
    // (V1+V2)+V3 = V1 + (V2 + V3), where + is concatenation without
    // communication.
    //
    // So all we need is the vectors listed in the correct order.
    //
    // We need only Block_to_dof_map_coarse to tell us which external dof types
    // are in this block, then Doftype_coarsen_map_fine to tell us which most
    // fine grain dofs to concatenate!
    //
    // All the mapping vectors are constructed to respect the ordering of
    // the dof types.

    // Get the most fine grain block to dof mapping.
    Vector<unsigned> most_fine_grain_dof = Block_to_dof_map_fine[b];

    // How many vectors do we need to concatenate?
    const unsigned n_dof_vec = most_fine_grain_dof.size();

    if (n_dof_vec == 1)
    // No need to concatenate, just extract the vector.
    {
      internal_get_block_vector(most_fine_grain_dof[0], v, w);
    }
    else
    // Need to concatenate dof-level vectors.
    {
      Vector<DoubleVector> dof_vector(n_dof_vec);

      // Get all the dof-level vectors in one go
      internal_get_block_vectors(most_fine_grain_dof, v, dof_vector);
      // Build w with the correct distribution.
      w.build(Block_distribution_pt[b], 0);

      // Concatenate the vectors.
      DoubleVectorHelpers::concatenate_without_communication(dof_vector, w);

      dof_vector.clear();
    }
  } // get_block_vector(...)

  //============================================================================
  /// Takes the n-th block ordered vector, b,  and copies its entries
  /// to the appropriate entries in the naturally ordered vector, v.
  /// Here n is the block number in the current block preconditioner.
  /// If the preconditioner is a subsidiary block preconditioner
  /// the other entries in v  that are not associated with it
  /// are left alone.
  ///
  /// This version works with the internal block types. This is legacy code
  /// but is kept alive, hence moved to private. Please use the
  /// function "return_block_vector(...)".
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::internal_return_block_vector(
    const unsigned& b, const DoubleVector& w, DoubleVector& v) const
  {
#ifdef PARANOID
    // the number of blocks
    const unsigned n_blocks = this->internal_nblock_types();

    // paranoid check that block i is in this block preconditioner
    if (b >= n_blocks)
    {
      std::ostringstream error_message;
      error_message
        << "Requested block  vector " << b
        << ", however this preconditioner has internal_nblock_types() "
        << "= " << internal_nblock_types() << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (!v.built())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*v.distribution_pt() != *this->master_distribution_pt())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must match the "
                    << " specified master_distribution_pt(). \n"
                    << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (!w.built())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the block vector w must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*w.distribution_pt() != *Internal_block_distribution_pt[b])
    {
      std::ostringstream error_message;
      error_message
        << "The distribution of the block vector w must match the "
        << " specified distribution at Internal_block_distribution_pt[b]";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // if + only one processor
    //    + more than one processor but matrix_pt is not distributed
    // then use the serial get_block method
    if (this->distribution_pt()->communicator_pt()->nproc() == 1 ||
        !this->distribution_pt()->distributed())
    {
      // length of vector
      unsigned n_row = this->internal_block_dimension(b);

      // copy back from the block vector to the naturally ordered vector
      double* v_pt = v.values_pt();
      const double* w_pt = w.values_pt();
      for (unsigned i = 0; i < n_row; i++)
      {
        v_pt[this->Global_index[b][i]] = w_pt[i];
      }
    }
    // otherwise use mpi
    else
    {
#ifdef OOMPH_HAS_MPI

      // my rank
      unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

      // the number of processors
      unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

      // determine the maximum number of rows to be sent or recv
      unsigned max_n_send_or_recv = 0;
      for (unsigned p = 0; p < nproc; p++)
      {
        max_n_send_or_recv =
          std::max(max_n_send_or_recv, Nrows_to_send_for_get_block(b, p));
        max_n_send_or_recv =
          std::max(max_n_send_or_recv, Nrows_to_recv_for_get_block(b, p));
      }

      // create a vectors of 1s (the size of the nblock for the mpi indexed
      // data types
      int* block_lengths = new int[max_n_send_or_recv];
      for (unsigned i = 0; i < max_n_send_or_recv; i++)
      {
        block_lengths[i] = 1;
      }

      // perform the sends and receives
      Vector<MPI_Request> requests;
      for (unsigned p = 0; p < nproc; p++)
      {
        // send and recv with other processors
        if (p != my_rank)
        {
          if (Nrows_to_recv_for_get_block(b, p) > 0)
          {
            // create the send datatype
            MPI_Datatype type_send;
            MPI_Type_indexed(Nrows_to_recv_for_get_block(b, p),
                             block_lengths,
                             Rows_to_recv_for_get_block(b, p),
                             MPI_DOUBLE,
                             &type_send);
            MPI_Type_commit(&type_send);

            // send
            MPI_Request send_req;
            MPI_Isend(const_cast<double*>(w.values_pt()),
                      1,
                      type_send,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &send_req);
            MPI_Type_free(&type_send);
            requests.push_back(send_req);
          }

          if (Nrows_to_send_for_get_block(b, p) > 0)
          {
            // create the recv datatype
            MPI_Datatype type_recv;
            MPI_Type_indexed(Nrows_to_send_for_get_block(b, p),
                             block_lengths,
                             Rows_to_send_for_get_block(b, p),
                             MPI_DOUBLE,
                             &type_recv);
            MPI_Type_commit(&type_recv);

            // recv
            MPI_Request recv_req;
            MPI_Irecv(v.values_pt(),
                      1,
                      type_recv,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &recv_req);
            MPI_Type_free(&type_recv);
            requests.push_back(recv_req);
          }
        }

        // communicate wih self
        else
        {
          const double* w_values_pt = w.values_pt();
          double* v_values_pt = v.values_pt();
          for (unsigned i = 0; i < Nrows_to_send_for_get_block(b, p); i++)
          {
            v_values_pt[Rows_to_send_for_get_block(b, p)[i]] =
              w_values_pt[Rows_to_recv_for_get_block(b, p)[i]];
          }
        }
      }

      // and then just wait
      unsigned c = requests.size();
      Vector<MPI_Status> stat(c);
      if (c)
      {
        MPI_Waitall(c, &requests[0], &stat[0]);
      }
      delete[] block_lengths;

#else
      // throw error
      std::ostringstream error_message;
      error_message << "The preconditioner is distributed and on more than one "
                    << "processor. MPI is required.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    }
  }

  //============================================================================
  /// Takes the n-th block ordered vector, b,  and copies its entries
  /// to the appropriate entries in the naturally ordered vector, v.
  /// Here n is the block number in the current block preconditioner.
  /// If the preconditioner is a subsidiary block preconditioner
  /// the other entries in v  that are not associated with it
  /// are left alone.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::return_block_vector(const unsigned& n,
                                                        const DoubleVector& b,
                                                        DoubleVector& v) const
  {
#ifdef PARANOID
    // the number of blocks
    const unsigned para_n_blocks = nblock_types();

    // paranoid check that block i is in this block preconditioner
    if (n >= para_n_blocks)
    {
      std::ostringstream err_msg;
      err_msg << "Requested block vector " << b
              << ", however this preconditioner has " << para_n_blocks
              << " block types.\n";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (!v.built())
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*v.distribution_pt() != *this->master_distribution_pt())
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the global vector v must match the "
              << " specified master_distribution_pt(). \n"
              << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (!b.built())
    {
      std::ostringstream err_msg;
      err_msg << "The distribution of the block vector b must be setup.";
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

#endif

    // Get the most fine grain dof
    Vector<unsigned> most_fine_grain_dof = Block_to_dof_map_fine[n];

    // How many dofs are in this block?
    const unsigned n_dof_vec = Block_to_dof_map_fine[n].size();

    if (n_dof_vec == 1)
    // There is only one dof, no need to split.
    {
      internal_return_block_vector(most_fine_grain_dof[0], b, v);
    }
    else
    // Need to split the vector up before we insert them all in one go.
    {
      Vector<DoubleVector> dof_vector(n_dof_vec);
      for (unsigned d = 0; d < n_dof_vec; d++)
      {
        dof_vector[d].build(
          internal_block_distribution_pt(most_fine_grain_dof[d]));
      }

      DoubleVectorHelpers::split_without_communication(b, dof_vector);

      // return to v
      internal_return_block_vectors(most_fine_grain_dof, dof_vector, v);
    }
  } // return_block_vector(...)

  //============================================================================
  /// Given the naturally ordered vector, v, return
  /// the vector rearranged in block order in w. This is a legacy function
  /// from the old block preconditioning framework. Kept alive in case it may
  /// be needed again.
  ///
  /// This uses the variables ending in "get_ordered". We no longer use this
  /// type of method. This function copy values from v and re-order them
  /// in "block order" and place them in w. Block order means that the
  /// values in w are the same as the concatenated block vectors.
  ///
  /// I.e. - v is naturally ordered.
  ///        v -> s_b, v is ordered into blocks vectors
  ///                  (requires communication)
  ///        concatenate_without_communication(s_{0,...,nblocks},w) gives w.
  ///
  /// But this function skips out the concatenation part and builds w directly
  /// from v.
  ///
  /// This is nice but the function is implemented in such a way that it
  /// always use all the (internal) blocks and concatenated with the
  /// identity ordering. I.e. if this preconditioner has 3 block types, then
  /// w will always be:
  /// concatenate_without_communication([s_0, s_1, s_2], w). There is easy
  /// way to change this.
  ///
  /// Furthermore, it does not take into account the new dof type coarsening
  /// feature. So this function will most likely produce the incorrect vector
  /// w from what the user intended. It still works, but w will be the
  /// concatenation of the most fine grain dof block vectors with the
  /// "natural" dof type ordering.
  ///
  /// This has been superseded by the function
  /// get_block_ordered_preconditioner_vector(...) which does the correct
  /// thing.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::
    internal_get_block_ordered_preconditioner_vector(const DoubleVector& v,
                                                     DoubleVector& w) const
  {
#ifdef PARANOID
    if (!v.built())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*v.distribution_pt() != *this->master_distribution_pt())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must match the "
                    << " specified master_distribution_pt(). \n"
                    << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    //  Cleared and resized w for reordered vector
    w.build(this->internal_preconditioner_matrix_distribution_pt(), 0.0);

    // if + only one processor
    //    + more than one processor but matrix_pt is not distributed
    // then use the serial get_block method
    if (this->distribution_pt()->communicator_pt()->nproc() == 1 ||
        !this->distribution_pt()->distributed())
    {
      // number of blocks
      unsigned nblock = this->Internal_nblock_types;

      // copy to w
      unsigned block_offset = 0;
      double* w_pt = w.values_pt();
      const double* v_pt = v.values_pt();
      for (unsigned b = 0; b < nblock; b++)
      {
        unsigned block_nrow = this->internal_block_dimension(b);
        for (unsigned i = 0; i < block_nrow; i++)
        {
          w_pt[block_offset + i] = v_pt[this->Global_index[b][i]];
        }
        block_offset += block_nrow;
      }
    }
    // otherwise use mpi
    else
    {
#ifdef OOMPH_HAS_MPI

      // my rank
      unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

      // the number of processors
      unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

      // determine the maximum number of rows to be sent or recv
      unsigned max_n_send_or_recv = 0;
      for (unsigned p = 0; p < nproc; p++)
      {
        max_n_send_or_recv =
          std::max(max_n_send_or_recv, Nrows_to_send_for_get_ordered[p]);
        max_n_send_or_recv =
          std::max(max_n_send_or_recv, Nrows_to_recv_for_get_ordered[p]);
      }

      // create a vectors of 1s (the size of the nblock for the mpi indexed
      // data types
      int* block_lengths = new int[max_n_send_or_recv];
      for (unsigned i = 0; i < max_n_send_or_recv; i++)
      {
        block_lengths[i] = 1;
      }

      // perform the sends and receives
      Vector<MPI_Request> requests;
      for (unsigned p = 0; p < nproc; p++)
      {
        // send and recv with other processors
        if (p != my_rank)
        {
          if (Nrows_to_send_for_get_ordered[p] > 0)
          {
            // create the send datatype
            MPI_Datatype type_send;
            MPI_Type_indexed(Nrows_to_send_for_get_ordered[p],
                             block_lengths,
                             Rows_to_send_for_get_ordered[p],
                             MPI_DOUBLE,
                             &type_send);
            MPI_Type_commit(&type_send);

            // send
            MPI_Request send_req;
            MPI_Isend(const_cast<double*>(v.values_pt()),
                      1,
                      type_send,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &send_req);
            MPI_Type_free(&type_send);
            requests.push_back(send_req);
          }

          if (Nrows_to_recv_for_get_ordered[p] > 0)
          {
            // create the recv datatype
            MPI_Datatype type_recv;
            MPI_Type_indexed(Nrows_to_recv_for_get_ordered[p],
                             block_lengths,
                             Rows_to_recv_for_get_ordered[p],
                             MPI_DOUBLE,
                             &type_recv);
            MPI_Type_commit(&type_recv);

            // recv
            MPI_Request recv_req;
            MPI_Irecv(w.values_pt(),
                      1,
                      type_recv,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &recv_req);
            MPI_Type_free(&type_recv);
            requests.push_back(recv_req);
          }
        }

        // communicate with self
        else
        {
          double* w_values_pt = w.values_pt();
          const double* v_values_pt = v.values_pt();
          for (unsigned i = 0; i < Nrows_to_send_for_get_ordered[p]; i++)
          {
            w_values_pt[Rows_to_recv_for_get_ordered[p][i]] =
              v_values_pt[Rows_to_send_for_get_ordered[p][i]];
          }
        }
      }

      // and then just wait
      unsigned c = requests.size();
      Vector<MPI_Status> stat(c);
      if (c)
      {
        MPI_Waitall(c, &requests[0], &stat[0]);
      }
      delete[] block_lengths;

#else
      // throw error
      std::ostringstream error_message;
      error_message << "The preconditioner is distributed and on more than one "
                    << "processor. MPI is required.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    }
  }

  //============================================================================
  /// Given the naturally ordered vector, v, return
  /// the vector rearranged in block order in w. This function calls
  /// get_concatenated_block_vector(...) with the identity block mapping.
  ///
  /// This function has been re-written to work with the new dof type
  /// coarsening feature. The old function is kept alive in
  /// internal_get_block_ordered_preconditioner_vector(...) and is moved to
  /// the private section of the code. The differences between the two are:
  ///
  /// 1) This function extracts all the block vectors (in one go) via the
  ///    function internal_get_block_vectors(...), and concatenates them.
  ///
  /// 2) The old function makes use of the variables ending in "get_ordered",
  ///    thus is slightly more efficient since it does not have to concatenate
  ///    any block vectors.
  ///
  /// 3) The old function no longer respect the new indirections if dof types
  ///    have been coarsened.
  ///
  /// 4) This function extracts the most fine grain dof-level vectors and
  ///    concatenates them. These dof-level vectors respect the re-ordering
  ///    caused by the coarsening of dof types. The overhead associated with
  ///    concatenating DoubleVectors without communication is very small.
  ///
  /// This function should be used.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::get_block_ordered_preconditioner_vector(
    const DoubleVector& v, DoubleVector& w)
  {
#ifdef PARANOID
    if (!v.built())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*v.distribution_pt() != *this->master_distribution_pt())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must match the "
                    << " specified master_distribution_pt(). \n"
                    << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Get the number of blocks.
    unsigned nblocks = this->nblock_types();

    // Fill in the identity mapping.
    Vector<unsigned> block_vec_number(nblocks, 0);
    for (unsigned b = 0; b < nblocks; b++)
    {
      block_vec_number[b] = b;
    }

    // Do the work.
    get_concatenated_block_vector(block_vec_number, v, w);
  } // get_block_ordered_preconditioner_vector(...)

  //============================================================================
  /// Takes the block ordered vector, w, and reorders it in the natural
  /// order. Reordered vector is returned in v. Note: If the preconditioner is
  /// a subsidiary preconditioner then only the components of the vector
  /// associated with the blocks of the subsidiary preconditioner will be
  /// included. Hence the length of v is master_nrow() whereas that of the
  /// vector w is of length this->nrow().
  ///
  /// This is the return function for the function
  /// internal_get_block_ordered_preconditioner_vector(...).
  /// Both internal_get_block_ordered_preconditioner_vector(...) and
  /// internal_return_block_ordered_preconditioner_vector(...) has been
  /// superseded by the functions
  ///
  /// get_block_ordered_preconditioner_vector(...) and
  /// return_block_ordered_preconditioner_vector(...),
  ///
  /// Thus this function is moved to the private section of the code.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::
    internal_return_block_ordered_preconditioner_vector(const DoubleVector& w,
                                                        DoubleVector& v) const
  {
#ifdef PARANOID
    if (!v.built())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*v.distribution_pt() != *this->master_distribution_pt())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must match the "
                    << " specified master_distribution_pt(). \n"
                    << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (!w.built())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the block vector w must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*w.distribution_pt() !=
        *this->internal_preconditioner_matrix_distribution_pt())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the block vector w must match the "
                    << " specified distribution at Distribution_pt[b]";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // if + only one processor
    //    + more than one processor but matrix_pt is not distributed
    // then use the serial get_block method
    if (this->distribution_pt()->communicator_pt()->nproc() == 1 ||
        !this->distribution_pt()->distributed())
    {
      // number of blocks
      unsigned nblock = this->Internal_nblock_types;

      // copy to w
      unsigned block_offset = 0;
      const double* w_pt = w.values_pt();
      double* v_pt = v.values_pt();
      for (unsigned b = 0; b < nblock; b++)
      {
        unsigned block_nrow = this->internal_block_dimension(b);
        for (unsigned i = 0; i < block_nrow; i++)
        {
          v_pt[this->Global_index[b][i]] = w_pt[block_offset + i];
        }
        block_offset += block_nrow;
      }
    }
    // otherwise use mpi
    else
    {
#ifdef OOMPH_HAS_MPI

      // my rank
      unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

      // the number of processors
      unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

      // determine the maximum number of rows to be sent or recv
      unsigned max_n_send_or_recv = 0;
      for (unsigned p = 0; p < nproc; p++)
      {
        max_n_send_or_recv =
          std::max(max_n_send_or_recv, Nrows_to_send_for_get_ordered[p]);
        max_n_send_or_recv =
          std::max(max_n_send_or_recv, Nrows_to_recv_for_get_ordered[p]);
      }

      // create a vectors of 1s (the size of the nblock for the mpi indexed
      // data types
      int* block_lengths = new int[max_n_send_or_recv];
      for (unsigned i = 0; i < max_n_send_or_recv; i++)
      {
        block_lengths[i] = 1;
      }

      // perform the sends and receives
      Vector<MPI_Request> requests;
      for (unsigned p = 0; p < nproc; p++)
      {
        // send and recv with other processors
        if (p != my_rank)
        {
          if (Nrows_to_recv_for_get_ordered[p] > 0)
          {
            // create the send datatype
            MPI_Datatype type_send;
            MPI_Type_indexed(Nrows_to_recv_for_get_ordered[p],
                             block_lengths,
                             Rows_to_recv_for_get_ordered[p],
                             MPI_DOUBLE,
                             &type_send);
            MPI_Type_commit(&type_send);

            // send
            MPI_Request send_req;
            MPI_Isend(const_cast<double*>(w.values_pt()),
                      1,
                      type_send,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &send_req);
            MPI_Type_free(&type_send);
            requests.push_back(send_req);
          }

          if (Nrows_to_send_for_get_ordered[p] > 0)
          {
            // create the recv datatype
            MPI_Datatype type_recv;
            MPI_Type_indexed(Nrows_to_send_for_get_ordered[p],
                             block_lengths,
                             Rows_to_send_for_get_ordered[p],
                             MPI_DOUBLE,
                             &type_recv);
            MPI_Type_commit(&type_recv);

            // recv
            MPI_Request recv_req;
            MPI_Irecv(v.values_pt(),
                      1,
                      type_recv,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &recv_req);
            MPI_Type_free(&type_recv);
            requests.push_back(recv_req);
          }
        }

        // communicate wih self
        else
        {
          const double* w_values_pt = w.values_pt();
          double* v_values_pt = v.values_pt();
          for (unsigned i = 0; i < Nrows_to_send_for_get_ordered[p]; i++)
          {
            v_values_pt[Rows_to_send_for_get_ordered[p][i]] =
              w_values_pt[Rows_to_recv_for_get_ordered[p][i]];
          }
        }
      }

      // and then just wait
      unsigned c = requests.size();
      Vector<MPI_Status> stat(c);
      if (c)
      {
        MPI_Waitall(c, &requests[0], &stat[0]);
      }
      delete[] block_lengths;

#else
      // throw error
      std::ostringstream error_message;
      error_message << "The preconditioner is distributed and on more than one "
                    << "processor. MPI is required.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    } // else use mpi
  } // function return_block_ordered_preconditioner_vector


  //============================================================================
  /// Takes the block ordered vector, w, and reorders it in natural
  /// order. Reordered vector is returned in v. Note: If the preconditioner is
  /// a subsidiary preconditioner then only the components of the vector
  /// associated with the blocks of the subsidiary preconditioner will be
  /// included. Hence the length of v is master_nrow() whereas that of the
  /// vector w is of length this->nrow().
  ///
  /// This is the return function for the function
  /// get_block_ordered_preconditioner_vector(...).
  ///
  /// It calls the function return_concatenated_block_vector(...) with the
  /// identity block number ordering.
  //============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::return_block_ordered_preconditioner_vector(
    const DoubleVector& w, DoubleVector& v) const
  {
#ifdef PARANOID
    if (!v.built())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*v.distribution_pt() != *this->master_distribution_pt())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the global vector v must match the "
                    << " specified master_distribution_pt(). \n"
                    << "i.e. Distribution_pt in the master preconditioner";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (!w.built())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the block vector w must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*w.distribution_pt() != *this->preconditioner_matrix_distribution_pt())
    {
      std::ostringstream error_message;
      error_message << "The distribution of the block vector w must match the "
                    << "concatenations of distributions in "
                    << "Block_distribution_pt.\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Split into the block vectors.
    const unsigned nblocks = nblock_types();
    Vector<unsigned> block_vec_number(nblocks, 0);
    for (unsigned b = 0; b < nblocks; b++)
    {
      block_vec_number[b] = b;
    }

    return_concatenated_block_vector(block_vec_number, w, v);
  } // function return_block_ordered_preconditioner_vector

  //=============================================================================
  /// Gets block (i,j) from the matrix pointed to by
  /// Matrix_pt and returns it in output_block. This is associated with the
  /// internal blocks. Please use the other get_block(...) function.
  //=============================================================================
  template<>
  void BlockPreconditioner<CRDoubleMatrix>::internal_get_block(
    const unsigned& block_i,
    const unsigned& block_j,
    CRDoubleMatrix& output_block) const
  {
#ifdef PARANOID
    // the number of blocks
    const unsigned n_blocks = this->internal_nblock_types();

    // paranoid check that block i is in this block preconditioner
    if (block_i >= n_blocks || block_j >= n_blocks)
    {
      std::ostringstream error_message;
      error_message
        << "Requested block (" << block_i << "," << block_j
        << "), however this preconditioner has internal_nblock_types() "
        << "= " << internal_nblock_types() << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check that the matrix is the same as that of the master
    if (is_subsidiary_block_preconditioner())
    {
      if (master_block_preconditioner_pt()->matrix_pt() != matrix_pt())
      {
        std::string err = "Master and subs should have same matrix.";
        throw OomphLibError(
          err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Cast the pointer
    CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt());

    // if + only one processor
    //    + more than one processor but matrix_pt is not distributed
    // then use the serial get_block method
    if (cr_matrix_pt->distribution_pt()->communicator_pt()->nproc() == 1 ||
        !cr_matrix_pt->distribution_pt()->distributed())
    {
      // pointers for the jacobian matrix is compressed row sparse format
      int* j_row_start;
      int* j_column_index;
      double* j_value;

      // sets pointers to jacobian matrix
      j_row_start = cr_matrix_pt->row_start();
      j_column_index = cr_matrix_pt->column_index();
      j_value = cr_matrix_pt->value();

      // get the block dimensions
      unsigned block_nrow = this->internal_block_dimension(block_i);
      unsigned block_ncol = this->internal_block_dimension(block_j);

      // allocate temporary storage for the component vectors of block (i,j)
      // temp_ptr is used to point to an element in each column - required as
      // cannot assume that order of block's rows in jacobian and the block
      // matrix will be the same
      int* temp_row_start = new int[block_nrow + 1];
      for (unsigned i = 0; i <= block_nrow; i++)
      {
        temp_row_start[i] = 0;
      }
      Vector<int> temp_ptr(block_nrow + 1);
      int block_nnz = 0;

      // get number of rows in source matrix
      unsigned master_nrow = this->master_nrow();

      // determine how many non zeros there are in the block (i,j)
      // also determines how many non zeros are stored in each row or column -
      // stored in temp_ptr temporarily
      for (unsigned k = 0; k < master_nrow; k++)
      {
        if (internal_block_number(k) == static_cast<int>(block_i))
        {
          for (int l = j_row_start[k]; l < j_row_start[k + 1]; l++)
          {
            if (internal_block_number(j_column_index[l]) ==
                static_cast<int>(block_j))
            {
              block_nnz++;
              temp_ptr[internal_index_in_block(k) + 1]++;
            }
          }
        }
      }

      // if the matrix is not empty
      int* temp_column_index = new int[block_nnz];
      double* temp_value = new double[block_nnz];
      if (block_nnz > 0)
      {
        // uses number of elements in each column of block to determine values
        // for the block column start (temp_row_start)
        temp_row_start[0] = 0;
        for (unsigned k = 1; k <= block_nrow; k++)
        {
          temp_row_start[k] = temp_row_start[k - 1] + temp_ptr[k];
          temp_ptr[k] = temp_row_start[k];
        }

        // copies the relevant elements of the jacobian to the correct entries
        // of the block matrix
        for (unsigned k = 0; k < master_nrow; k++)
        {
          if (internal_block_number(k) == static_cast<int>(block_i))
          {
            for (int l = j_row_start[k]; l < j_row_start[k + 1]; l++)
            {
              if (internal_block_number(j_column_index[l]) ==
                  static_cast<int>(block_j))
              {
                int kk = temp_ptr[internal_index_in_block(k)]++;
                temp_value[kk] = j_value[l];
                temp_column_index[kk] =
                  internal_index_in_block(j_column_index[l]);
              }
            }
          }
        }
      }


      // Fill in the compressed row matrix ??ds Note: I kept the calls to
      // build as close as I could to before (had to replace new(dist) with
      // .build(dist) ).
      output_block.build(Internal_block_distribution_pt[block_i]);
      output_block.build_without_copy(
        block_ncol, block_nnz, temp_value, temp_column_index, temp_row_start);

#ifdef PARANOID
      // checks to see if block matrix has been set up correctly
      //   block_matrix_test(matrix_pt,block_i,block_j,block_pt);
      if (Run_block_matrix_test)
      {
        // checks to see if block matrix has been set up correctly
        block_matrix_test(block_i, block_j, &output_block);
      }
#endif
    }


    // otherwise we are dealing with a distributed matrix
    else
    {
#ifdef OOMPH_HAS_MPI
      // number of processors
      unsigned nproc = this->distribution_pt()->communicator_pt()->nproc();

      // my rank
      unsigned my_rank = this->distribution_pt()->communicator_pt()->my_rank();

      // sets pointers to jacobian matrix
      int* j_row_start = cr_matrix_pt->row_start();
      int* j_column_index = cr_matrix_pt->column_index();
      double* j_value = cr_matrix_pt->value();

      // number of non zeros in each row to be sent
      Vector<int*> nnz_send(nproc, 0);

      // number of non zeros in each row to be received
      Vector<int*> nnz_recv(nproc, 0);

      // storage for data to be sent
      Vector<int*> column_index_for_proc(nproc, 0);
      Vector<double*> values_for_proc(nproc, 0);

      // number of non zeros to be sent to each processor
      Vector<unsigned> total_nnz_send(nproc, 0);

      // number of rows of the block matrix on this processor
      unsigned nrow_local =
        Internal_block_distribution_pt[block_i]->nrow_local();

      // resize the nnz storage and compute nnz_send
      // and send and recv the nnz
      Vector<MPI_Request> send_req;
      Vector<MPI_Request> recv1_req;
      for (unsigned p = 0; p < nproc; p++)
      {
        int nrow_send = Nrows_to_send_for_get_block(block_i, p);
        int nrow_recv = Nrows_to_recv_for_get_block(block_i, p);

        // assemble nnz recv
        nnz_recv[p] = new int[nrow_recv];

        // assemble the storage to send
        if (nrow_send > 0 && p != my_rank)
        {
          nnz_send[p] = new int[nrow_send];
        }

        // compute the number of nnzs in each row and the total number
        // of nnzs
        for (int i = 0; i < nrow_send; i++)
        {
          unsigned row = Rows_to_send_for_get_block(block_i, p)[i];
          int c = 0;
          for (int r = j_row_start[row]; r < j_row_start[row + 1]; r++)
          {
            if (internal_block_number(j_column_index[r]) == int(block_j))
            {
              c++;
            }
          }
          if (p != my_rank)
          {
            nnz_send[p][i] = c;
          }
          else
          {
            nnz_recv[p][i] = c;
          }
          total_nnz_send[p] += c;
        }

        // send
        if (p != my_rank)
        {
          if (nrow_send)
          {
            MPI_Request req;
            MPI_Isend(nnz_send[p],
                      nrow_send,
                      MPI_INT,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &req);
            send_req.push_back(req);
          }

          // recv
          if (nrow_recv)
          {
            MPI_Request req;
            MPI_Irecv(nnz_recv[p],
                      nrow_recv,
                      MPI_INT,
                      p,
                      0,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &req);
            recv1_req.push_back(req);
          }
        }
      }

      // next assemble the values and row_start data to be sent for each
      // processor
      for (unsigned p = 0; p < nproc; p++)
      {
        int nrow_send = Nrows_to_send_for_get_block(block_i, p);

        // assemble the storage for the values and column indices to be sent
        if (p != my_rank)
        {
          if (total_nnz_send[p] > 0)
          {
            values_for_proc[p] = new double[total_nnz_send[p]];
            column_index_for_proc[p] = new int[total_nnz_send[p]];

            // copy the values and column indices to the storage
            unsigned ptr = 0;
            for (int i = 0; i < nrow_send; i++)
            {
              unsigned row = Rows_to_send_for_get_block(block_i, p)[i];
              for (int r = j_row_start[row]; r < j_row_start[row + 1]; r++)
              {
                if (internal_block_number(j_column_index[r]) == int(block_j))
                {
                  values_for_proc[p][ptr] = j_value[r];
                  column_index_for_proc[p][ptr] =
                    internal_index_in_block(j_column_index[r]);
                  ptr++;
                }
              }
            }

            // create the datatypes
            MPI_Datatype types[2];
            MPI_Type_contiguous(total_nnz_send[p], MPI_DOUBLE, &types[0]);
            MPI_Type_commit(&types[0]);
            MPI_Type_contiguous(total_nnz_send[p], MPI_INT, &types[1]);
            MPI_Type_commit(&types[1]);

            // get the start address of the vectors
            MPI_Aint displacement[2];
            MPI_Get_address(values_for_proc[p], &displacement[0]);
            MPI_Get_address(column_index_for_proc[p], &displacement[1]);

            // compute the displacements
            displacement[1] -= displacement[0];
            displacement[0] -= displacement[0];

            // compute the block lengths
            int length[2];
            length[0] = length[1] = 1;

            // build the struct data type
            MPI_Datatype final_type;
            MPI_Type_create_struct(2, length, displacement, types, &final_type);
            MPI_Type_commit(&final_type);
            MPI_Type_free(&types[0]);
            MPI_Type_free(&types[1]);

            // and send
            MPI_Request req;
            MPI_Isend(values_for_proc[p],
                      1,
                      final_type,
                      p,
                      1,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &req);
            send_req.push_back(req);
            MPI_Type_free(&final_type);
          }
        }
      }

      // wait for the recv to complete (the row_start recv which actually
      // contains the number of nnzs in each row)
      int c_recv = recv1_req.size();
      if (c_recv != 0)
      {
        MPI_Waitall(c_recv, &recv1_req[0], MPI_STATUS_IGNORE);
      }

      // compute the total number of nnzs to be received
      Vector<int> total_nnz_recv_from_proc(nproc);
      int local_block_nnz = 0;
      for (unsigned p = 0; p < nproc; p++)
      {
        // compute the total nnzs
        for (unsigned i = 0; i < Nrows_to_recv_for_get_block(block_i, p); i++)
        {
          total_nnz_recv_from_proc[p] += nnz_recv[p][i];
        }
        local_block_nnz += total_nnz_recv_from_proc[p];
      }

      // compute the offset for each block of nnzs (a matrix row) in the
      // values_recv and column_index_recv vectors

      // fisrt determine how many blocks of rows are to be recv
      Vector<int> n_recv_block(nproc, 0);
      for (unsigned p = 0; p < nproc; p++)
      {
        if (Nrows_to_recv_for_get_block(block_i, p) > 0)
        {
          n_recv_block[p] = 1;
        }
        for (unsigned i = 1; i < Nrows_to_recv_for_get_block(block_i, p); i++)
        {
          if (Rows_to_recv_for_get_block(block_i, p)[i] !=
              Rows_to_recv_for_get_block(block_i, p)[i - 1] + 1)
          {
            n_recv_block[p]++;
          }
        }
      }

      // next assemble row start recv
      int* row_start_recv = new int[nrow_local + 1];
      for (unsigned i = 0; i <= nrow_local; i++)
      {
        row_start_recv[i] = 0;
      }
      for (unsigned p = 0; p < nproc; p++)
      {
        for (unsigned i = 0; i < Nrows_to_recv_for_get_block(block_i, p); i++)
        {
          row_start_recv[Rows_to_recv_for_get_block(block_i, p)[i]] =
            nnz_recv[p][i];
        }
      }
      int g = row_start_recv[0];
      row_start_recv[0] = 0;
      for (unsigned i = 1; i < nrow_local; i++)
      {
        int temp_g = g;
        g = row_start_recv[i];
        row_start_recv[i] = row_start_recv[i - 1] + temp_g;
      }
      row_start_recv[nrow_local] = row_start_recv[nrow_local - 1] + g;

      // next assemble the offset and the number of nzs in each recv block
      Vector<int*> offset_recv_block(nproc, 0);
      Vector<int*> nnz_recv_block(nproc, 0);
      for (unsigned p = 0; p < nproc; p++)
      {
        if (Nrows_to_recv_for_get_block(block_i, p) > 0)
        {
          offset_recv_block[p] = new int[n_recv_block[p]];
          offset_recv_block[p][0] = 0;
          nnz_recv_block[p] = new int[n_recv_block[p]];
          for (int i = 0; i < n_recv_block[p]; i++)
          {
            nnz_recv_block[p][i] = 0;
          }
          unsigned ptr = 0;
          nnz_recv_block[p][ptr] += nnz_recv[p][0];
          offset_recv_block[p][0] =
            row_start_recv[Rows_to_recv_for_get_block(block_i, p)[0]];
          for (unsigned i = 1; i < Nrows_to_recv_for_get_block(block_i, p); i++)
          {
            if (Rows_to_recv_for_get_block(block_i, p)[i] !=
                Rows_to_recv_for_get_block(block_i, p)[i - 1] + 1)
            {
              ptr++;
              offset_recv_block[p][ptr] =
                row_start_recv[Rows_to_recv_for_get_block(block_i, p)[i]];
            }
            nnz_recv_block[p][ptr] += nnz_recv[p][i];
          }
        }
        delete[] nnz_recv[p];
      }

      // post the receives
      int* column_index_recv = new int[local_block_nnz];
      double* values_recv = new double[local_block_nnz];
      Vector<MPI_Request> recv2_req;
      for (unsigned p = 0; p < nproc; p++)
      {
        if (p != my_rank)
        {
          if (total_nnz_recv_from_proc[p] != 0)
          {
            // create the datatypes
            MPI_Datatype types[2];
            MPI_Type_indexed(n_recv_block[p],
                             nnz_recv_block[p],
                             offset_recv_block[p],
                             MPI_DOUBLE,
                             &types[0]);
            MPI_Type_commit(&types[0]);
            MPI_Type_indexed(n_recv_block[p],
                             nnz_recv_block[p],
                             offset_recv_block[p],
                             MPI_INT,
                             &types[1]);
            MPI_Type_commit(&types[1]);

            // compute the displacements
            MPI_Aint displacements[2];
            MPI_Get_address(values_recv, &displacements[0]);
            MPI_Get_address(column_index_recv, &displacements[1]);
            displacements[1] -= displacements[0];
            displacements[0] -= displacements[0];

            // compute the block lengths
            int length[2];
            length[0] = length[1] = 1;

            // create the final datatype
            MPI_Datatype final_type;
            MPI_Type_create_struct(
              2, length, displacements, types, &final_type);
            MPI_Type_commit(&final_type);
            MPI_Type_free(&types[0]);
            MPI_Type_free(&types[1]);

            // and the recv
            MPI_Request req;
            MPI_Irecv(values_recv,
                      1,
                      final_type,
                      p,
                      1,
                      this->distribution_pt()->communicator_pt()->mpi_comm(),
                      &req);
            recv2_req.push_back(req);
            MPI_Type_free(&final_type);
          }
        }
        else
        {
          // next send the values and column indices to self
          unsigned block_ptr = 0;
          unsigned counter = 0;
          int nrow_send = Nrows_to_send_for_get_block(block_i, my_rank);
          if (nrow_send > 0)
          {
            unsigned offset = offset_recv_block[my_rank][0];
            for (int i = 0; i < nrow_send; i++)
            {
              if (i > 0)
              {
                if (Rows_to_recv_for_get_block(block_i, p)[i] !=
                    Rows_to_recv_for_get_block(block_i, p)[i - 1] + 1)
                {
                  counter = 0;
                  block_ptr++;
                  offset = offset_recv_block[my_rank][block_ptr];
                }
              }
              unsigned row = Rows_to_send_for_get_block(block_i, my_rank)[i];
              for (int r = j_row_start[row]; r < j_row_start[row + 1]; r++)
              {
                if (internal_block_number(j_column_index[r]) == int(block_j))
                {
                  values_recv[offset + counter] = j_value[r];
                  column_index_recv[offset + counter] =
                    internal_index_in_block(j_column_index[r]);
                  counter++;
                }
              }
            }
          }
        }
      }

      // wait for the recv to complete (for the column_index and the values_
      c_recv = recv2_req.size();
      if (c_recv != 0)
      {
        MPI_Waitall(c_recv, &recv2_req[0], MPI_STATUS_IGNORE);
      }

      // Fill in the compressed row matrix
      output_block.build(Internal_block_distribution_pt[block_i]);
      output_block.build_without_copy(this->internal_block_dimension(block_j),
                                      local_block_nnz,
                                      values_recv,
                                      column_index_recv,
                                      row_start_recv);

      // wait for the send to complete (nnz / row_start)
      int c_send = send_req.size();
      if (c_send)
      {
        MPI_Waitall(c_send, &send_req[0], MPI_STATUS_IGNORE);
      }

      // delete temp storage used for assembling data for communication
      for (unsigned p = 0; p < nproc; p++)
      {
        delete[] nnz_send[p];
        delete[] column_index_for_proc[p];
        delete[] values_for_proc[p];
        delete[] offset_recv_block[p];
        delete[] nnz_recv_block[p];
      }
#else
      // throw error
      std::ostringstream error_message;
      error_message << "The matrix is distributed and on more than one "
                    << "processor. MPI is required.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    }
  }

  //=============================================================================
  /// Gets dof-level block (i,j).
  /// If Replacement_dof_block_pt(i,j) is not null, then the replacement
  /// block is returned via a deep copy.
  ///
  /// Otherwise if this is the uppermost block preconditioner then it calls
  /// internal_get_block(i,j), else if it is a subsidiary
  /// block preconditioner, it will call it's master block preconditioners'
  /// get_dof_level_block function.
  //=============================================================================
  template<>
  void BlockPreconditioner<CRDoubleMatrix>::get_dof_level_block(
    const unsigned& block_i,
    const unsigned& block_j,
    CRDoubleMatrix& output_block,
    const bool& ignore_replacement_block) const
  {
#ifdef PARANOID
    // the number of dof types.
    unsigned para_ndofs = ndof_types();

    // paranoid check that block i is in this block preconditioner
    if (block_i >= para_ndofs || block_j >= para_ndofs)
    {
      std::ostringstream err_msg;
      err_msg << "Requested dof block (" << block_i << "," << block_j
              << "), however this preconditioner has ndof_types() "
              << "= " << para_ndofs << std::endl;
      throw OomphLibError(
        err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    CRDoubleMatrix* tmp_block_pt =
      Replacement_dof_block_pt.get(block_i, block_j);

    if ((tmp_block_pt == 0) || ignore_replacement_block)
    {
      // Getting the block from parent preconditioner
      const unsigned ndof_in_parent_i =
        Doftype_coarsen_map_coarse[block_i].size();
      const unsigned ndof_in_parent_j =
        Doftype_coarsen_map_coarse[block_j].size();

      if (ndof_in_parent_i == 1 && ndof_in_parent_j == 1)
      {
        unsigned parent_dof_i = Doftype_coarsen_map_coarse[block_i][0];
        unsigned parent_dof_j = Doftype_coarsen_map_coarse[block_j][0];

        if (is_master_block_preconditioner())
        {
          internal_get_block(parent_dof_i, parent_dof_j, output_block);
        }
        else
        {
          parent_dof_i = Doftype_in_master_preconditioner_coarse[parent_dof_i];
          parent_dof_j = Doftype_in_master_preconditioner_coarse[parent_dof_j];

          master_block_preconditioner_pt()->get_dof_level_block(
            parent_dof_i, parent_dof_j, output_block, ignore_replacement_block);
        }
      }
      else
      {
        DenseMatrix<CRDoubleMatrix*> tmp_blocks_pt(
          ndof_in_parent_i, ndof_in_parent_j, 0);

        Vector<Vector<unsigned>> new_block(
          ndof_in_parent_i, Vector<unsigned>(ndof_in_parent_j, 0));

        for (unsigned dof_i = 0; dof_i < ndof_in_parent_i; dof_i++)
        {
          unsigned parent_dof_i = Doftype_coarsen_map_coarse[block_i][dof_i];
          if (is_subsidiary_block_preconditioner())
          {
            parent_dof_i =
              Doftype_in_master_preconditioner_coarse[parent_dof_i];
          }

          for (unsigned dof_j = 0; dof_j < ndof_in_parent_j; dof_j++)
          {
            unsigned parent_dof_j = Doftype_coarsen_map_coarse[block_j][dof_j];

            tmp_blocks_pt(dof_i, dof_j) = new CRDoubleMatrix;

            new_block[dof_i][dof_j] = 1;

            if (is_master_block_preconditioner())
            {
              internal_get_block(
                parent_dof_i, parent_dof_j, *tmp_blocks_pt(dof_i, dof_j));
            }
            else
            {
              parent_dof_j =
                Doftype_in_master_preconditioner_coarse[parent_dof_j];

              master_block_preconditioner_pt()->get_dof_level_block(
                parent_dof_i,
                parent_dof_j,
                *tmp_blocks_pt(dof_i, dof_j),
                ignore_replacement_block);
            }
          }
        }

        Vector<LinearAlgebraDistribution*> tmp_row_dist_pt(ndof_in_parent_i, 0);

        for (unsigned parent_dof_i = 0; parent_dof_i < ndof_in_parent_i;
             parent_dof_i++)
        {
          unsigned mapped_dof_i =
            Doftype_coarsen_map_coarse[block_i][parent_dof_i];

          if (is_master_block_preconditioner())
          {
            tmp_row_dist_pt[parent_dof_i] =
              Internal_block_distribution_pt[mapped_dof_i];
          }
          else
          {
            mapped_dof_i =
              Doftype_in_master_preconditioner_coarse[mapped_dof_i];

            tmp_row_dist_pt[parent_dof_i] =
              master_block_preconditioner_pt()->dof_block_distribution_pt(
                mapped_dof_i);
          }
        }

        Vector<LinearAlgebraDistribution*> tmp_col_dist_pt(ndof_in_parent_j, 0);

        for (unsigned parent_dof_j = 0; parent_dof_j < ndof_in_parent_j;
             parent_dof_j++)
        {
          unsigned mapped_dof_j =
            Doftype_coarsen_map_coarse[block_j][parent_dof_j];

          if (is_master_block_preconditioner())
          {
            tmp_col_dist_pt[parent_dof_j] =
              Internal_block_distribution_pt[mapped_dof_j];
          }
          else
          {
            mapped_dof_j =
              Doftype_in_master_preconditioner_coarse[mapped_dof_j];
            tmp_col_dist_pt[parent_dof_j] =
              master_block_preconditioner_pt()->dof_block_distribution_pt(
                mapped_dof_j);
          }
        }

        CRDoubleMatrixHelpers::concatenate_without_communication(
          tmp_row_dist_pt, tmp_col_dist_pt, tmp_blocks_pt, output_block);

        for (unsigned dof_i = 0; dof_i < ndof_in_parent_i; dof_i++)
        {
          for (unsigned dof_j = 0; dof_j < ndof_in_parent_j; dof_j++)
          {
            if (new_block[dof_i][dof_j])
            {
              delete tmp_blocks_pt(dof_i, dof_j);
            }
          }
        }
      }
    }
    else
    {
      CRDoubleMatrixHelpers::deep_copy(tmp_block_pt, output_block);
    }
  }

  //=============================================================================
  /// test function to check that every element in the block matrix
  /// (block_i,block_j) matches the corresponding element in the original matrix
  //=============================================================================
  template<typename MATRIX>
  void BlockPreconditioner<MATRIX>::block_matrix_test(
    const unsigned& block_i,
    const unsigned& block_j,
    const MATRIX* block_matrix_pt) const
  {
    // boolean flag to indicate whether test is passed
    bool check = true;

    // number of rows in matrix
    unsigned n_row = matrix_pt()->nrow();

    // number of columns in matrix
    unsigned n_col = matrix_pt()->ncol();

    // loop over rows of original matrix
    for (unsigned i = 0; i < n_row; i++)
    {
      // if this coefficient is associated with a block in this block
      // preconditioner
      if (static_cast<int>(block_i) == this->internal_block_number(i))
      {
        // loop over columns of original matrix
        for (unsigned j = 0; j < n_col; j++)
        {
          // if the coeeficient is associated with a block in this block
          // preconditioner
          if (static_cast<int>(block_j) == this->internal_block_number(j))
          {
            // check whether elements in original matrix and matrix of block
            // pointers match
            if (matrix_pt()->operator()(i, j) !=
                block_matrix_pt->operator()(internal_index_in_block(i),
                                            internal_index_in_block(j)))
            {
              check = false;
            }
          }
        }
      }
    }

    // throw error
    if (!check)
    {
      std::ostringstream error_message;
      error_message << "The require elements have not been successfully copied"
                    << " from the original matrix to the block matrices";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  template class BlockPreconditioner<CRDoubleMatrix>;

} // namespace oomph
