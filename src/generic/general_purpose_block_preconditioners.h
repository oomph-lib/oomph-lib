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

// Include guards
#ifndef OOMPH_GENERAL_BLOCK_PRECONDITIONERS
#define OOMPH_GENERAL_BLOCK_PRECONDITIONERS


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// c++ include
// #include<list>

// oomph-lib includes
#include "matrices.h"
// #include "mesh.h"
// #include "problem.h"
#include "block_preconditioner.h"
#include "SuperLU_preconditioner.h"
#include "preconditioner_array.h"
#include "matrix_vector_product.h"


namespace oomph
{
  namespace PreconditionerCreationFunctions
  {
    /// Helper function to create a SuperLu preconditioner (for use as
    /// the default subsididary preconditioner creator in
    /// GeneralPurposeBlockPreconditioners).
    inline Preconditioner* create_super_lu_preconditioner()
    {
      return new SuperLUPreconditioner;
    }
  } // namespace PreconditionerCreationFunctions


  //============================================================================
  /// Base class for general purpose block preconditioners. Deals with
  /// setting subsidiary preconditioners and dof to block maps.
  /// Subsidiary preconditioners can be set in two ways:
  /// 1) A pointer to a subsidiary preconditioner for block i can be passed
  /// to set_subsidiary_preconditioner_pt(prec, i).
  /// 2) A default subsidiary preconditioner can be set up by providing a
  /// function pointer to a function which creates a preconditioner. During
  /// setup() all unset subsidiary preconditioner pointers will be filled in
  /// using this function. By default this uses SuperLU.
  //============================================================================
  template<typename MATRIX>
  class GeneralPurposeBlockPreconditioner : public BlockPreconditioner<MATRIX>
  {
  public:
    /// typedef for a function that allows other preconditioners to be
    /// employed to solve the subsidiary linear systems.
    /// The function should return a pointer to the required subsidiary
    /// preconditioner generated using new. This preconditioner is responsible
    /// for the destruction of the subsidiary preconditioners.
    typedef Preconditioner* (*SubsidiaryPreconditionerFctPt)();

    /// constructor
    GeneralPurposeBlockPreconditioner()
      : BlockPreconditioner<MATRIX>(),
        Subsidiary_preconditioner_creation_function_pt(
          &PreconditionerCreationFunctions::create_super_lu_preconditioner)
    {
      // Make sure that the Gp_mesh_pt container is size zero.
      Gp_mesh_pt.resize(0);
    }

    /// Destructor: clean up memory then delete all subsidiary
    /// preconditioners.
    virtual ~GeneralPurposeBlockPreconditioner()
    {
      this->clean_up_memory();

      for (unsigned j = 0, nj = Subsidiary_preconditioner_pt.size(); j < nj;
           j++)
      {
        delete Subsidiary_preconditioner_pt[j];
      }
    }

    /// ??ds I think clean_up_memory is supposed to clear out any stuff
    /// that doesn't need to be stored between solves. Call clean up on any
    /// non-null subsidiary preconditioners.
    virtual void clean_up_memory()
    {
      // Call clean up in any subsidiary precondtioners that are set.
      for (unsigned j = 0, nj = Subsidiary_preconditioner_pt.size(); j < nj;
           j++)
      {
        if (Subsidiary_preconditioner_pt[j] != 0)
        {
          Subsidiary_preconditioner_pt[j]->clean_up_memory();
        }
      }

      // Clean up the block preconditioner base class stuff
      this->clear_block_preconditioner_base();
    }

    /// Broken copy constructor
    GeneralPurposeBlockPreconditioner(
      const GeneralPurposeBlockPreconditioner&) = delete;

    /// Broken assignment operator
    void operator=(const GeneralPurposeBlockPreconditioner&) = delete;

    /// access function to set the subsidiary preconditioner function.
    void set_subsidiary_preconditioner_function(
      SubsidiaryPreconditionerFctPt sub_prec_fn)
    {
      Subsidiary_preconditioner_creation_function_pt = sub_prec_fn;
    }

    /// Reset the subsidiary preconditioner function to its default
    void reset_subsidiary_preconditioner_function_to_default()
    {
      Subsidiary_preconditioner_creation_function_pt =
        &PreconditionerCreationFunctions::create_super_lu_preconditioner;
    }
    /// Set the subsidiary preconditioner to use for block i. The
    /// subsidiary preconditioner should have been created using new (the
    /// general purpose block preconditioner will delete it later). If null
    /// the general purpose block preconditioner will use the
    /// Subsidiary_preconditioner_creation_function_pt to create the
    /// preconditioner during setup().
    void set_subsidiary_preconditioner_pt(Preconditioner* prec,
                                          const unsigned& i)
    {
      // If the vector is currently too small to hold that many
      // preconditioners then expand it and fill with nulls.
      if (Subsidiary_preconditioner_pt.size() < i + 1)
      {
        Subsidiary_preconditioner_pt.resize(i + 1, 0);
      }
      // Note: the size of the vector is checked by
      // fill_in_subsidiary_preconditioners(..)  when we know what size it
      // should be.

      // I'm assuming that the number of preconditioners is always "small"
      // compared to Jacobian size, so a resize doesn't waste much time.

      // Put the pointer in the vector
      Subsidiary_preconditioner_pt[i] = prec;
    }

    /// Get the subsidiary precondtioner pointer in block i (is
    /// allowed to be null if not yet set).
    Preconditioner* subsidiary_preconditioner_pt(const unsigned& i) const
    {
      return Subsidiary_preconditioner_pt[i];
    }

    /// Specify a DOF to block map
    void set_dof_to_block_map(Vector<unsigned>& dof_to_block_map)
    {
      Dof_to_block_map = dof_to_block_map;
    }

    /// Adds a mesh to be used by the
    /// block preconditioning framework for classifying DOF types. Optional
    /// boolean argument (default: false) allows the mesh to contain multiple
    /// element types.
    void add_mesh(const Mesh* mesh_pt,
                  const bool& allow_multiple_element_type_in_mesh = false)
    {
#ifdef PARANOID
      // Check that the mesh pointer is not null.
      if (mesh_pt == 0)
      {
        std::ostringstream err_msg;
        err_msg << "The mesh_pt is null, please point it to a mesh.\n";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // Push back the mesh pointer and the boolean in a pair.
      Gp_mesh_pt.push_back(
        std::make_pair(mesh_pt, allow_multiple_element_type_in_mesh));
    }

    /// Returns the number of meshes currently set in the
    /// GeneralPurposeBlockPreconditioner base class.
    unsigned gp_nmesh()
    {
      return Gp_mesh_pt.size();
    }

  protected:
    /// Set the mesh in the block preconditioning framework.
    void gp_preconditioner_set_all_meshes()
    {
      const unsigned nmesh = gp_nmesh();
#ifdef PARANOID
      if (nmesh == 0)
      {
        std::ostringstream err_msg;
        err_msg << "There are no meshes set.\n"
                << "Have you remembered to call add_mesh(...)?\n";
        throw OomphLibError(
          err_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      this->set_nmesh(nmesh);
      for (unsigned mesh_i = 0; mesh_i < nmesh; mesh_i++)
      {
        this->set_mesh(
          mesh_i, Gp_mesh_pt[mesh_i].first, Gp_mesh_pt[mesh_i].second);
      }
    }

    /// Modified block setup for general purpose block preconditioners
    void gp_preconditioner_block_setup()
    {
      if (Dof_to_block_map.size() > 0)
      {
        BlockPreconditioner<MATRIX>::block_setup(Dof_to_block_map);
      }
      else
      {
        BlockPreconditioner<MATRIX>::block_setup();
      }
    }

    /// Create any subsidiary preconditioners needed. Usually
    /// nprec_needed = nblock_types, except for the ExactBlockPreconditioner
    /// which only requires one preconditioner.
    void fill_in_subsidiary_preconditioners(const unsigned& nprec_needed)
    {
      // If it's empty then fill it in with null pointers.
      if (Subsidiary_preconditioner_pt.empty())
      {
        Subsidiary_preconditioner_pt.assign(nprec_needed, 0);
      }
      else
      {
        // Otherwise check we have the right number of them
#ifdef PARANOID
        if (Subsidiary_preconditioner_pt.size() != nprec_needed)
        {
          using namespace StringConversion;
          std::string error_msg = "Wrong number of precondtioners in";
          error_msg += "Subsidiary_preconditioner_pt, should have ";
          error_msg += to_string(nprec_needed) + " but we actually have ";
          error_msg += to_string(Subsidiary_preconditioner_pt.size());
          throw OomphLibError(
            error_msg, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }


      // Now replace any null pointers with new preconditioners
      for (unsigned j = 0, nj = Subsidiary_preconditioner_pt.size(); j < nj;
           j++)
      {
        if (Subsidiary_preconditioner_pt[j] == 0)
        {
          Subsidiary_preconditioner_pt[j] =
            (*Subsidiary_preconditioner_creation_function_pt)();
        }
      }
    }

    /// List of preconditioners to use for the blocks to be solved.
    Vector<Preconditioner*> Subsidiary_preconditioner_pt;

    /// Function to create any subsidiary preconditioners not set in
    /// Subsidiary_preconditioner_pt.
    SubsidiaryPreconditionerFctPt
      Subsidiary_preconditioner_creation_function_pt;

  private:
    /// the set of dof to block maps for this preconditioner
    Vector<unsigned> Dof_to_block_map;

    /// Vector of mesh pointers and a boolean indicating if we allow multiple
    /// element types in the same mesh.
    Vector<std::pair<const Mesh*, bool>> Gp_mesh_pt;
  };


  //=============================================================================
  /// Block diagonal preconditioner. By default SuperLU is used to solve
  /// the subsidiary systems, but other preconditioners can be used by setting
  /// them using passing a pointer to a function of type
  /// SubsidiaryPreconditionerFctPt to the method
  /// subsidiary_preconditioner_function_pt().
  //=============================================================================
  template<typename MATRIX>
  class BlockDiagonalPreconditioner
    : public GeneralPurposeBlockPreconditioner<MATRIX>
  {
  public:
    /// constructor - when the preconditioner is used as a master preconditioner
    BlockDiagonalPreconditioner() : GeneralPurposeBlockPreconditioner<MATRIX>()
    {
      // by default we do not use two level parallelism
      Use_two_level_parallelisation = false;

      // null the Preconditioner array pt
      Preconditioner_array_pt = 0;

      // Don't doc by default
      Doc_time_during_preconditioner_solve = false;
    }

    /// Destructor - delete the preconditioner matrices
    virtual ~BlockDiagonalPreconditioner()
    {
      this->clean_up_memory();
    }

    /// clean up the memory
    virtual void clean_up_memory()
    {
      if (Use_two_level_parallelisation)
      {
        delete Preconditioner_array_pt;
        Preconditioner_array_pt = 0;
      }

      // Clean up the base class too
      GeneralPurposeBlockPreconditioner<MATRIX>::clean_up_memory();
    }

    /// Broken copy constructor
    BlockDiagonalPreconditioner(const BlockDiagonalPreconditioner&) = delete;

    /// Broken assignment operator
    void operator=(const BlockDiagonalPreconditioner&) = delete;

    /// Apply preconditioner to r
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    /// Setup the preconditioner
    virtual void setup();

    /// Use two level parallelisation
    void enable_two_level_parallelisation()
    {
#ifndef OOMPH_HAS_MPI
      throw OomphLibError("Cannot do any parallelism since we don't have MPI.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
#endif
      Use_two_level_parallelisation = true;
    }

    /// Don't use two-level parallelisation
    void disable_two_level_parallelisation()
    {
      Use_two_level_parallelisation = false;
    }

    /// Enable Doc timings in application of block sub-preconditioners
    void enable_doc_time_during_preconditioner_solve()
    {
      Doc_time_during_preconditioner_solve = true;
    }

    /// Disable Doc timings in application of block sub-preconditioners
    void disable_doc_time_during_preconditioner_solve()
    {
      Doc_time_during_preconditioner_solve = false;
    }

    void fill_in_subsidiary_preconditioners(const unsigned& nprec_needed)
    {
#ifdef PARANOID
      if ((Use_two_level_parallelisation) &&
          !this->Subsidiary_preconditioner_pt.empty())
      {
        std::string err_msg =
          "Two level parallelism diagonal block preconditioners cannot have";
        err_msg +=
          " any preset preconditioners (due to weird memory management";
        err_msg += " in the PreconditionerArray, you could try fixing it).";
        throw OomphLibError(
          err_msg, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Now call the real function
      GeneralPurposeBlockPreconditioner<
        MATRIX>::fill_in_subsidiary_preconditioners(nprec_needed);
    }

  protected:
    /// This is a helper function to allow us to implement AntiDiagonal
    /// preconditioner by only changing this function. Get the second index
    /// for block number i. Obviously for a diagonal preconditioner we want
    /// the blocks (i,i), (for anti diagonal we will want blocks (i, nblock -
    /// i), see that class).
    virtual unsigned get_other_diag_ds(const unsigned& i,
                                       const unsigned& nblock) const
    {
      return i;
    }


  private:
    /// pointer for the PreconditionerArray
    PreconditionerArray* Preconditioner_array_pt;

    /// Use two level parallelism using the PreconditionerArray
    bool Use_two_level_parallelisation;

    /// Doc timings in application of block sub-preconditioners?
    bool Doc_time_during_preconditioner_solve;
  };


  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////


  //=============================================================================
  /// General purpose block triangular preconditioner
  /// By default this is Upper triangular.
  /// By default SuperLUPreconditioner (or SuperLUDistPreconditioner) is used to
  /// solve the subsidiary systems, but other preconditioners can be used by
  /// setting them using passing a pointer to a function of type
  /// SubsidiaryPreconditionerFctPt to the method
  /// subsidiary_preconditioner_function_pt().
  //=============================================================================
  template<typename MATRIX>
  class BlockTriangularPreconditioner
    : public GeneralPurposeBlockPreconditioner<MATRIX>
  {
  public:
    /// Constructor. (By default this preconditioner is upper triangular).
    BlockTriangularPreconditioner()
      : GeneralPurposeBlockPreconditioner<MATRIX>()
    {
      // default to upper triangular
      Upper_triangular = true;
    }

    /// Destructor - delete the preconditioner matrices
    virtual ~BlockTriangularPreconditioner()
    {
      this->clean_up_memory();
    }

    /// clean up the memory
    virtual void clean_up_memory()
    {
      // Delete anything in Off_diagonal_matrix_vector_products
      for (unsigned i = 0, ni = Off_diagonal_matrix_vector_products.nrow();
           i < ni;
           i++)
      {
        for (unsigned j = 0, nj = Off_diagonal_matrix_vector_products.ncol();
             j < nj;
             j++)
        {
          delete Off_diagonal_matrix_vector_products(i, j);
          Off_diagonal_matrix_vector_products(i, j) = 0;
        }
      }

      // Clean up the base class too
      GeneralPurposeBlockPreconditioner<MATRIX>::clean_up_memory();
    }

    /// Broken copy constructor
    BlockTriangularPreconditioner(const BlockTriangularPreconditioner&) =
      delete;

    /// Broken assignment operator
    void operator=(const BlockTriangularPreconditioner&) = delete;

    /// Apply preconditioner to r
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    /// Setup the preconditioner
    void setup();

    /// Use as an upper triangular preconditioner
    void upper_triangular()
    {
      Upper_triangular = true;
    }

    /// Use as a lower triangular preconditioner
    void lower_triangular()
    {
      Upper_triangular = false;
    }

  private:
    /// Matrix of matrix vector product operators for the off diagonals
    DenseMatrix<MatrixVectorProduct*> Off_diagonal_matrix_vector_products;

    /// Boolean indicating upper or lower triangular
    bool Upper_triangular;
  };


  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////


  //=============================================================================
  /// Exact block preconditioner - block preconditioner assembled from all
  /// blocks associated with the preconditioner and solved by SuperLU.
  //=============================================================================
  template<typename MATRIX>
  class ExactBlockPreconditioner
    : public GeneralPurposeBlockPreconditioner<MATRIX>
  {
  public:
    /// constructor
    ExactBlockPreconditioner() : GeneralPurposeBlockPreconditioner<MATRIX>() {}

    /// Destructor
    virtual ~ExactBlockPreconditioner() {}

    /// Broken copy constructor
    ExactBlockPreconditioner(const ExactBlockPreconditioner&) = delete;

    /// Broken assignment operator
    void operator=(const ExactBlockPreconditioner&) = delete;

    /// Apply preconditioner to r
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    /// Setup the preconditioner
    void setup();

    /// Access for the preconditioner pointer used to solve the
    /// system (stored in the vector of pointers in the base class);
    Preconditioner*& preconditioner_pt()
    {
      return this->Subsidiary_preconditioner_pt[0];
    }
  };


  // =================================================================
  /// Block "anti-diagonal" preconditioner, i.e. same as block
  /// diagonal but along the other diagonal of the matrix (top-right to
  /// bottom-left).
  // =================================================================
  template<typename MATRIX>
  class BlockAntiDiagonalPreconditioner
    : public BlockDiagonalPreconditioner<MATRIX>
  {
  protected:
    /// This is a helper function to allow us to implement BlockAntiDiagonal
    /// using BlockDiagonal preconditioner and only changing this
    /// function. Get the second index for block number i. Obviously for a
    /// diagonal preconditioner we want the blocks (i,i). For anti diagonal
    /// we will want blocks (i, nblock - i - 1).
    unsigned get_other_diag_ds(const unsigned& i, const unsigned& nblock) const
    {
      return nblock - i - 1;
    }
  };


  // =================================================================
  /// Preconditioner that doesn't actually do any preconditioning, it just
  /// allows access to the Jacobian blocks. This is pretty hacky but oh well..
  // =================================================================
  template<typename MATRIX>
  class DummyBlockPreconditioner
    : public GeneralPurposeBlockPreconditioner<MATRIX>
  {
  public:
    /// Constructor
    DummyBlockPreconditioner() : GeneralPurposeBlockPreconditioner<MATRIX>() {}

    /// Destructor
    ~DummyBlockPreconditioner() {}

    /// Broken copy constructor
    DummyBlockPreconditioner(const DummyBlockPreconditioner&) = delete;

    /// Broken assignment operator
    void operator=(const DummyBlockPreconditioner&) = delete;

    /// Apply preconditioner to r (just copy r to z).
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z)
    {
      z.build(r);
    }

    /// Setup the preconditioner
    void setup()
    {
      // Set up the block look up schemes
      this->block_setup();
    }
  };

} // namespace oomph
#endif
