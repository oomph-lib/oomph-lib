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
#ifndef OOMPH_PSEUDO_ELASTIC_SUBSIDIARY_PRECONDITIONER
#define OOMPH_PSEUDO_ELASTIC_SUBSIDIARY_PRECONDITIONER

// includes
#include "../generic/problem.h"
#include "../generic/block_preconditioner.h"
#include "../generic/preconditioner.h"
#include "../generic/SuperLU_preconditioner.h"
#include "../generic/matrix_vector_product.h"
#include "../generic/general_purpose_preconditioners.h"
#include "../generic/general_purpose_block_preconditioners.h"
#ifdef OOMPH_HAS_HYPRE
#include "../generic/hypre_solver.h"
#endif
#ifdef OOMPH_HAS_TRILINOS
#include "../generic/trilinos_solver.h"
#endif
namespace oomph
{
  //=============================================================================
  /// Functions to create instances of optimal subsidiary operators for
  /// the PseudoElasticPreconditioner. By default we use hypre for the
  /// the elastic blocks but can use Trilinos ML too.
  //=============================================================================
  namespace Pseudo_Elastic_Preconditioner_Subsidiary_Operator_Helper
  {
#ifdef OOMPH_HAS_HYPRE

    /// Hypre AMG w/ GS smoothing for the augmented elastic
    /// subsidiary linear systems
    Preconditioner* get_elastic_preconditioner_hypre();

    /// AMG w/ GS smoothing for the augmented elastic subsidiary linear
    /// systems -- calls Hypre version to stay consistent with previous default
    Preconditioner* get_elastic_preconditioner();

#endif

#ifdef OOMPH_HAS_TRILINOS

    /// TrilinosML smoothing for the augmented elastic
    /// subsidiary linear systems
    Preconditioner* get_elastic_preconditioner_trilinos_ml();

    /// CG with diagonal preconditioner for the lagrange multiplier
    /// subsidiary linear systems.
    Preconditioner* get_lagrange_multiplier_preconditioner();
#endif
  } // namespace Pseudo_Elastic_Preconditioner_Subsidiary_Operator_Helper


  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////


  //=============================================================================
  /// A subsidiary preconditioner for the pseudo-elastic FSI
  /// preconditioner. Also a stand-alone preconditioner for the problem of
  /// non-linear elasticity subject to prescribed displacement by Lagrange
  /// multiplier.
  /// \b Enumeration of Elastic DOF types in the Pseudo-Elastic Elements
  /// The method get_dof_types_for_unknowns() must be implemented such that
  /// DOFs subject be Lagrange multiplier and DOFs NOT subject to Lagrange
  /// multiplier have different labels. For example in a 3D problem there are
  /// 6 DOF types and the following labelling must be implemented:
  /// 0 - x displacement (without lagr mult traction)
  /// 1 - y displacement (without lagr mult traction)
  /// 2 - z displacement (without lagr mult traction)
  /// 4 - x displacement (with lagr mult traction)
  /// 5 - y displacement (with lagr mult traction)
  /// 6 - z displacement (with lagr mult traction)
  //=============================================================================
  class PseudoElasticPreconditioner : public BlockPreconditioner<CRDoubleMatrix>
  {
    /// PseudoElasticFSIPreconditioner is a friend to access the private
    /// *_preconditioner_solve(...) method
    friend class PseudoElasticFSIPreconditioner;

  public:
    /// This preconditioner includes the option to use subsidiary
    /// operators other than SuperLUPreconditioner for this problem.
    /// This is the typedef of a function that should return an instance
    /// of a subsidiary preconditioning operator.  This preconditioner is
    /// responsible for the destruction of the subsidiary preconditioners.
    typedef Preconditioner* (*SubsidiaryPreconditionerFctPt)();

    /// The augmented elasticity system can be preconditioned in one
    /// of four ways.
    /// 0 - Exact preconditioner
    /// 1 - Block diagonal preconditioning
    /// 2 - Block upper triangular preconditioner
    /// 3 - Block lower triangular preconditioner
    /// We group together the different components of the displacement vector
    /// field for the block decomposition.
    enum Elastic_preconditioner_type
    {
      Exact_block_preconditioner,
      Block_diagonal_preconditioner,
      Block_lower_triangular_preconditioner,
      Block_upper_triangular_preconditioner
    };

    /// Default (and only) constructor.
    PseudoElasticPreconditioner()
    {
      // null pointers
      Lagrange_multiplier_subsidiary_preconditioner_function_pt = 0;
      Elastic_subsidiary_preconditioner_function_pt = 0;
      Elastic_preconditioner_pt = 0;

      // set defaults
      Use_inf_norm_of_s_scaling = true;
      E_preconditioner_type = Exact_block_preconditioner;

      // resize the Mesh_pt
      Lagrange_multiplier_mesh_pt = 0;
      Elastic_mesh_pt = 0;
    }

    /// destructor
    virtual ~PseudoElasticPreconditioner()
    {
      this->clean_up_memory();
    }

    /// Broken copy constructor
    PseudoElasticPreconditioner(const PseudoElasticPreconditioner&) = delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=
     (const PseudoElasticPreconditioner&) = delete;*/

    /// Setup method for the PseudoElasticPreconditioner.
    void setup();

    /// Apply the preconditioner. Method implemented in two
    /// other methods (elastic and lagrange multiplier subsidiary
    /// preocnditioner) for the PseudoElasticFSIPreconditioner
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z)
    {
      this->elastic_preconditioner_solve(r, z);
      this->lagrange_multiplier_preconditioner_solve(r, z);
    }

    /// Access function to mesh containing the block-preconditionable
    /// elastic elements
    void set_elastic_mesh(Mesh* mesh_pt)
    {
      Elastic_mesh_pt = mesh_pt;
    }

    /// Access function to mesh containing the block-preconditionable
    /// lagrange multiplier elements
    void set_lagrange_multiplier_mesh(Mesh* mesh_pt)
    {
      Lagrange_multiplier_mesh_pt = mesh_pt;
    }

    /// Call to use the inf norm of S as scaling
    void enable_inf_norm_of_s_scaling()
    {
      Use_inf_norm_of_s_scaling = true;
    }

    /// Call to use no scaling
    void disable_inf_norm_of_s_scaling()
    {
      Use_inf_norm_of_s_scaling = false;
    }

    /// By default the Lagrange multiplier subsidiary systems are
    /// preconditioner with SuperLUPreconditioner. For a different
    /// preconditioner, pass a function to this
    /// method returning a different subsidiary operator.
    void set_lagrange_multiplier_subsidiary_preconditioner(
      SubsidiaryPreconditionerFctPt prec_fn)
    {
      Lagrange_multiplier_subsidiary_preconditioner_function_pt = prec_fn;
    }

    /// By default the elastic subsidiary systems are
    /// preconditioner with SuperLUPreconditioner. For a different
    /// preconditioner, pass a function to this
    /// method returning a different subsidiary operator.
    void set_elastic_subsidiary_preconditioner(
      SubsidiaryPreconditionerFctPt prec_fn)
    {
      Elastic_subsidiary_preconditioner_function_pt = prec_fn;
    }

    /// Set the type of preconditioner applied to the elastic:
    /// 0 - Exact preconditioner
    /// 1 - Block diagonal preconditioning
    /// 2 - Block upper triangular preconditioner
    /// 3 - Block lower triangular preconditioner
    /// We group together the different components of the displacement vector
    /// field for the block decomposition.
    Elastic_preconditioner_type& elastic_preconditioner_type()
    {
      return E_preconditioner_type;
    }

    /// Clears the memory.
    void clean_up_memory();

  private:
    /// Apply the elastic subsidiary preconditioner.
    void elastic_preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    ///  Apply the lagrange multiplier subsidiary preconditioner.
    void lagrange_multiplier_preconditioner_solve(const DoubleVector& r,
                                                  DoubleVector& z);

    /// The scaling. Defaults to infinity norm of S.
    double Scaling;

    /// boolean indicating whether the inf-norm of S should be used as
    /// scaling. Default = true;
    bool Use_inf_norm_of_s_scaling;

    /// An unsigned indicating which method should be used for
    /// preconditioning the solid component.
    Elastic_preconditioner_type E_preconditioner_type;

    /// the dimension of the problem
    unsigned Dim;

    /// storage for the preconditioner for the solid system
    Preconditioner* Elastic_preconditioner_pt;

    /// lagrange multiplier preconditioner pt
    Vector<Preconditioner*> Lagrange_multiplier_preconditioner_pt;

    /// The Lagrange multiplier subsidiary preconditioner function pointer
    SubsidiaryPreconditionerFctPt
      Lagrange_multiplier_subsidiary_preconditioner_function_pt;

    /// The solid subsidiary preconditioner function pointer
    SubsidiaryPreconditionerFctPt Elastic_subsidiary_preconditioner_function_pt;

    /// Pointer to the mesh containing the solid elements
    Mesh* Elastic_mesh_pt;

    /// Pointer to the mesh containing the Lagrange multiplier elements
    Mesh* Lagrange_multiplier_mesh_pt;
  };

  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////


  //=============================================================================
  /// A subsidiary preconditioner for the pseudo-elastic FSI
  /// preconditioner. Also a stand-alone preconditioner for the problem of
  /// non-linear elasticity subject to prescribed displacement by Lagrange
  /// multiplier..
  /// \b Enumeration of Elastic DOF types in the Pseudo-Elastic Elements
  /// The method get_dof_types_for_unknowns() must be implemented such that
  /// DOFs subject be Lagrange multiplier and DOFs NOT subject to Lagrange
  /// multiplier have different labels. For example in a 3D problem there are
  /// 6 DOF types and the following labelling must be implemented:
  /// 0 - x displacement (without lagr mult traction)
  /// 1 - y displacement (without lagr mult traction)
  /// 2 - z displacement (without lagr mult traction)
  /// 4 - x displacement (with lagr mult traction)
  /// 5 - y displacement (with lagr mult traction)
  /// 6 - z displacement (with lagr mult traction)
  //=============================================================================
  class PseudoElasticPreconditionerOld
    : public BlockPreconditioner<CRDoubleMatrix>
  {
    /// PseudoElasticFSIPreconditioner is a friend to access the private
    /// *_preconditioner_solve(...) method
    friend class PseudoElasticFSIPreconditioner;

  public:
    /// This preconditioner includes the option to use subsidiary
    /// operators other than SuperLUPreconditioner for this problem.
    /// This is the typedef of a function that should return an instance
    /// of a subsidiary preconditioning operator.  This preconditioner is
    /// responsible for the destruction of the subsidiary preconditioners.
    typedef Preconditioner* (*SubsidiaryPreconditionerFctPt)();

    /// The augmented elasticity system can be preconditioned in one
    /// of four ways.
    /// 0 - Exact preconditioner
    /// 1 - Block diagonal preconditioning
    /// 2 - Block upper triangular preconditioner
    /// 3 - Block lower triangular preconditioner
    /// We group together the different components of the displacement vector
    /// field for the block decomposition.
    enum Elastic_preconditioner_type
    {
      Exact_block_preconditioner,
      Block_diagonal_preconditioner,
      Block_lower_triangular_preconditioner,
      Block_upper_triangular_preconditioner
    };

    /// Default (and only) constructor.
    PseudoElasticPreconditionerOld()
    {
      // null pointers
      Lagrange_multiplier_subsidiary_preconditioner_function_pt = 0;
      Elastic_subsidiary_preconditioner_function_pt = 0;
      Elastic_preconditioner_pt = 0;

      // set defaults
      Use_inf_norm_of_s_scaling = true;
      E_preconditioner_type = Exact_block_preconditioner;

      // resize the Mesh_pt
      Lagrange_multiplier_mesh_pt = 0;
      Elastic_mesh_pt = 0;
    }

    /// destructor
    virtual ~PseudoElasticPreconditionerOld()
    {
      this->clean_up_memory();
    }

    /// Broken copy constructor
    PseudoElasticPreconditionerOld(const PseudoElasticPreconditionerOld&) =
      delete;

    /// Broken assignment operator
    /*void operator=(const PseudoElasticPreconditionerOld&) = delete;*/

    /// Setup method for the PseudoElasticPreconditionerOld.
    void setup();

    /// Apply the preconditioner. Method implemented in two
    /// other methods (elastic and lagrange multiplier subsidiary
    /// preocnditioner) for the PseudoElasticFSIPreconditioner
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z)
    {
      this->elastic_preconditioner_solve(r, z);
      this->lagrange_multiplier_preconditioner_solve(r, z);
    }

    /// Access function to mesh containing the block-preconditionable
    /// elastic elements
    void set_elastic_mesh(Mesh* mesh_pt)
    {
      Elastic_mesh_pt = mesh_pt;
    }

    /// Access function to mesh containing the block-preconditionable
    /// lagrange multiplier elements
    void set_lagrange_multiplier_mesh(Mesh* mesh_pt)
    {
      Lagrange_multiplier_mesh_pt = mesh_pt;
    }

    /// Call to use the inf norm of S as scaling
    void enable_inf_norm_of_s_scaling()
    {
      Use_inf_norm_of_s_scaling = true;
    }

    /// Call to use no scaling
    void disable_inf_norm_of_s_scaling()
    {
      Use_inf_norm_of_s_scaling = false;
    }

    /// By default the Lagrange multiplier subsidiary systems are
    /// preconditioner with SuperLUPreconditioner. For a different
    /// preconditioner, pass a function to this
    /// method returning a different subsidiary operator.
    void set_lagrange_multiplier_subsidiary_preconditioner(
      SubsidiaryPreconditionerFctPt prec_fn)
    {
      Lagrange_multiplier_subsidiary_preconditioner_function_pt = prec_fn;
    }

    /// By default the elastic subsidiary systems are
    /// preconditioner with SuperLUPreconditioner. For a different
    /// preconditioner, pass a function to this
    /// method returning a different subsidiary operator.
    void set_elastic_subsidiary_preconditioner(
      SubsidiaryPreconditionerFctPt prec_fn)
    {
      Elastic_subsidiary_preconditioner_function_pt = prec_fn;
    }

    /// Set the type of preconditioner applied to the elastic:
    /// 0 - Exact preconditioner
    /// 1 - Block diagonal preconditioning
    /// 2 - Block upper triangular preconditioner
    /// 3 - Block lower triangular preconditioner
    /// We group together the different components of the displacement vector
    /// field for the block decomposition.
    Elastic_preconditioner_type& elastic_preconditioner_type()
    {
      return E_preconditioner_type;
    }

    /// Clears the memory.
    void clean_up_memory();

  private:
    /// Apply the elastic subsidiary preconditioner.
    void elastic_preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    ///  Apply the lagrange multiplier subsidiary preconditioner.
    void lagrange_multiplier_preconditioner_solve(const DoubleVector& r,
                                                  DoubleVector& z);

    /// The scaling. Defaults to infinity norm of S.
    double Scaling;

    /// boolean indicating whether the inf-norm of S should be used as
    /// scaling. Default = true;
    bool Use_inf_norm_of_s_scaling;

    /// An unsigned indicating which method should be used for
    /// preconditioning the solid component.
    Elastic_preconditioner_type E_preconditioner_type;

    /// the dimension of the problem
    unsigned Dim;

    /// storage for the preconditioner for the solid system
    Preconditioner* Elastic_preconditioner_pt;

    /// lagrange multiplier preconditioner pt
    Vector<Preconditioner*> Lagrange_multiplier_preconditioner_pt;

    /// The Lagrange multiplier subsidary preconditioner function pointer
    SubsidiaryPreconditionerFctPt
      Lagrange_multiplier_subsidiary_preconditioner_function_pt;

    /// The solid subsidiary preconditioner function pointer
    SubsidiaryPreconditionerFctPt Elastic_subsidiary_preconditioner_function_pt;

    /// Pointer to the mesh containing the solid elements
    Mesh* Elastic_mesh_pt;

    /// Pointer to the mesh containing the Lagrange multiplier elements
    Mesh* Lagrange_multiplier_mesh_pt;
  };


  /// /////////////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////////////


  //=============================================================================
  /// Subsidiary helper preconditioner for the PseudoElasticPreconditioner.
  /// Required to construct the augmented elastic system prior to
  /// preconditioning.
  /// NOTE:
  /// 1. This is only intended to be used as a subsidiary preconditioner within
  /// the PseudoElasticPreconditioner.
  /// 2. If this preconditioner has N DOF types then the first N/2 are assumed
  /// to be ordinary solid DOF types, and the second N/2 are the solid DOF types
  /// with lagrange multiplier tractions applied.
  /// 3. By default this preconditioner uses a superlu preconditioner.
  //=============================================================================
  class PseudoElasticPreconditionerSubsidiaryPreconditionerOld
    : public BlockPreconditioner<CRDoubleMatrix>
  {
  public:
    /// typedef for a function that allows other preconditioners to be
    /// emplyed to solve the subsidiary linear systems.
    /// The function should return a pointer to the requred subsidiary
    /// preconditioner generated using new. This preconditioner is responsible
    /// for the destruction of the subsidiary preconditioners.
    typedef Preconditioner* (*SubsidiaryPreconditionerFctPt)();

    /// Constructor
    PseudoElasticPreconditionerSubsidiaryPreconditionerOld()
    {
      Scaling = 1.0;
      Preconditioner_pt = 0;
      Subsidiary_preconditioner_function_pt = 0;
    }

    /// Destructor
    ~PseudoElasticPreconditionerSubsidiaryPreconditionerOld()
    {
      this->clean_up_memory();
    }

    /// Broken copy constructor
    PseudoElasticPreconditionerSubsidiaryPreconditionerOld(
      const PseudoElasticPreconditionerSubsidiaryPreconditionerOld&) = delete;

    /// Broken assignment operator
    /*void operator=(const
     PseudoElasticPreconditionerSubsidiaryPreconditionerOld&) = delete;*/

    // Setup the preconditioner
    void setup();

    // Apply the preconditioner
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    /// Specify the scaling. Default is 1.0  Must be called before
    /// setup(...).
    double& scaling()
    {
      return Scaling;
    }

    /// access function to set the subsidiary preconditioner function.
    void set_subsidiary_preconditioner_function(
      SubsidiaryPreconditionerFctPt sub_prec_fn)
    {
      Subsidiary_preconditioner_function_pt = sub_prec_fn;
    }

  private:
    /// clears the memory
    void clean_up_memory()
    {
      delete Preconditioner_pt;
      Preconditioner_pt = 0;
    }

    // the augmentation scaling
    double Scaling;

    /// the preconditioner pt
    Preconditioner* Preconditioner_pt;

    /// the SubisidaryPreconditionerFctPt
    SubsidiaryPreconditionerFctPt Subsidiary_preconditioner_function_pt;
  }; // end of PseudoElasticPreconditionerSubsidiaryPreconditionerOld


  /// /////////////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////////////


  //=============================================================================
  /// Subsidiary helper preconditioner for the PseudoElasticPreconditioner.
  /// Required for block preconditioner of the augmented elastic subsidiary
  /// problem.
  /// NOTE:
  /// 1. This is only intended to be used as a subsidiary preconditioner within
  /// the PseudoElasticPreconditioner.
  /// 2. If this preconditioner has N DOF types then the first N/2 are assumed
  /// to be ordinary solid DOF types, and the second N/2 are the solid DOF types
  /// with lagrange multiplier tractions applied.
  /// 3. By default this preconditioner uses a superlu preconditioner.
  //=============================================================================
  class PseudoElasticPreconditionerSubsidiaryBlockPreconditionerOld
    : public BlockPreconditioner<CRDoubleMatrix>
  {
  public:
    /// This preconditioner includes the option to use subsidiary
    /// operators other than SuperLUPreconditioner for this problem.
    /// This is the typedef of a function that should return an instance
    /// of a subsidiary preconditioning operator.  This preconditioner is
    /// responsible for the destruction of the subsidiary preconditioners.
    typedef Preconditioner* (*SubsidiaryPreconditionerFctPt)();

    /// Constructor. (By default this preconditioner is upper triangular).
    PseudoElasticPreconditionerSubsidiaryBlockPreconditionerOld()
      : BlockPreconditioner<CRDoubleMatrix>()
    {
      // null the subsidiary preconditioner function pointer
      Subsidiary_preconditioner_function_pt = 0;

      // default to block diagonal
      Method = 0;

      // default scaling = 1.0
      Scaling = 1.0;
    };

    /// Destructor
    ~PseudoElasticPreconditionerSubsidiaryBlockPreconditionerOld()
    {
      this->clean_up_memory();
    }

    /// Broken copy constructor
    PseudoElasticPreconditionerSubsidiaryBlockPreconditionerOld(
      const PseudoElasticPreconditionerSubsidiaryBlockPreconditionerOld&) =
      delete;

    /// Broken assignment operator
    /*void operator=
     (const PseudoElasticPreconditionerSubsidiaryBlockPreconditionerOld&) =
     delete;*/

    /// clean up the memory
    void clean_up_memory();

    /// Setup the preconditioner
    void setup();

    /// Apply preconditioner to r
    void preconditioner_solve(const DoubleVector& res, DoubleVector& z);

    /// access function to set the subsidiary preconditioner function.
    void set_subsidiary_preconditioner_function(
      SubsidiaryPreconditionerFctPt sub_prec_fn)
    {
      Subsidiary_preconditioner_function_pt = sub_prec_fn;
    };

    /// use as a block diagonal preconditioner
    void use_block_diagonal_approximation()
    {
      Method = 0;
    }

    /// Use as an upper triangular preconditioner
    void use_upper_triangular_approximation()
    {
      Method = 1;
    }

    /// Use as a lower triangular preconditioner
    void use_lower_triangular_approximation()
    {
      Method = 2;
    }

    /// Specify the scaling. Default is 1.0  Must be set before
    /// setup(...).
    double& scaling()
    {
      return Scaling;
    }

  private:
    /// Vector of SuperLU preconditioner pointers for storing the
    /// preconditioners for each diagonal block
    Vector<PseudoElasticPreconditionerSubsidiaryPreconditionerOld*>
      Diagonal_block_preconditioner_pt;

    /// Matrix of matrix vector product operators for the off diagonals
    DenseMatrix<MatrixVectorProduct*> Off_diagonal_matrix_vector_products;

    /// the preconditioning method.
    /// 0 - block diagonal
    /// 1 - upper triangular
    /// 2 - lower triangular
    unsigned Method;

    /// The SubisidaryPreconditionerFctPt
    SubsidiaryPreconditionerFctPt Subsidiary_preconditioner_function_pt;

    /// The scaling. default 1.0.
    double Scaling;
  };


  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////


  // /*
  //=============================================================================
  /// A helper class for PseudoElasticPreconditioner.
  ///  Note that this is NOT actually a functioning preconditioner.
  /// We simply derive from this class to get access to the blocks.
  //=============================================================================
  class PseudoElasticPreconditionerScalingHelperOld
    : public BlockPreconditioner<CRDoubleMatrix>
  {
  public:
    /// The constructor.
    /// NOTE:
    /// 1. master_prec_pt should point to the
    /// PseudoElasticPreconditioner.
    /// 2. matrix_pt should point to the jacobian.
    /// 3. The vector dof_list should contain the full list of
    /// DOFS associated with the solid subsidiary system.
    /// 4. "solid_mesh_pt" should be a pointer to the solid mesh used in the
    ///   master preconditioner.
    PseudoElasticPreconditionerScalingHelperOld(
      BlockPreconditioner<CRDoubleMatrix>* master_prec_pt,
      CRDoubleMatrix* matrix_pt,
      Vector<unsigned>& dof_list,
      const Mesh* const solid_mesh_pt,
      const OomphCommunicator* comm_pt)
    {
      // turn into a subisiary preconditioner
      this->turn_into_subsidiary_block_preconditioner(master_prec_pt, dof_list);

      // all dofs are of the same block type
      Vector<unsigned> dof_to_block_map(dof_list.size(), 0);

      // store the matrix_pt
      set_matrix_pt(matrix_pt);

      // set the mesh
      this->set_nmesh(1);
      this->set_mesh(0, solid_mesh_pt);

      // set the communicator pointer
      this->set_comm_pt(comm_pt);

      // call block_setup(...)
      this->block_setup(dof_to_block_map);
    }

    /// Destructor.
    ~PseudoElasticPreconditionerScalingHelperOld()
    {
      this->clear_block_preconditioner_base();
    }

    /// Broken copy constructor
    PseudoElasticPreconditionerScalingHelperOld(
      const PseudoElasticPreconditionerScalingHelperOld&) = delete;

    /// Broken assignment operator
    /*void operator=(const PseudoElasticPreconditionerScalingHelperOld&) =
     * delete;*/

    /// returns the infinite norm of S
    double s_inf_norm()
    {
      CRDoubleMatrix* m_pt = new CRDoubleMatrix;
      this->get_block(0, 0, *m_pt);
      double s_inf_norm = m_pt->inf_norm();
      delete m_pt;
      return s_inf_norm;
    }

    // broken preconditioner setup
    void setup()
    {
      std::ostringstream error_message;
      error_message << "This method is intentionally broken. This class is not "
                       "a functioning "
                    << "preconditioner.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // broken preconditioner solve
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z)
    {
      std::ostringstream error_message;
      error_message << "This method is intentionally broken. This class is not "
                       "a functioning "
                    << "preconditioner.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


  }; // end of PseudoElasticPreconditionerScalingHelperOld
  // */

} // namespace oomph
#endif
