// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_LAGRANGE_ENFORCED_FLOW_PRECONDITIONERS_HEADER
#define OOMPH_LAGRANGE_ENFORCED_FLOW_PRECONDITIONERS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// oomphlib headers
#include "../generic/matrices.h"
#include "../generic/assembly_handler.h"
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
  //==========================================================================
  /// Namespace for subsidiary preconditioner creation helper functions
  //==========================================================================
  namespace Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper
  {
    /// CG with diagonal preconditioner for W-block subsidiary linear
    /// systems.
    extern Preconditioner* get_w_cg_preconditioner();

    /// Hypre Boomer AMG setting for the augmented momentum block
    /// of a 2D Navier-Stokes problem using the simple form of the viscous
    /// term (for serial code).
    extern Preconditioner* boomer_amg_for_2D_momentum_simple_visc();

    /// Hypre Boomer AMG setting for the augmented momentum block
    /// of a 2D Navier-Stokes problem using the stress divergence form of the
    /// viscous term (for serial code).
    extern Preconditioner* boomer_amg_for_2D_momentum_stressdiv_visc();

    /// Hypre Boomer AMG setting for the augmented momentum block
    /// of a 3D Navier-Stokes problem (for serial code).
    extern Preconditioner* boomer_amg_for_3D_momentum();

    /// Hypre Boomer AMG setting for the augmented momentum block
    /// of a 3D Navier-Stokes problem (for serial code).
    extern Preconditioner* boomer_amg2v22_for_3D_momentum();

    /// Hypre Boomer AMG setting for the 2D Poisson problem
    /// (for serial code).
    extern Preconditioner* boomer_amg_for_2D_poisson_problem();

    /// Hypre Boomer AMG setting for the 3D Poisson problem
    /// (for serial code).
    extern Preconditioner* boomer_amg_for_3D_poisson_problem();

  } // namespace
    // Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper


  //==========================================================================
  /// The preconditioner for the Lagrange multiplier constrained
  /// Navier-Stokes equations. The velocity components are constrained by
  /// Lagrange multiplier, which are applied via OOMPH-LIB's FACE elements.
  ///
  ///
  /// The linearised Jacobian takes the block form:
  ///
  /// | F_ns | L^T |
  /// |------------|
  /// |   L  | 0   |
  ///
  /// where L correspond to the constrained block,
  /// F_ns is the Navier-Stokes block with the following block structure
  ///
  /// |  F | G^T |
  /// |----------|
  /// |  D |  0  |
  ///
  /// Here F is the momentum block, G the discrete gradient operator,
  /// and D the discrete divergence operator. For unstabilised elements,
  /// we have D = G^T and in much of the literature the divergence matrix is
  /// denoted by B.
  ///
  /// The Lagrange enforced flow preconditioner takes the form:
  /// | F_aug |    |
  /// |-------|----|
  /// |       | Wd |
  ///
  /// where  F_aug = F_ns + L^T*inv(Wd)*L is an augmented Navier-Stokes block
  /// and Wd=(1/Scaling_sigma)*diag(LL^T).
  ///
  /// In our implementation of the preconditioner, the linear systems
  /// associated with the (1,1) block can either be solved "exactly",
  /// using SuperLU (in its incarnation as an exact preconditioner;
  /// this is the default) or by any other Preconditioner (inexact solver)
  /// specified via the access functions
  ///
  /// LagrangeEnforcedFlowPreconditioner::set_navier_stokes_preconditioner(...)
  ///
  /// Access to the elements is provided via meshes. However, a Vector of
  /// meshes is taken, each mesh contains a different type of block
  /// preconditionable element. This allows the (re-)classification of the
  /// constrained velocity DOF types.
  ///
  /// The first mesh in the Vector Mesh_pt must be the 'bulk' mesh.
  /// The rest are assumed to contain FACEELMENTS.
  ///
  /// Thus, the most general block structure (in 3D) is:
  ///
  ///  0 1 2 3   4 5 6 7  8  ..x   x+0 x+1 x+2 x+3 x+4
  /// [u v w p] [u v w l1 l2 ...] [u   v   w   l1  l2 ...] ...
  ///   Bulk       Surface 1             Surface 2         ...
  ///
  /// where the DOF types in [] are the DOF types associated with each mesh.
  ///
  /// For example, consider a unit cube domain [0,1]^3 with parallel outflow
  /// imposed (in mesh 0) and tangential flow imposed (in mesh 1), then there
  /// are 13 DOF types and our implementation respects the following
  /// (natural) DOF type order:
  ///
  ///   bulk          mesh 0           mesh 1
  /// [0 1 2 3] [4  5  6   7   8 ] [9  10 11 12 ]
  /// [u v w p] [up vp wp Lp1 Lp2] [ut vt wt Lt1]
  ///
  /// Via the appropriate mapping, the block_setup(...) function will
  /// re-order the DOF types into the following block types:
  ///
  ///  0 1  2  3 4  5   6 7  8     9   10  11  12   <- Block type
  ///  0 4  9  1 5  10  2 6  11    3    7   8  12   <- DOF type
  /// [u up ut v vp vt  w wp wt ] [p] [Lp1 Lp2 Lt1]
  ///
  //==========================================================================
  class LagrangeEnforcedFlowPreconditioner
    : public BlockPreconditioner<CRDoubleMatrix>
  {
  public:
    /// This preconditioner includes the option to use subsidiary
    /// operators other than SuperLUPreconditioner for this problem.
    /// This is the typedef of a function that should return an instance
    /// of a subsidiary preconditioning operator.  This preconditioner is
    /// responsible for the destruction of the subsidiary preconditioners.
    typedef Preconditioner* (*SubsidiaryPreconditionerFctPt)();

    /// Constructor - initialise variables.
    LagrangeEnforcedFlowPreconditioner() : BlockPreconditioner<CRDoubleMatrix>()
    {
      // The null pointer.
      Navier_stokes_preconditioner_pt = 0;

      // By default, the linear systems associated with the diagonal blocks
      // are solved "exactly" using SuperLU (in its incarnation as an exact
      // preconditioner. This is not a block preconditioner.
      Navier_stokes_preconditioner_is_block_preconditioner = false;

      // Flag to indicate to use SuperLU or not.
      Using_superlu_ns_preconditioner = true;

      // Empty vector of meshes and set the number of meshes to zero.
      My_mesh_pt.resize(0, 0);
      My_nmesh = 0;

      // The number of DOF types within the meshes.
      My_ndof_types_in_mesh.resize(0, 0);

      // Initialise other variables.
      Use_norm_f_for_scaling_sigma = true;
      Scaling_sigma = 0.0;
      N_lagrange_doftypes = 0;
      N_fluid_doftypes = 0;
      N_velocity_doftypes = 0;
      Preconditioner_has_been_setup = false;
    } // constructor

    /// Destructor
    virtual ~LagrangeEnforcedFlowPreconditioner()
    {
      this->clean_up_memory();
    }

    /// Broken copy constructor
    LagrangeEnforcedFlowPreconditioner(
      const LagrangeEnforcedFlowPreconditioner&) = delete;

    /// Broken assignment operator
    void operator=(const LagrangeEnforcedFlowPreconditioner&) = delete;

    /// Setup method for the LagrangeEnforcedFlowPreconditioner.
    void setup();

    /// Apply the preconditioner.
    /// r is the residual (rhs), z will contain the solution.
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    /// Set the meshes,
    /// the first mesh in the vector must be the bulk mesh.
    void set_meshes(const Vector<Mesh*>& mesh_pt);

    /// Set flag to use the infinite norm of the Navier-Stokes F matrix
    /// as the scaling sigma. This is the default behaviour. Note: the norm of
    /// the NS F matrix positive, however, we actually use the negative of
    /// the norm. This is because the underlying Navier-Stokes Jacobian is
    /// multiplied by -1. Ask Andrew/Matthias for more detail.
    void use_norm_f_for_scaling_sigma()
    {
      Use_norm_f_for_scaling_sigma = true;
    }

    /// Access function to set the scaling sigma.
    /// Note: this also sets the flag to use the infinite norm of
    /// the Navier-Stokes F matrix as the scaling sigma to false.
    /// Warning is given if trying to set scaling sigma to be equal to
    /// or greater than zero.
    void set_scaling_sigma(const double& scaling_sigma)
    {
      // Check if scaling sigma is zero or positive.
#ifdef PARANOID
      if (scaling_sigma == 0.0)
      {
        std::ostringstream warning_stream;
        warning_stream << "WARNING: \n"
                       << "Setting scaling_sigma = 0.0 may cause values.\n";
        OomphLibWarning(warning_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
      }
      if (scaling_sigma > 0.0)
      {
        std::ostringstream warning_stream;
        warning_stream << "WARNING: " << std::endl
                       << "The scaling (scaling_sigma) is positive: "
                       << Scaling_sigma << "\n"
                       << "Performance may be degraded.\n";
        OomphLibWarning(warning_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
      }
#endif

      Scaling_sigma = scaling_sigma;
      Use_norm_f_for_scaling_sigma = false;
    }

    /// Read (const) function to get the scaling sigma.
    double scaling_sigma() const
    {
      return Scaling_sigma;
    }

    /// Set a new Navier-Stokes matrix preconditioner
    /// (inexact solver)
    void set_navier_stokes_preconditioner(
      Preconditioner* new_ns_preconditioner_pt = 0);

    /// Set Navier-Stokes matrix preconditioner (inexact
    /// solver) to SuperLU
    void set_superlu_for_navier_stokes_preconditioner()
    {
      if (!Using_superlu_ns_preconditioner)
      {
        delete Navier_stokes_preconditioner_pt;
        Navier_stokes_preconditioner_pt = new SuperLUPreconditioner;
        Using_superlu_ns_preconditioner = true;
      }
    }

    /// Clears the memory.
    void clean_up_memory();

  private:
    /// Control flag is true if the preconditioner has been setup
    /// (used so we can wipe the data when the preconditioner is
    /// called again)
    bool Preconditioner_has_been_setup;

    /// Scaling for the augmentation: Scaling_sigma*(LL^T)
    double Scaling_sigma;

    /// Flag to indicate if we want to use the infinite norm of the
    /// Navier-Stokes momentum block for the scaling sigma.
    bool Use_norm_f_for_scaling_sigma;

    /// Inverse W values
    Vector<Vector<double>> Inv_w_diag_values;

    /// Pointer to the 'preconditioner' for the Navier-Stokes block
    Preconditioner* Navier_stokes_preconditioner_pt;

    /// Flag to indicate if the preconditioner for the Navier-Stokes
    /// block is a block preconditioner or not.
    bool Navier_stokes_preconditioner_is_block_preconditioner;

    /// Flag to indicate whether the default NS preconditioner is used
    bool Using_superlu_ns_preconditioner;

    /// Storage for the meshes. In our implementation, the first mesh
    /// must always be the Navier-Stokes (bulk) mesh, followed by surface
    /// meshes.
    Vector<Mesh*> My_mesh_pt;

    /// The number of DOF types in each mesh. This is used create
    /// various lookup lists.
    Vector<unsigned> My_ndof_types_in_mesh;

    /// The number of meshes. This is used to create various lookup
    /// lists.
    unsigned My_nmesh;

    /// The number of Lagrange multiplier DOF types.
    unsigned N_lagrange_doftypes;

    /// The number of fluid DOF types (including pressure).
    unsigned N_fluid_doftypes;

    /// The number of velocity DOF types.
    unsigned N_velocity_doftypes;

  }; // end of LagrangeEnforcedFlowPreconditioner class

} // namespace oomph
#endif
