// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_PSEUDO_ELASTIC_FSI_PRECONDITIONER
#define OOMPH_PSEUDO_ELASTIC_FSI_PRECONDITIONER

// includes
#include "../generic/problem.h"
#include "../generic/block_preconditioner.h"
#include "../generic/preconditioner.h"
#include "../generic/SuperLU_preconditioner.h"
#include "../generic/matrix_vector_product.h"
#include "../navier_stokes/navier_stokes_preconditioners.h"
#include "../generic/general_purpose_block_preconditioners.h"
#include "pseudo_elastic_preconditioner.h"

namespace oomph
{
  //============================================================================
  /// Preconditioner for FSI problems with pseudo-elastic fluid node
  /// updates.
  /// Note:
  /// NavierStokesSchurComplementPreconditioner is applied to the Navier Stokes
  /// subsidiary system.
  /// Default solid preconditioner is SuperLUPreconditioner.
  /// \b Enumeration of Elastic DOF types in the Pseudo-Elastic Elements
  /// The method get_dof_types_for_unknowns() must be implemented such that
  /// DOFs subject be Lagrange multiplier and DOFs NOT subject to Lagrange
  /// multiplier have different labels. For example in a 3D problem there are
  /// 6 DOF types and the following labelling must be implemented:
  /// 0 - x displacement (without lagr mult traction)
  /// 1 - y displacement (without lagr mult traction)
  /// 2 - z displacement (without lagr mult traction)
  /// 3 - x displacement (with lagr mult traction)
  /// 4 - y displacement (with lagr mult traction)
  /// 5 - z displacement (with lagr mult traction)
  //============================================================================
  class PseudoElasticFSIPreconditioner
    : public BlockPreconditioner<CRDoubleMatrix>
  {
  public:
    /// constructor - just set defaults. Specify the spatial
    /// dimension of the fluid and a (non-const) problem pointer needed for
    /// the underlying NavierStokesSchurComplementPreconditioner.
    PseudoElasticFSIPreconditioner(const unsigned& dim, Problem* problem_pt)
      : Dim(dim)
    {
      Use_navier_stokes_schur_complement_preconditioner = true;

      // set the number of meshes
      this->set_nmesh(3);

      // null pointers
      Fluid_and_pseudo_elastic_mesh_pt = 0;
      Solid_mesh_pt = 0;
      Lagrange_multiplier_mesh_pt = 0;

      // create the pseudo solid preconditioner
      Pseudo_elastic_preconditioner_pt = new PseudoElasticPreconditioner();

      // using Schur complement preconditioner for NS
      Navier_stokes_preconditioner_pt = new SuperLUPreconditioner;
      Navier_stokes_schur_complement_preconditioner_pt =
        new NavierStokesSchurComplementPreconditioner(problem_pt);

      // set defaults
      Using_default_solid_preconditioner = true;

      // default super lu
      Solid_preconditioner_pt = new SuperLUPreconditioner;

      // create the matrix vector product operatrs
      Solid_fluid_matvec_pt = new MatrixVectorProduct;
      Solid_pseudo_elastic_matvec_pt = new MatrixVectorProduct;
      Fluid_pseudo_elastic_matvec_pt = new MatrixVectorProduct;
      Lagrange_solid_matvec_pt = new MatrixVectorProduct;
    }

    // destructor
    virtual ~PseudoElasticFSIPreconditioner()
    {
      // clean the memory
      this->clean_up_memory();

      // delete the pseudo solid preconditioner
      delete Pseudo_elastic_preconditioner_pt;

      // delete the navier stokes preconditioner
      delete Navier_stokes_preconditioner_pt;
      delete Navier_stokes_schur_complement_preconditioner_pt;

      // delete the solid preconditioner if default
      if (Using_default_solid_preconditioner)
      {
        delete Solid_preconditioner_pt;
      }

      // delete the matrix vector product operators
      delete Fluid_pseudo_elastic_matvec_pt;
      delete Solid_fluid_matvec_pt;
      delete Solid_pseudo_elastic_matvec_pt;
      delete Lagrange_solid_matvec_pt;
    }

    /// Broken copy constructor
    PseudoElasticFSIPreconditioner(const PseudoElasticFSIPreconditioner&) =
      delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const PseudoElasticFSIPreconditioner&) =
      delete;*/

    /// clean up memory method
    void clean_up_memory();

    /// Setup the precoonditioner.
    void setup();

    ///  Apply the preconditioner
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    /// specify the mesh containing the combined fluid/pseudo solid elements
    void set_fluid_and_pseudo_elastic_mesh_pt(Mesh* mesh_pt)
    {
      Fluid_and_pseudo_elastic_mesh_pt = mesh_pt;
    }

    /// specify the mesh containing the solid elements
    void set_solid_mesh_pt(Mesh* mesh_pt)
    {
      Solid_mesh_pt = mesh_pt;
    }

    /// specify the mesh containing the lagrange multiplier elements
    void set_lagrange_multiplier_mesh_pt(Mesh* mesh_pt)
    {
      Lagrange_multiplier_mesh_pt = mesh_pt;
    }

    /// speicify a non default solid preconditioner. This preconditioner
    /// will not delete it
    void set_solid_preconditioner(Preconditioner* prec_pt)
    {
      if (Using_default_solid_preconditioner)
      {
        delete Solid_preconditioner_pt;
      }
      Solid_preconditioner_pt = prec_pt;
      Using_default_solid_preconditioner = false;
    }

    /// Access function to the pseudo elastic subsidiary preconditioner
    PseudoElasticPreconditioner* const pseudo_elastic_preconditioner_pt()
    {
      return Pseudo_elastic_preconditioner_pt;
    }

    /// Access function to the Navier Stokes Schur complement preconditioner.
    NavierStokesSchurComplementPreconditioner* const navier_stokes_schur_complement_preconditioner_pt()
    {
      return Navier_stokes_schur_complement_preconditioner_pt;
    }

    /// Call to use the Navier Stokes Schur complement
    /// preconditioner.
    void enable_navier_stokes_schur_complement_preconditioner()
    {
      Use_navier_stokes_schur_complement_preconditioner = true;
    }

    /// Call to use the SuperLUPreconditioner is used for the
    /// Navier Stokes subsidiary system.
    void disable_navier_stokes_schur_complement_preconditioner()
    {
      Use_navier_stokes_schur_complement_preconditioner = false;
    }

  private:
    /// pointer to the pseudo solid preconditioner
    PseudoElasticPreconditioner* Pseudo_elastic_preconditioner_pt;

    /// pointer to the navier stokes precondtioner
    Preconditioner* Navier_stokes_preconditioner_pt;

    /// Navier Stokes Schur complement preconditioner.
    NavierStokesSchurComplementPreconditioner*
      Navier_stokes_schur_complement_preconditioner_pt;

    /// pointer to the solid preconditioner
    Preconditioner* Solid_preconditioner_pt;

    /// boolean flag to indicate whether default Solid preconditioner
    /// is used
    bool Using_default_solid_preconditioner;

    /// boolean flag to indicate whether the Solid preconditioner is a
    /// block preconditioner
    bool Solid_preconditioner_is_block_preconditioner;

    /// fluid onto pseudosolid matrix vector operator
    MatrixVectorProduct* Fluid_pseudo_elastic_matvec_pt;

    /// solid onto fluid matrix vector operatio
    MatrixVectorProduct* Solid_fluid_matvec_pt;

    /// solid onto pseudo solid matrix vector operatio
    MatrixVectorProduct* Solid_pseudo_elastic_matvec_pt;

    // lagrange onto solid matric vector product
    MatrixVectorProduct* Lagrange_solid_matvec_pt;

    /// Mesh containing the combined fluid and pseudo solid element
    Mesh* Fluid_and_pseudo_elastic_mesh_pt;

    /// Mesh containing the solid elements
    Mesh* Solid_mesh_pt;

    /// Mesh containing the lagrange multiplier elements
    Mesh* Lagrange_multiplier_mesh_pt;

    /// the dimension of the fluid
    unsigned Dim;

    /// If true the Navier Stokes Schur complement preconditioner
    /// is used. Otherwise SuperLUPreconditioner is used for the
    /// Navier Stokes subsidiary system.
    bool Use_navier_stokes_schur_complement_preconditioner;

  }; // end of class FSILagrangeMultiplierPreconditioner

} // namespace oomph
#endif
