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

#include "pseudo_elastic_fsi_preconditioner.h"

namespace oomph
{
  //=============================================================================
  /// clean up memory method
  //=============================================================================
  void PseudoElasticFSIPreconditioner::clean_up_memory()
  {
    // wipe the preconditioner
    Pseudo_elastic_preconditioner_pt->clean_up_memory();
    Navier_stokes_preconditioner_pt->clean_up_memory();
    Navier_stokes_schur_complement_preconditioner_pt->clean_up_memory();
    Solid_preconditioner_pt->clean_up_memory();

    // clean the subsidiary matvec operators
    Fluid_pseudo_elastic_matvec_pt->clean_up_memory();
    Solid_fluid_matvec_pt->clean_up_memory();
    Solid_pseudo_elastic_matvec_pt->clean_up_memory();
    Lagrange_solid_matvec_pt->clean_up_memory();
  }

  //=============================================================================
  /// \short Setup the precoonditioner.
  //=============================================================================
  void PseudoElasticFSIPreconditioner::setup()
  {
    // clean the memory
    this->clean_up_memory();

#ifdef PARANOID
    // paranoid check that the meshes have been set
    if (Fluid_and_pseudo_elastic_mesh_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "The fluid and pseudo elastic mesh must be set.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (Solid_mesh_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "The solid mesh must be set.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (Lagrange_multiplier_mesh_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "The Lagrange multiplier mesh must be set.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // add the meshes
    this->set_mesh(0, Fluid_and_pseudo_elastic_mesh_pt);
    this->set_mesh(1, Solid_mesh_pt);
    this->set_mesh(2, Lagrange_multiplier_mesh_pt);

    // determine the number of fluid dofs
    unsigned nfluid_dof = Dim + 1;

    // determine the number of pseudo solid dofs
    unsigned npseudo_elastic_dof = this->ndof_types_in_mesh(0) - nfluid_dof;

    // determine the number of solid dofs
    unsigned nsolid_dof = this->ndof_types_in_mesh(1);

    // determine the number of lagrange multiplier dofs
    unsigned nlagr_mult_dof = this->ndof_types_in_mesh(2);

    // total number of dof types
    unsigned ntotal_dof =
      nfluid_dof + npseudo_elastic_dof + nsolid_dof + nlagr_mult_dof;

    // setup the block lookup scheme
    // block 0 - fluid
    //       1 - solid
    //       2 - pseudo solid inc. lagrange mult
    Vector<unsigned> dof_to_block_map(ntotal_dof, 0);
    int c = nfluid_dof;
    for (unsigned i = 0; i < npseudo_elastic_dof; i++)
    {
      dof_to_block_map[c] = 2;
      c++;
    }
    for (unsigned i = 0; i < nsolid_dof; i++)
    {
      dof_to_block_map[c] = 1;
      c++;
    }
    for (unsigned i = 0; i < nlagr_mult_dof; i++)
    {
      dof_to_block_map[c] = 3;
      c++;
    }

    // Recast Jacobian matrix to CRDoubleMatrix
    CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt());
#ifdef PARANOID
    if (cr_matrix_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "FSIPreconditioner only works with"
                    << " CRDoubleMatrix matrices" << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Call block setup for this preconditioner
    this->block_setup(dof_to_block_map);

    // SETUP THE PRECONDITIONERS
    // =========================

    // setup the navier stokes preconditioner
    if (Use_navier_stokes_schur_complement_preconditioner)
    {
      Navier_stokes_schur_complement_preconditioner_pt->set_navier_stokes_mesh(
        Fluid_and_pseudo_elastic_mesh_pt);

      Vector<unsigned> ns_dof_list(nfluid_dof, 0);
      for (unsigned i = 0; i < nfluid_dof; i++)
      {
        ns_dof_list[i] = i;
      }

      Navier_stokes_schur_complement_preconditioner_pt
        ->turn_into_subsidiary_block_preconditioner(this, ns_dof_list);

      Navier_stokes_schur_complement_preconditioner_pt->Preconditioner::setup(
        matrix_pt());
    }
    else
    {
      CRDoubleMatrix* ns_matrix_pt = new CRDoubleMatrix;
      this->get_block(0, 0, *ns_matrix_pt);

      Navier_stokes_preconditioner_pt->setup(ns_matrix_pt);
      delete ns_matrix_pt;
      ns_matrix_pt = 0;
    }

    // next the solid preconditioner
    if (dynamic_cast<BlockPreconditioner<CRDoubleMatrix>*>(
          Solid_preconditioner_pt) != 0)
    {
      Solid_preconditioner_is_block_preconditioner = true;
      GeneralPurposeBlockPreconditioner<CRDoubleMatrix>*
        solid_block_preconditioner_pt =
          dynamic_cast<GeneralPurposeBlockPreconditioner<CRDoubleMatrix>*>(
            Solid_preconditioner_pt);

      if (solid_block_preconditioner_pt != 0)
      {
        unsigned offset = nfluid_dof + npseudo_elastic_dof;
        Vector<unsigned> solid_prec_dof_list(nsolid_dof);
        for (unsigned i = 0; i < nsolid_dof; i++)
        {
          solid_prec_dof_list[i] = offset + i;
        }
        solid_block_preconditioner_pt
          ->turn_into_subsidiary_block_preconditioner(this,
                                                      solid_prec_dof_list);
        solid_block_preconditioner_pt->setup(cr_matrix_pt);
      }
      else
      {
        std::ostringstream error_message;
        error_message << "If the (real) solid preconditioner is a "
                      << "BlockPreconditioner then is must be a "
                      << "GeneralPurposeBlockPreconditioner";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
    // otherwise it is a solid preconditioner
    else
    {
      Solid_preconditioner_is_block_preconditioner = false;
      CRDoubleMatrix* s_matrix_pt = new CRDoubleMatrix;
      this->get_block(1, 1, *s_matrix_pt);
      Solid_preconditioner_pt->setup(s_matrix_pt);
      delete s_matrix_pt;
      s_matrix_pt = 0;
    }

    // next the pseudo solid preconditioner
    unsigned ndof_for_pseudo_elastic_prec = Dim * 3;
    Vector<unsigned> pseudo_elastic_prec_dof_list(ndof_for_pseudo_elastic_prec,
                                                  0);
    for (unsigned i = 0; i < Dim * 2; i++)
    {
      pseudo_elastic_prec_dof_list[i] = nfluid_dof + i;
    }
    for (unsigned i = 0; i < Dim; i++)
    {
      pseudo_elastic_prec_dof_list[i + Dim * 2] =
        nfluid_dof + npseudo_elastic_dof + nsolid_dof + i;
    }
    Pseudo_elastic_preconditioner_pt->turn_into_subsidiary_block_preconditioner(
      this, pseudo_elastic_prec_dof_list);
    Pseudo_elastic_preconditioner_pt->set_elastic_mesh(
      this->Fluid_and_pseudo_elastic_mesh_pt);
    Pseudo_elastic_preconditioner_pt->set_lagrange_multiplier_mesh(
      this->Lagrange_multiplier_mesh_pt);
    Pseudo_elastic_preconditioner_pt->Preconditioner::setup(matrix_pt());

    // SETUP THE MATRIX VECTOR PRODUCT OPERATORS
    // =========================================

    // setup the fluid pseudo-solid matvec operator
    CRDoubleMatrix* fp_matrix_pt = new CRDoubleMatrix;
    get_block(0, 2, *fp_matrix_pt);
    //  Fluid_pseudo_elastic_matvec_pt->setup(fp_matrix_pt);
    this->setup_matrix_vector_product(
      Fluid_pseudo_elastic_matvec_pt, fp_matrix_pt, 2);
    delete fp_matrix_pt;
    fp_matrix_pt = 0;

    // setup the solid fluid matvec operator
    CRDoubleMatrix* sf_matrix_pt = new CRDoubleMatrix;
    get_block(1, 0, *sf_matrix_pt);
    //  Solid_fluid_matvec_pt->setup(sf_matrix_pt);
    this->setup_matrix_vector_product(Solid_fluid_matvec_pt, sf_matrix_pt, 0);
    delete sf_matrix_pt;
    sf_matrix_pt = 0;

    // setup the solid pseudo-solid matvec operator
    CRDoubleMatrix* sp_matrix_pt = new CRDoubleMatrix;
    get_block(1, 2, *sp_matrix_pt);
    //  Solid_pseudo_elastic_matvec_pt->setup(sp_matrix_pt);
    this->setup_matrix_vector_product(
      Solid_pseudo_elastic_matvec_pt, sp_matrix_pt, 2);
    delete sp_matrix_pt;
    sp_matrix_pt = 0;

    // build the lagrange solid matvec operator
    CRDoubleMatrix* ls_matrix_pt = new CRDoubleMatrix;
    get_block(3, 1, *ls_matrix_pt);
    //  Lagrange_solid_matvec_pt->setup(ls_matrix_pt);
    this->setup_matrix_vector_product(
      Lagrange_solid_matvec_pt, ls_matrix_pt, 1);
    delete ls_matrix_pt;
    ls_matrix_pt = 0;
  }

  //=============================================================================
  /// \short Apply the preconditioner
  //=============================================================================
  void PseudoElasticFSIPreconditioner::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
  {
    // apply the "pseudo solid" component of the pseudo solid preconditioner
    Pseudo_elastic_preconditioner_pt->elastic_preconditioner_solve(r, z);

    // apply the fluid on pseudo solid matrix vector product operator
    DoubleVector x;
    this->get_block_vector(2, z, x);
    DoubleVector y;
    Fluid_pseudo_elastic_matvec_pt->multiply(x, y);
    DoubleVector w;
    Solid_pseudo_elastic_matvec_pt->multiply(x, w);
    x.clear();
    this->get_block_vector(0, r, x);
    x -= y;
    y.clear();

    // storage for a copy of z
    DoubleVector z_copy;

    // apply the ns preconditioner
    if (Use_navier_stokes_schur_complement_preconditioner)
    {
      z_copy.build(z);
      this->return_block_vector(0, x, z_copy);
      x.clear();
      Navier_stokes_schur_complement_preconditioner_pt->preconditioner_solve(
        z_copy, z);
      z_copy.clear();
    }
    else
    {
      Navier_stokes_preconditioner_pt->preconditioner_solve(x, y);
      x.clear();
      this->return_block_vector(0, y, z);
      y.clear();
    }

    // apply the solid onto fluid matrix vector product operator
    this->get_block_vector(0, z, x);
    Solid_fluid_matvec_pt->multiply(x, y);
    x.clear();
    this->get_block_vector(1, r, x);
    x -= y;
    y.clear();

    // apply the result of the solid onto pseudo solid matrix vector product
    x -= w;
    w.clear();

    // apply the solid preconditioner
    if (Solid_preconditioner_is_block_preconditioner)
    {
      DoubleVector z_copy(z);
      this->return_block_vector(1, x, z_copy);
      x.clear();
      (dynamic_cast<GeneralPurposeBlockPreconditioner<CRDoubleMatrix>*>(
         Solid_preconditioner_pt))
        ->preconditioner_solve(z_copy, z);
      this->get_block_vector(1, z, y);
    }
    else
    {
      Solid_preconditioner_pt->preconditioner_solve(x, y);
      x.clear();
      this->return_block_vector(1, y, z);
    }

    // apply the lagrange multiplier solid matrix vector product operator
    Lagrange_solid_matvec_pt->multiply(y, x);
    y.clear();
    this->get_block_vector(3, r, y);
    y -= x;
    x.clear();
    z_copy.build(z);
    this->return_block_vector(3, y, z_copy);

    // apply the lagrange multiplier compenent of the pseudo solid
    // preconditioner
    Pseudo_elastic_preconditioner_pt->lagrange_multiplier_preconditioner_solve(
      z_copy, z);
    z_copy.clear();
  }
} // namespace oomph
