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
#include "solid_preconditioners.h"


namespace oomph
{
  //===========================================================================
  /// Setup the least-squares commutator solid preconditioner. This
  /// extracts blocks corresponding to the position/displacement and pressure
  /// unknowns, creates the matrices actually needed in the application of the
  /// preconditioner and deletes what can be deleted... Note that
  /// this preconditioner needs a CRDoubleMatrix.
  //============================================================================
  void PressureBasedSolidLSCPreconditioner::setup()
  {
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // NOTE: In the interest of minimising memory usage, several containers
    //       are recycled, therefore their content/meaning changes
    //       throughout this function. The code is carefully annotated
    //       but you'll have to read it line by line!
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // make sure any old data is deleted
    clean_up_memory();

#ifdef PARANOID
    // paranoid check that the solid mesh pt has been set
    if (Solid_mesh_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "The solid mesh pointer must be set.\n"
                    << "Use method set_solid_mesh(...)";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // set the mesh
    this->set_mesh(0, Solid_mesh_pt);

    // Get blocks
    // ----------

    // Set up block look up schemes (done automatically in the
    // BlockPreconditioner base class, based on the information
    // provided in the block-preconditionable elements in the problem)

    // this preconditioner has two types of block:
    // type 0: displacement/positions - corresponding to DOFs 0 to n-2
    // type 1: pressure - corresponding to DOF n-1
    double t_block_start = TimingHelpers::timer();
    unsigned ndof_types = 0;
    if (this->is_subsidiary_block_preconditioner())
    {
      ndof_types = this->ndof_types();
    }
    else
    {
      ndof_types = this->ndof_types_in_mesh(0);
    }
    Vector<unsigned> dof_to_block_map(ndof_types);
    dof_to_block_map[ndof_types - 1] = 1;
    this->block_setup(dof_to_block_map);
    double t_block_finish = TimingHelpers::timer();
    double block_setup_time = t_block_finish - t_block_start;
    if (Doc_time)
    {
      oomph_info << "Time for block_setup(...) [sec]: " << block_setup_time
                 << "\n";
    }

    // determine whether the F preconditioner is a block preconditioner (and
    // therefore a subsidiary preconditioner)
    BlockPreconditioner<CRDoubleMatrix>* F_block_preconditioner_pt =
      dynamic_cast<BlockPreconditioner<CRDoubleMatrix>*>(F_preconditioner_pt);
    F_preconditioner_is_block_preconditioner = true;
    if (F_block_preconditioner_pt == 0)
    {
      F_preconditioner_is_block_preconditioner = false;
    }

    // Get B (the divergence block)
    double t_get_B_start = TimingHelpers::timer();
    CRDoubleMatrix* b_pt = new CRDoubleMatrix;
    this->get_block(1, 0, *b_pt);
    double t_get_B_finish = TimingHelpers::timer();
    if (Doc_time)
    {
      double get_B_time = t_get_B_finish - t_get_B_start;
      oomph_info << "Time to get B [sec]: " << get_B_time << "\n";
    }

    // the pointer for f
    CRDoubleMatrix* f_pt = new CRDoubleMatrix;

    // the pointer for the p matrix
    CRDoubleMatrix* p_matrix_pt = new CRDoubleMatrix;

    // the pointer for bt
    CRDoubleMatrix* bt_pt = new CRDoubleMatrix;

    // if BFBt is to be formed
    if (Form_BFBt_product)
    {
      // If using scaling then replace B with Bt
      CRDoubleMatrix* ivmm_pt = 0;
      if (P_matrix_using_scaling)
      {
        // Assemble the ivmm diagonal
        double ivmm_assembly_start_t = TimingHelpers::timer();
        ivmm_pt = assemble_mass_matrix_diagonal();
        double ivmm_assembly_finish_t = TimingHelpers::timer();
        if (Doc_time)
        {
          double ivmm_assembly_time =
            ivmm_assembly_finish_t - ivmm_assembly_start_t;
          oomph_info << "Time to assemble inverse mass matrix [sec]: "
                     << ivmm_assembly_time << "\n";
        }

        // assemble BQ (stored in B)
        double t_BQ_start = TimingHelpers::timer();
        CRDoubleMatrix* temp_matrix_pt = new CRDoubleMatrix;
        b_pt->multiply(*ivmm_pt, *temp_matrix_pt);
        delete b_pt;
        b_pt = 0;
        b_pt = temp_matrix_pt;
        double t_BQ_finish = TimingHelpers::timer();
        if (Doc_time)
        {
          double t_BQ_time = t_BQ_finish - t_BQ_start;
          oomph_info << "Time to generate BQ [sec]: " << t_BQ_time << std::endl;
        }
      }

      // Get Bt
      double t_get_Bt_start = TimingHelpers::timer();
      this->get_block(0, 1, *bt_pt);
      double t_get_Bt_finish = TimingHelpers::timer();
      if (Doc_time)
      {
        double t_get_Bt_time = t_get_Bt_finish - t_get_Bt_start;
        oomph_info << "Time to get Bt [sec]: " << t_get_Bt_time << std::endl;
      }

      // now form the P matrix by multiplying B (which if using scaling will be
      // BQ) with Bt
      double t_P_start = TimingHelpers::timer();
      b_pt->multiply(*bt_pt, *p_matrix_pt);
      double t_P_finish = TimingHelpers::timer();
      if (Doc_time)
      {
        double t_P_time = t_P_finish - t_P_start;
        oomph_info << "Time to generate P matrix [sec]: " << t_P_time
                   << std::endl;
      }

      // Multiply auxiliary matrix by diagonal of mass matrix if
      // required
      if (P_matrix_using_scaling)
      {
        CRDoubleMatrix* temp_matrix_pt = new CRDoubleMatrix;
        double t_QBt_start = TimingHelpers::timer();
        ivmm_pt->multiply(*bt_pt, *temp_matrix_pt);
        delete bt_pt;
        bt_pt = 0;
        bt_pt = temp_matrix_pt;
        double t_QBt_finish = TimingHelpers::timer();

        // Output times
        if (Doc_time)
        {
          double t_QBt_time = t_QBt_finish - t_QBt_start;
          oomph_info << "Time to generate QBt [sec]: " << t_QBt_time
                     << std::endl;
        }
      }

      // Clean up memory
      delete ivmm_pt;

      // get block 0 0
      double t_get_F_start = TimingHelpers::timer();
      this->get_block(0, 0, *f_pt);
      double t_get_F_finish = TimingHelpers::timer();
      if (Doc_time)
      {
        double t_get_F_time = t_get_F_finish - t_get_F_start;
        oomph_info << "Time to get F [sec]: " << t_get_F_time << std::endl;
      }

      // Auxiliary matrix for intermediate results
      double t_aux_matrix_start = TimingHelpers::timer();
      CRDoubleMatrix* aux_matrix_pt = new CRDoubleMatrix;
      f_pt->multiply(*bt_pt, *aux_matrix_pt);
      double t_aux_matrix_finish = TimingHelpers::timer();
      if (Doc_time)
      {
        double t_aux_time = t_aux_matrix_finish - t_aux_matrix_start;
        oomph_info << "Time to generate FQBt [sec]: " << t_aux_time
                   << std::endl;
      }

      // can now delete Bt (or QBt if using scaling)
      delete bt_pt;
      bt_pt = 0;

      // and if F requires a block preconditioner then we can delete F
      if (F_preconditioner_is_block_preconditioner)
      {
        delete f_pt;
      }

      // now form BFBt
      double t_E_matrix_start = TimingHelpers::timer();
      CRDoubleMatrix* e_matrix_pt = new CRDoubleMatrix;
      b_pt->multiply(*aux_matrix_pt, *e_matrix_pt);
      delete aux_matrix_pt;
      delete b_pt;
      double t_E_matrix_finish = TimingHelpers::timer();
      if (Doc_time)
      {
        double t_E_time = t_E_matrix_finish - t_E_matrix_start;
        oomph_info << "Time to generate E (B*(F*Bt)) [sec]: " << t_E_time
                   << std::endl;
      }
      double t_E_matvec_start = TimingHelpers::timer();
      E_mat_vec_pt = new MatrixVectorProduct;
      // E_mat_vec_pt->setup(e_matrix_pt);
      this->setup_matrix_vector_product(E_mat_vec_pt, e_matrix_pt, 1);
      double t_E_matvec_finish = TimingHelpers::timer();
      delete e_matrix_pt;
      if (Doc_time)
      {
        double t_E_time = t_E_matvec_finish - t_E_matvec_start;
        oomph_info << "Time to build E (BFBt) matrix vector operator E [sec]: "
                   << t_E_time << std::endl;
      }

      // and rebuild Bt
      t_get_Bt_start = TimingHelpers::timer();
      bt_pt = new CRDoubleMatrix;
      this->get_block(0, 1, *bt_pt);
      t_get_Bt_finish = TimingHelpers::timer();
      if (Doc_time)
      {
        double t_get_Bt_time = t_get_Bt_finish - t_get_Bt_start;
        oomph_info << "Time to get Bt [sec]: " << t_get_Bt_time << std::endl;
      }
    }


    /// //////////////////////////////////////////////////////////////////////////


    // otherwise we are not forming BFBt
    else
    {
      // get the inverse mass matrix (Q)
      CRDoubleMatrix* ivmm_pt = 0;
      if (P_matrix_using_scaling)
      {
        double ivmm_assembly_start_t = TimingHelpers::timer();
        ivmm_pt = assemble_mass_matrix_diagonal();
        double ivmm_assembly_finish_t = TimingHelpers::timer();
        if (Doc_time)
        {
          double ivmm_assembly_time =
            ivmm_assembly_finish_t - ivmm_assembly_start_t;
          oomph_info << "Time to assemble Q (inverse diagonal "
                     << "mass matrix) [sec]: " << ivmm_assembly_time << "\n";
        }
      }

      // Get Bt
      double t_get_Bt_start = TimingHelpers::timer();
      this->get_block(0, 1, *bt_pt);
      double t_get_Bt_finish = TimingHelpers::timer();
      if (Doc_time)
      {
        double t_get_Bt_time = t_get_Bt_finish - t_get_Bt_start;
        oomph_info << "Time to get Bt [sec]: " << t_get_Bt_time << std::endl;
      }

      // next QBt
      if (P_matrix_using_scaling)
      {
        double t_QBt_matrix_start = TimingHelpers::timer();
        CRDoubleMatrix* qbt_pt = new CRDoubleMatrix;
        ivmm_pt->multiply(*bt_pt, *qbt_pt);
        delete bt_pt;
        bt_pt = 0;
        bt_pt = qbt_pt;
        double t_QBt_matrix_finish = TimingHelpers::timer();
        if (Doc_time)
        {
          double t_QBt_time = t_QBt_matrix_finish - t_QBt_matrix_start;
          oomph_info << "Time to generate QBt [sec]: " << t_QBt_time
                     << std::endl;
        }
        delete ivmm_pt;
      }

      // form P
      double t_p_matrix_start = TimingHelpers::timer();
      b_pt->multiply(*bt_pt, *p_matrix_pt);
      double t_p_matrix_finish = TimingHelpers::timer();
      if (Doc_time)
      {
        double t_p_time = t_p_matrix_finish - t_p_matrix_start;
        oomph_info << "Time to generate P [sec]: " << t_p_time << std::endl;
      }
      delete b_pt;
      b_pt = 0;

      // build the matvec operator for QBt
      double t_QBt_MV_start = TimingHelpers::timer();
      QBt_mat_vec_pt = new MatrixVectorProduct;
      // QBt_mat_vec_pt->setup(bt_pt);
      this->setup_matrix_vector_product(QBt_mat_vec_pt, bt_pt, 1);
      double t_QBt_MV_finish = TimingHelpers::timer();
      if (Doc_time)
      {
        double t_p_time = t_QBt_MV_finish - t_QBt_MV_start;
        oomph_info << "Time to build QBt matrix vector operator [sec]: "
                   << t_p_time << std::endl;
      }
      delete bt_pt;
      bt_pt = 0;

      // get F
      double t_get_F_start = TimingHelpers::timer();
      this->get_block(0, 0, *f_pt);
      double t_get_F_finish = TimingHelpers::timer();
      if (Doc_time)
      {
        double t_get_F_time = t_get_F_finish - t_get_F_start;
        oomph_info << "Time to get F [sec]: " << t_get_F_time << std::endl;
      }

      // form the matrix vector product helper
      double t_F_MV_start = TimingHelpers::timer();
      F_mat_vec_pt = new MatrixVectorProduct;
      // F_mat_vec_pt->setup(f_pt);
      this->setup_matrix_vector_product(F_mat_vec_pt, f_pt, 0);
      double t_F_MV_finish = TimingHelpers::timer();
      if (Doc_time)
      {
        double t_F_MV_time = t_F_MV_finish - t_F_MV_start;
        oomph_info << "Time to build F Matrix Vector Operator [sec]: "
                   << t_F_MV_time << std::endl;
      }

      // if F is a block preconditioner then we can delete the F matrix
      if (F_preconditioner_is_block_preconditioner)
      {
        delete f_pt;
        f_pt = 0;
      }

      // and rebuild Bt
      t_get_Bt_start = TimingHelpers::timer();
      bt_pt = new CRDoubleMatrix;
      this->get_block(0, 1, *bt_pt);
      t_get_Bt_finish = TimingHelpers::timer();
      if (Doc_time)
      {
        double t_get_Bt_time = t_get_Bt_finish - t_get_Bt_start;
        oomph_info << "Time to get Bt [sec]: " << t_get_Bt_time << std::endl;
      }
    }


    /// //////////////////////////////////////////////////////////////////////////


    // form the matrix vector operator for Bt
    double t_Bt_MV_start = TimingHelpers::timer();
    Bt_mat_vec_pt = new MatrixVectorProduct;
    // Bt_mat_vec_pt->setup(bt_pt);
    this->setup_matrix_vector_product(Bt_mat_vec_pt, bt_pt, 1);
    double t_Bt_MV_finish = TimingHelpers::timer();
    if (Doc_time)
    {
      double t_Bt_MV_time = t_Bt_MV_finish - t_Bt_MV_start;
      oomph_info << "Time to build Bt Matrix Vector Operator [sec]: "
                 << t_Bt_MV_time << std::endl;
    }
    delete bt_pt;
    bt_pt = 0;

    // if the P preconditioner has not been setup
    if (P_preconditioner_pt == 0)
    {
      P_preconditioner_pt = new SuperLUPreconditioner;
      Using_default_p_preconditioner = true;
    }

    // std::stringstream junk;
    // junk << "p_matrix" << comm_pt()->my_rank()
    //      << ".dat";
    // p_matrix_pt->sparse_indexed_output_with_offset(junk.str());
    // oomph_info << "Done output of " << junk.str() << std::endl;

    // Setup the preconditioner for the Pressure matrix
    double t_p_prec_start = TimingHelpers::timer();
    P_preconditioner_pt->setup(p_matrix_pt);
    delete p_matrix_pt;
    p_matrix_pt = 0;
    p_matrix_pt = 0;
    double t_p_prec_finish = TimingHelpers::timer();
    if (Doc_time)
    {
      double t_p_prec_time = t_p_prec_finish - t_p_prec_start;
      oomph_info << "P sub-preconditioner setup time [sec]: " << t_p_prec_time
                 << "\n";
    }

    // Set up solver for solution of system with momentum matrix
    // ----------------------------------------------------------

    // if the F preconditioner has not been setup
    if (F_preconditioner_pt == 0)
    {
      F_preconditioner_pt = new SuperLUPreconditioner;
      Using_default_f_preconditioner = true;
    }

    // if F is a block preconditioner
    double t_f_prec_start = TimingHelpers::timer();
    if (F_preconditioner_is_block_preconditioner)
    {
      unsigned ndof_types = this->ndof_types();
      ndof_types--;
      Vector<unsigned> dof_map(ndof_types);
      for (unsigned i = 0; i < ndof_types; i++)
      {
        dof_map[i] = i;
      }
      F_block_preconditioner_pt->turn_into_subsidiary_block_preconditioner(
        this, dof_map);
      F_block_preconditioner_pt->setup(matrix_pt());
    }
    // otherwise F is not a block preconditioner
    else
    {
      F_preconditioner_pt->setup(f_pt);
      delete f_pt;
      f_pt = 0;
    }
    double t_f_prec_finish = TimingHelpers::timer();
    if (Doc_time)
    {
      double t_f_prec_time = t_f_prec_finish - t_f_prec_start;
      oomph_info << "F sub-preconditioner setup time [sec]: " << t_f_prec_time
                 << "\n";
    }

    // Remember that the preconditioner has been setup so
    // the stored information can be wiped when we
    // come here next...
    Preconditioner_has_been_setup = true;
  }


  //=======================================================================
  /// Apply preconditioner to r.
  //=======================================================================
  void PressureBasedSolidLSCPreconditioner::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
  {
#ifdef PARANOID
    if (Preconditioner_has_been_setup == false)
    {
      std::ostringstream error_message;
      error_message << "setup must be called before using preconditioner_solve";
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

    double t_start_overall = TimingHelpers::timer();
    double t_start = TimingHelpers::timer();
    double t_end = 0;

    // if z is not setup then give it the same distribution
    if (!z.built())
    {
      z.build(r.distribution_pt(), 0.0);
    }

    // Step 1 - apply approximate Schur inverse to pressure unknowns (block 1)
    // -----------------------------------------------------------------------

    // Working vectors
    DoubleVector temp_vec;
    DoubleVector another_temp_vec;

    // Copy pressure values from residual vector to temp_vec:
    // Loop over all entries in the global vector (this one
    // includes displacement/position and pressure dofs in some random fashion)
    this->get_block_vector(1, r, temp_vec);


    if (Doc_time)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "LSC prec solve: Time for get block vector: "
                 << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }

    // NOTE: The vector temp_vec now contains the vector r_p.

    // Solve first pressure Poisson system
#ifdef PARANOID
    // check a solver has been set
    if (P_preconditioner_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "P_preconditioner_pt has not been set.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // use some Preconditioner's preconditioner_solve function
    P_preconditioner_pt->preconditioner_solve(temp_vec, another_temp_vec);


    if (Doc_time)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "LSC prec solve: First P solve [nrow="
                 << P_preconditioner_pt->nrow() << "]: " << t_end - t_start
                 << std::endl;
      t_start = TimingHelpers::timer();
    }


    // NOTE: The vector another_temp_vec now contains the vector P^{-1} r_p

    // Multiply another_temp_vec by matrix E and stick the result into temp_vec
    temp_vec.clear();
    if (Form_BFBt_product)
    {
      E_mat_vec_pt->multiply(another_temp_vec, temp_vec);
    }
    else
    {
      QBt_mat_vec_pt->multiply(another_temp_vec, temp_vec);
      another_temp_vec.clear();
      F_mat_vec_pt->multiply(temp_vec, another_temp_vec);
      temp_vec.clear();
      QBt_mat_vec_pt->multiply_transpose(another_temp_vec, temp_vec);
    }


    if (Doc_time)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "LSC prec solve: E matrix vector product: "
                 << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }

    // NOTE: The vector temp_vec now contains E P^{-1} r_p

    // Solve second pressure Poisson system using preconditioner_solve
    another_temp_vec.clear();
    P_preconditioner_pt->preconditioner_solve(temp_vec, another_temp_vec);


    if (Doc_time)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "LSC prec solve: Second P solve [nrow="
                 << P_preconditioner_pt->nrow() << "]: " << t_end - t_start
                 << std::endl;
      t_start = TimingHelpers::timer();
    }


    // NOTE: The vector another_temp_vec now contains  z_p = P^{-1} E P^{-1} r_p
    //       as required (apart from the sign which we'll fix in the
    //       next step.

    // Now copy another_temp_vec (i.e. z_p) back into the global vector z.
    // Loop over all entries in the global results vector z:
    temp_vec.build(another_temp_vec.distribution_pt(), 0.0);
    temp_vec -= another_temp_vec;
    return_block_vector(1, temp_vec, z);

    // Step 2 - apply preconditioner to displacement/positon unknowns (block 0)
    // ------------------------------------------------------------------------

    // Recall that another_temp_vec (computed above) contains the
    // negative of the solution of the Schur complement systen, -z_p.
    // Multiply by G (stored in Block_matrix_pt(0,1) and store
    // result in temp_vec (vector resizes itself).
    temp_vec.clear();
    Bt_mat_vec_pt->multiply(another_temp_vec, temp_vec);


    if (Doc_time)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "LSC prec solve: G matrix vector product: "
                 << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }

    // NOTE: temp_vec now contains -G z_p

    // The vector another_temp_vec is no longer needed -- re-use it to store
    // displacement/position quantities:
    another_temp_vec.clear();

    // Loop over all enries in the global vector and find the
    // entries associated with the displacement/position:
    get_block_vector(0, r, another_temp_vec);
    another_temp_vec += temp_vec;

    // NOTE:  The vector another_temp_vec now contains r_u - G z_p

    // Solve momentum system
#ifdef PARANOID
    // check a solver has been set
    if (F_preconditioner_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "F_preconditioner_pt has not been set.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // use some Preconditioner's preconditioner solve
    // and return
    if (F_preconditioner_is_block_preconditioner)
    {
      return_block_vector(0, another_temp_vec, z);
      F_preconditioner_pt->preconditioner_solve(z, z);
    }
    else
    {
      F_preconditioner_pt->preconditioner_solve(another_temp_vec, temp_vec);
      return_block_vector(0, temp_vec, z);
    }

    if (Doc_time)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "LSC prec solve: F solve [nrow="
                 << P_preconditioner_pt->nrow() << "]: " << t_end - t_start
                 << std::endl;
      oomph_info << "LSC prec solve: Overall " << t_end - t_start_overall
                 << std::endl;
    }
  }


  //========================================================================
  /// Helper function to assemble the diagonal of the
  /// mass matrix from the elemental contributions defined in
  /// SolidElementWithDiagonalMassMatrix::get_mass_matrix_diagonal(...).
  //========================================================================
  CRDoubleMatrix* PressureBasedSolidLSCPreconditioner::
    assemble_mass_matrix_diagonal()
  {
    // determine the rows required by this processor
    unsigned first_row = this->block_distribution_pt(0)->first_row();
    unsigned nrow_local = this->block_distribution_pt(0)->nrow_local();
    unsigned nrow = this->block_distribution_pt(0)->nrow();

    // create storage for the diagonals
    double* m_values = new double[nrow_local];
    for (unsigned i = 0; i < nrow_local; i++)
    {
      m_values[i] = 0;
    }

    // if the problem is distributed
    bool distributed = false;
#ifdef OOMPH_HAS_MPI
    if (any_mesh_distributed() || this->master_distribution_pt()->distributed())
    {
      distributed = true;
    }
#endif

    // next we get the diagonal mass matrix data
    if (distributed)
    {
#ifdef OOMPH_HAS_MPI
      // the number of processors
      unsigned nproc = this->comm_pt()->nproc();

      // and my rank
      unsigned my_rank = this->comm_pt()->my_rank();

      // determine the rows for which we have lookup rows
      // if the problem is NOT distributed then we only classify global equation
      // on this processor to avoid duplication (as every processor holds
      // every element)
      unsigned first_lookup_row = 0;
      unsigned last_lookup_row = 0;
      first_lookup_row = this->master_distribution_pt()->first_row();
      last_lookup_row =
        first_lookup_row + this->master_distribution_pt()->nrow_local() - 1;

      // find number of local elements
      unsigned n_el = Solid_mesh_pt->nelement();

      // the diagonal mass matrix contributions that have been
      // classified and should be sent to another processor
      Vector<double>* classified_contributions_send = new Vector<double>[nproc];

      // the corresponding block indices
      Vector<unsigned>* classified_indices_send = new Vector<unsigned>[nproc];

      // the maitrix contributions that cannot be classified by this processor
      // and therefore must be sent to another for classification
      Vector<double>* unclassified_contributions_send =
        new Vector<double>[nproc];

      // the corresponding global indices that require classification
      Vector<unsigned>* unclassified_indices_send = new Vector<unsigned>[nproc];

      // get the master distribution pt
      const LinearAlgebraDistribution* master_distribution_pt =
        this->master_distribution_pt();

      // get the displacement/position distribution pt
      const LinearAlgebraDistribution* displ_dist_pt =
        this->block_distribution_pt(0);

      // get the contribution for each element
      for (unsigned e = 0; e < n_el; e++)
      {
        // check that the element is not halo d
        if (!Solid_mesh_pt->element_pt(e)->is_halo())
        {
          // find number of degrees of freedom in the element
          // (this is slightly too big because it includes the
          // pressure dofs but this doesn't matter)
          unsigned el_dof = Solid_mesh_pt->element_pt(e)->ndof();

          // allocate local storage for the element's contribution to the
          // mass matrix diagonal
          Vector<double> el_vmm_diagonal(el_dof);

          SolidElementWithDiagonalMassMatrix* cast_el_pt =
            dynamic_cast<SolidElementWithDiagonalMassMatrix*>(
              Solid_mesh_pt->element_pt(e));


          if (cast_el_pt == 0)
          {
#ifdef PARANOID
            std::ostringstream error_message;
            error_message << "Failed cast to "
                          << "SolidElementWithDiagonalMassMatrix*\n"
                          << "Element is of type: "
                          << typeid(*(Solid_mesh_pt->element_pt(e))).name()
                          << "\n"
                          << typeid(Solid_mesh_pt->element_pt(e)).name()
                          << std::endl;
            OomphLibWarning(error_message.str(),
                            "PressureBasedSolidLSCPreconditioner::assemble_"
                            "mass_matrix_diagonal()",
                            OOMPH_EXCEPTION_LOCATION);
#endif
          }
          else
          {
            cast_el_pt->get_mass_matrix_diagonal(el_vmm_diagonal);
          }

          // get the contribution for each dof
          for (unsigned i = 0; i < el_dof; i++)
          {
            // Get the equation number
            unsigned eqn_number = Solid_mesh_pt->element_pt(e)->eqn_number(i);

            // if I have lookup information on this processor
            if (eqn_number >= first_lookup_row && eqn_number <= last_lookup_row)
            {
              // bypass non displacement/position DOFs
              if (this->block_number(eqn_number) == 0)
              {
                // get the index in the block
                unsigned index = this->index_in_block(eqn_number);

                // determine which processor requires the block index
                for (unsigned p = 0; p < nproc; p++)
                {
                  if (index >= displ_dist_pt->first_row(p) &&
                      (index < (displ_dist_pt->first_row(p) +
                                displ_dist_pt->nrow_local(p))))
                  {
                    // if it is required by this processor then add the
                    // contribution
                    if (p == my_rank)
                    {
                      m_values[index - first_row] += el_vmm_diagonal[i];
                    }
                    // other wise store it for communication
                    else
                    {
                      classified_contributions_send[p].push_back(
                        el_vmm_diagonal[i]);
                      classified_indices_send[p].push_back(index);
                    }
                  }
                }
              }
            }

            // if we do not have the lookup information on this processor
            // then we send the mass matrix contribution to a processor
            // which we know does have the lookup information
            // the assumption: the processor for which the master block
            // preconditioner distribution will definitely hold the lookup
            // data for eqn_number (although others may)
            else if (any_mesh_distributed())
            {
              // determine which processor requires the block index
              unsigned p = 0;
              while (!(eqn_number >= master_distribution_pt->first_row(p) &&
                       (eqn_number < (master_distribution_pt->first_row(p) +
                                      master_distribution_pt->nrow_local(p)))))
              {
                p++;
              }

              // store the data
              unclassified_contributions_send[p].push_back(el_vmm_diagonal[i]);
              unclassified_indices_send[p].push_back(eqn_number);
            }
          }
        }
      }

      // next the unclassified contributions are communicated to
      // processors that can classify them

      // first determine how many unclassified rows are to be sent to
      // each processor
      unsigned* n_unclassified_send = new unsigned[nproc];
      for (unsigned p = 0; p < nproc; p++)
      {
        if (p == my_rank)
        {
          n_unclassified_send[p] = 0;
        }
        else
        {
          n_unclassified_send[p] = unclassified_contributions_send[p].size();
        }
      }

      // then all-to-all number of unclassified to be sent / recv
      unsigned* n_unclassified_recv = new unsigned[nproc];
      MPI_Alltoall(n_unclassified_send,
                   1,
                   MPI_UNSIGNED,
                   n_unclassified_recv,
                   1,
                   MPI_UNSIGNED,
                   this->comm_pt()->mpi_comm());

      // the base displacement for the sends
      MPI_Aint base_displacement;
      MPI_Get_address(m_values, &base_displacement);

      // allocate storage for the data to be recieved
      // and post the sends and recvs
      Vector<double*> unclassified_contributions_recv(nproc);
      Vector<unsigned*> unclassified_indices_recv(nproc);
      Vector<MPI_Request> unclassified_recv_requests;
      Vector<MPI_Request> unclassified_send_requests;
      Vector<unsigned> unclassified_recv_proc;
      for (unsigned p = 0; p < nproc; p++)
      {
        if (p != my_rank)
        {
          // recv
          if (n_unclassified_recv[p] > 0)
          {
            unclassified_contributions_recv[p] =
              new double[n_unclassified_recv[p]];
            unclassified_indices_recv[p] = new unsigned[n_unclassified_recv[p]];

            // data for the struct data type
            MPI_Datatype recv_types[2];
            MPI_Aint recv_displacements[2];
            int recv_sz[2];

            // contributions
            MPI_Type_contiguous(
              n_unclassified_recv[p], MPI_DOUBLE, &recv_types[0]);
            MPI_Type_commit(&recv_types[0]);
            MPI_Get_address(unclassified_contributions_recv[p],
                            &recv_displacements[0]);
            recv_displacements[0] -= base_displacement;
            recv_sz[0] = 1;

            // indices
            MPI_Type_contiguous(
              n_unclassified_recv[p], MPI_UNSIGNED, &recv_types[1]);
            MPI_Type_commit(&recv_types[1]);
            MPI_Get_address(unclassified_indices_recv[p],
                            &recv_displacements[1]);
            recv_displacements[1] -= base_displacement;
            recv_sz[1] = 1;

            // build the final recv type
            MPI_Datatype final_recv_type;
            MPI_Type_create_struct(
              2, recv_sz, recv_displacements, recv_types, &final_recv_type);
            MPI_Type_commit(&final_recv_type);

            // and recv
            MPI_Request req;
            MPI_Irecv(
              m_values, 1, final_recv_type, p, 0, comm_pt()->mpi_comm(), &req);
            unclassified_recv_requests.push_back(req);
            unclassified_recv_proc.push_back(p);
            MPI_Type_free(&recv_types[0]);
            MPI_Type_free(&recv_types[1]);
            MPI_Type_free(&final_recv_type);
          }

          // send
          if (n_unclassified_send[p] > 0)
          {
            // data for the struct data type
            MPI_Datatype send_types[2];
            MPI_Aint send_displacements[2];
            int send_sz[2];

            // contributions
            MPI_Type_contiguous(
              n_unclassified_send[p], MPI_DOUBLE, &send_types[0]);
            MPI_Type_commit(&send_types[0]);
            MPI_Get_address(&unclassified_contributions_send[p][0],
                            &send_displacements[0]);
            send_displacements[0] -= base_displacement;
            send_sz[0] = 1;

            // indices
            MPI_Type_contiguous(
              n_unclassified_send[p], MPI_UNSIGNED, &send_types[1]);
            MPI_Type_commit(&send_types[1]);
            MPI_Get_address(&unclassified_indices_send[p][0],
                            &send_displacements[1]);
            send_displacements[1] -= base_displacement;
            send_sz[1] = 1;

            // build the final send type
            MPI_Datatype final_send_type;
            MPI_Type_create_struct(
              2, send_sz, send_displacements, send_types, &final_send_type);
            MPI_Type_commit(&final_send_type);

            // and send
            MPI_Request req;
            MPI_Isend(
              m_values, 1, final_send_type, p, 0, comm_pt()->mpi_comm(), &req);
            unclassified_send_requests.push_back(req);
            MPI_Type_free(&send_types[0]);
            MPI_Type_free(&send_types[1]);
            MPI_Type_free(&final_send_type);
          }
        }
      }

      // next classify the data as it is received
      unsigned n_unclassified_recv_req = unclassified_recv_requests.size();
      while (n_unclassified_recv_req > 0)
      {
        // get the processor number and remove the completed request
        // for the vector of requests
        int req_num;
        MPI_Waitany(n_unclassified_recv_req,
                    &unclassified_recv_requests[0],
                    &req_num,
                    MPI_STATUS_IGNORE);
        unsigned p = unclassified_recv_proc[req_num];
        unclassified_recv_requests.erase(unclassified_recv_requests.begin() +
                                         req_num);
        unclassified_recv_proc.erase(unclassified_recv_proc.begin() + req_num);
        n_unclassified_recv_req--;

        // next classify the dofs
        // and store them for sending to other processors if required
        unsigned n_recv = n_unclassified_recv[p];
        for (unsigned i = 0; i < n_recv; i++)
        {
          unsigned eqn_number = unclassified_indices_recv[p][i];
          // bypass non displacement/position DOFs
          if (this->block_number(eqn_number) == 0)
          {
            // get the index in the block
            unsigned index = this->index_in_block(eqn_number);

            // determine which processor requires the block index
            for (unsigned pp = 0; pp < nproc; pp++)
            {
              if (index >= displ_dist_pt->first_row(pp) &&
                  (index < (displ_dist_pt->first_row(pp) +
                            displ_dist_pt->nrow_local(pp))))
              {
                // if it is required by this processor then add the
                // contribution
                if (pp == my_rank)
                {
                  m_values[index - first_row] +=
                    unclassified_contributions_recv[p][i];
                }
                // other wise store it for communication
                else
                {
                  double v = unclassified_contributions_recv[p][i];
                  classified_contributions_send[pp].push_back(v);
                  classified_indices_send[pp].push_back(index);
                }
              }
            }
          }
        }

        // clean up
        delete[] unclassified_contributions_recv[p];
        delete[] unclassified_indices_recv[p];
      }
      delete[] n_unclassified_recv;

      // now all indices have been classified

      // next the classified contributions are communicated to
      // processors that require them

      // first determine how many classified rows are to be sent to
      // each processor
      unsigned* n_classified_send = new unsigned[nproc];
      for (unsigned p = 0; p < nproc; p++)
      {
        if (p == my_rank)
        {
          n_classified_send[p] = 0;
        }
        else
        {
          n_classified_send[p] = classified_contributions_send[p].size();
        }
      }

      // then all-to-all com number of classified to be sent / recv
      unsigned* n_classified_recv = new unsigned[nproc];
      MPI_Alltoall(n_classified_send,
                   1,
                   MPI_UNSIGNED,
                   n_classified_recv,
                   1,
                   MPI_UNSIGNED,
                   this->comm_pt()->mpi_comm());

      // allocate storage for the data to be recieved
      // and post the sends and recvs
      Vector<double*> classified_contributions_recv(nproc);
      Vector<unsigned*> classified_indices_recv(nproc);
      Vector<MPI_Request> classified_recv_requests;
      Vector<MPI_Request> classified_send_requests;
      Vector<unsigned> classified_recv_proc;
      for (unsigned p = 0; p < nproc; p++)
      {
        if (p != my_rank)
        {
          // recv
          if (n_classified_recv[p] > 0)
          {
            classified_contributions_recv[p] = new double[n_classified_recv[p]];
            classified_indices_recv[p] = new unsigned[n_classified_recv[p]];

            // data for the struct data type
            MPI_Datatype recv_types[2];
            MPI_Aint recv_displacements[2];
            int recv_sz[2];

            // contributions
            MPI_Type_contiguous(
              n_classified_recv[p], MPI_DOUBLE, &recv_types[0]);
            MPI_Type_commit(&recv_types[0]);
            MPI_Get_address(classified_contributions_recv[p],
                            &recv_displacements[0]);
            recv_displacements[0] -= base_displacement;
            recv_sz[0] = 1;

            // indices
            MPI_Type_contiguous(
              n_classified_recv[p], MPI_UNSIGNED, &recv_types[1]);
            MPI_Type_commit(&recv_types[1]);
            MPI_Get_address(classified_indices_recv[p], &recv_displacements[1]);
            recv_displacements[1] -= base_displacement;
            recv_sz[1] = 1;

            // build the final recv type
            MPI_Datatype final_recv_type;
            MPI_Type_create_struct(
              2, recv_sz, recv_displacements, recv_types, &final_recv_type);
            MPI_Type_commit(&final_recv_type);

            // and recv
            MPI_Request req;
            MPI_Irecv(
              m_values, 1, final_recv_type, p, 0, comm_pt()->mpi_comm(), &req);
            classified_recv_requests.push_back(req);
            classified_recv_proc.push_back(p);
            MPI_Type_free(&recv_types[0]);
            MPI_Type_free(&recv_types[1]);
            MPI_Type_free(&final_recv_type);
          }

          // send
          if (n_classified_send[p] > 0)
          {
            // data for the struct data type
            MPI_Datatype send_types[2];
            MPI_Aint send_displacements[2];
            int send_sz[2];

            // contributions
            MPI_Type_contiguous(
              n_classified_send[p], MPI_DOUBLE, &send_types[0]);
            MPI_Type_commit(&send_types[0]);
            MPI_Get_address(&classified_contributions_send[p][0],
                            &send_displacements[0]);
            send_displacements[0] -= base_displacement;
            send_sz[0] = 1;

            // indices
            MPI_Type_contiguous(
              n_classified_send[p], MPI_UNSIGNED, &send_types[1]);
            MPI_Type_commit(&send_types[1]);
            MPI_Get_address(&classified_indices_send[p][0],
                            &send_displacements[1]);
            send_displacements[1] -= base_displacement;
            send_sz[1] = 1;

            // build the final send type
            MPI_Datatype final_send_type;
            MPI_Type_create_struct(
              2, send_sz, send_displacements, send_types, &final_send_type);
            MPI_Type_commit(&final_send_type);

            // and send
            MPI_Request req;
            MPI_Isend(
              m_values, 1, final_send_type, p, 0, comm_pt()->mpi_comm(), &req);
            classified_send_requests.push_back(req);
            MPI_Type_free(&send_types[0]);
            MPI_Type_free(&send_types[1]);
            MPI_Type_free(&final_send_type);
          }
        }
      }

      // next classify the data as it is received
      unsigned n_classified_recv_req = classified_recv_requests.size();
      while (n_classified_recv_req > 0)
      {
        // get the processor number and remove the completed request
        // for the vector of requests
        int req_num;
        MPI_Waitany(n_classified_recv_req,
                    &classified_recv_requests[0],
                    &req_num,
                    MPI_STATUS_IGNORE);
        unsigned p = classified_recv_proc[req_num];
        classified_recv_requests.erase(classified_recv_requests.begin() +
                                       req_num);
        classified_recv_proc.erase(classified_recv_proc.begin() + req_num);
        n_classified_recv_req--;

        // next classify the dofs
        // and store them for sending to other processors if required
        unsigned n_recv = n_classified_recv[p];
        for (unsigned i = 0; i < n_recv; i++)
        {
          m_values[classified_indices_recv[p][i] - first_row] +=
            classified_contributions_recv[p][i];
        }

        // clean up
        delete[] classified_contributions_recv[p];
        delete[] classified_indices_recv[p];
      }

      // wait for the unclassified sends to complete
      unsigned n_unclassified_send_req = unclassified_send_requests.size();
      if (n_unclassified_send_req > 0)
      {
        MPI_Waitall(n_unclassified_send_req,
                    &unclassified_send_requests[0],
                    MPI_STATUS_IGNORE);
      }
      delete[] unclassified_contributions_send;
      delete[] unclassified_indices_send;
      delete[] n_unclassified_send;

      // wait for the classified sends to complete
      unsigned n_classified_send_req = classified_send_requests.size();
      if (n_classified_send_req > 0)
      {
        MPI_Waitall(n_classified_send_req,
                    &classified_send_requests[0],
                    MPI_STATUS_IGNORE);
      }
      delete[] classified_indices_send;
      delete[] classified_contributions_send;
      delete[] n_classified_recv;
      delete[] n_classified_send;
#endif
    }

    // or if the problem is not distributed
    else
    {
      // find number of elements
      unsigned n_el = Solid_mesh_pt->nelement();

      // get the contribution for each element
      for (unsigned e = 0; e < n_el; e++)
      {
        // find number of degrees of freedom in the element
        // (this is slightly too big because it includes the
        // pressure dofs but this doesn't matter)
        unsigned el_dof = Solid_mesh_pt->element_pt(e)->ndof();

        // allocate local storage for the element's contribution to the
        //  mass matrix diagonal
        Vector<double> el_vmm_diagonal(el_dof);

        SolidElementWithDiagonalMassMatrix* cast_el_pt =
          dynamic_cast<SolidElementWithDiagonalMassMatrix*>(
            Solid_mesh_pt->element_pt(e));

        if (cast_el_pt == 0)
        {
#ifdef PARANOID
          // #pragma clang diagnostic push
          // #pragma clang diagnostic ignored
          // "-Wpotentially-evaluated-expression"
          std::ostringstream error_message;
          error_message << "Failed cast to "
                        << "SolidElementWithDiagonalMassMatrix*\n"
                        << "Element is of type: "
                        << typeid(*(Solid_mesh_pt->element_pt(e))).name()
                        << "\n"
                        << typeid(Solid_mesh_pt->element_pt(e)).name()
                        << std::endl;
          OomphLibWarning(error_message.str(),
                          "PressureBasedSolidLSCPreconditioner::assemble_mass_"
                          "matrix_diagonal()",
                          OOMPH_EXCEPTION_LOCATION);
//#pragma clang diagnostic pop
#endif
        }
        else
        {
          cast_el_pt->get_mass_matrix_diagonal(el_vmm_diagonal);
        }

        // get the contribution for each dof
        for (unsigned i = 0; i < el_dof; i++)
        {
          // Get the equation number
          unsigned eqn_number = Solid_mesh_pt->element_pt(e)->eqn_number(i);

          // bypass non displacement/position DOFs
          if (this->block_number(eqn_number) == 0)
          {
            // get the index in the block
            unsigned index = this->index_in_block(eqn_number);

            // if it is required on this processor
            if (index >= first_row && index < first_row + nrow_local)
            {
              m_values[index - first_row] += el_vmm_diagonal[i];
            }
          }
        }
      }
    }

    // create column index and row start
    int* m_column_index = new int[nrow_local];
    int* m_row_start = new int[nrow_local + 1];
    for (unsigned i = 0; i < nrow_local; i++)
    {
      m_values[i] = 1 / m_values[i];
      m_column_index[i] = first_row + i;
      m_row_start[i] = i;
    }
    m_row_start[nrow_local] = nrow_local;

    // build the matrix
    CRDoubleMatrix* m_pt = new CRDoubleMatrix(this->block_distribution_pt(0));
    m_pt->build_without_copy(
      nrow, nrow_local, m_values, m_column_index, m_row_start);

    // return the matrix;
    return m_pt;
  }

  //=========================================================================
  /// Helper function to delete preconditioner data.
  //=========================================================================
  void PressureBasedSolidLSCPreconditioner::clean_up_memory()
  {
    if (Preconditioner_has_been_setup)
    {
      // delete matvecs
      delete Bt_mat_vec_pt;
      Bt_mat_vec_pt = 0;

      if (!Form_BFBt_product)
      {
        delete F_mat_vec_pt;
        F_mat_vec_pt = 0;
        delete QBt_mat_vec_pt;
        QBt_mat_vec_pt = 0;
      }
      else
      {
        delete E_mat_vec_pt;
        E_mat_vec_pt = 0;
      }

      // delete stuff from displacement/position solve
      if (Using_default_f_preconditioner)
      {
        delete F_preconditioner_pt;
        F_preconditioner_pt = 0;
      }

      // delete stuff from Schur complement approx
      if (Using_default_p_preconditioner)
      {
        delete P_preconditioner_pt;
        P_preconditioner_pt = 0;
      }
    }
  }
} // namespace oomph
