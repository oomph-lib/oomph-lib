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
#include "trilinos_preconditioners.h"

namespace oomph
{
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////
  // functions for the TrilinosPreconditionerBase class
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////


  //=============================================================================
  /// Static double that accumulates the preconditioner
  /// solve time of all instantiations of this class. Reset
  /// it manually, e.g. after every Newton solve.
  //=============================================================================
  double TrilinosPreconditionerBase::Cumulative_preconditioner_solve_time = 0.0;


  //=============================================================================
  /// Function to set up a preconditioner for the linear system
  /// defined by matrix_pt. This function must be called before using
  /// preconditioner_solve.
  /// \b NOTE 1. matrix_pt must point to an object of class CRDoubleMatrix or
  /// DistributedCRDoubleMatrix
  /// This method should be called by oomph-lib solvers and preconditioners
  //=============================================================================
  void TrilinosPreconditionerBase::setup()
  {
    // clean up the memory
    clean_up_memory();

#ifdef PARANOID
    // check the matrix is square
    if (matrix_pt()->nrow() != matrix_pt()->ncol())
    {
      std::ostringstream error_message;
      error_message << "Preconditioners require a square matrix. "
                    << "Matrix is " << matrix_pt()->nrow() << " by "
                    << matrix_pt()->ncol() << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // get a pointer to the cr double matrix
    CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt());

#ifdef PARANOID
    if (cr_matrix_pt == 0)
    {
      throw OomphLibError("Matrix must be of type CRDoubleMatrix",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // build the distribution of this preconditioner
    this->build_distribution(cr_matrix_pt->distribution_pt());

    // create the matrices
    Epetra_matrix_pt = TrilinosEpetraHelpers::create_distributed_epetra_matrix(
      cr_matrix_pt, this->distribution_pt());

    // set up preconditioner
    setup_trilinos_preconditioner(Epetra_matrix_pt);
  }

  //===========================================================================
  /// Function to setup a preconditioner for the linear system defined
  /// by the oomph-lib oomph_matrix_pt and Epetra epetra_matrix_pt matrices.
  /// This method is called by Trilinos solvers.
  //===========================================================================
  void TrilinosPreconditionerBase::setup(Epetra_CrsMatrix* epetra_matrix_pt)
  {
    // clean up old data
    clean_up_memory();

    // first try CRDoubleMatrix
    CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt());
    if (cr_matrix_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "TrilinosSolver only work with "
                    << "DistributedCRDoubleMatrix matrices" << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // setup the specific preconditioner
    setup_trilinos_preconditioner(epetra_matrix_pt);
  }


  //=============================================================================
  /// preconditioner_solve - applies the preconditioner to the vector r
  /// taking distributed oomph-lib vectors (DistributedVector<double>) as
  /// arguments.
  //=============================================================================
  void TrilinosPreconditionerBase::preconditioner_solve(const DoubleVector& r,
                                                        DoubleVector& z)
  {
    // Get ready to do cumulative timings
    double t_start = TimingHelpers::timer();

#ifdef PARANOID
    // check preconditioner data exists
    if (Epetra_preconditioner_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "preconditioner_solve requires that solver data has "
                    << "been set up" << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    if (*this->distribution_pt() != r.distribution_pt())
    {
      std::ostringstream error_message;
      error_message << "The rhs vector and the matrix must have the same "
                    << "distribution.\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // create Epetra version of r
    Epetra_Vector* epetra_r_pt =
      TrilinosEpetraHelpers::create_distributed_epetra_vector(r);

    // create an empty Epetra vector for z
    Epetra_Vector* epetra_z_pt =
      TrilinosEpetraHelpers::create_distributed_epetra_vector(
        this->distribution_pt());

    // Apply preconditioner
    Epetra_preconditioner_pt->ApplyInverse(*epetra_r_pt, *epetra_z_pt);

    // Copy result to z
    z.build(this->distribution_pt(), 0.0);
    TrilinosEpetraHelpers::copy_to_oomphlib_vector(epetra_z_pt, z);

    // clean up memory
    delete epetra_r_pt;
    delete epetra_z_pt;


    // Add to cumulative solve time
    double t_end = TimingHelpers::timer();
    Cumulative_preconditioner_solve_time += (t_end - t_start);
  }


  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////
  // function for the TrilinosML class
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////


  //=============================================================================
  /// (Static) default number of V cycles (one to be consistent
  /// with previous default). (It's an int because Trilinos wants it to be!)
  //=============================================================================
  int TrilinosMLPreconditioner::Default_n_cycles = 1;

  //=============================================================================
  // Function to set up the ML preconditioner.
  //=============================================================================
  void TrilinosMLPreconditioner::setup_trilinos_preconditioner(
    Epetra_CrsMatrix* epetra_matrix_pt)
  {
    // doc setup time
    oomph_info << "Setting up TrilinosML, ";
    double t_start = TimingHelpers::timer();


    // create the preconditioner
    Epetra_preconditioner_pt = new ML_Epetra::MultiLevelPreconditioner(
      *epetra_matrix_pt, ML_parameters, true);

    // doc setup time
    oomph_info << "time for setup [s] : " << TimingHelpers::timer() - t_start
               << std::endl;

    // oomph_info << "In here\n";
    // ML_Epetra::MultiLevelPreconditioner* tmp_pt=0;
    // tmp_pt=dynamic_cast<ML_Epetra::MultiLevelPreconditioner*>(
    //  Epetra_preconditioner_pt);
    // if (tmp_pt!=0)
    //  {
    //   oomph_info << "Doing test\n";
    //   ML_parameters.print(*(oomph_info.stream_pt()));
    //   oomph_info.stream_pt()->flush();
    //   // tmp_pt->TestSmoothers();
    //   oomph_info << "Done test\n";
    //  }
  }


  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////
  // function for the TrilinosIFPACK
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////


  //=============================================================================
  // Function to set up the IFPACK preconditioner.
  //=============================================================================
  void TrilinosIFPACKPreconditioner::setup_trilinos_preconditioner(
    Epetra_CrsMatrix* epetra_matrix_pt)
  {
    // set up a parameter list
    Teuchos::ParameterList ifpack_parameters;
    ifpack_parameters.set("fact: level-of-fill", ILU_fill_level);
    ifpack_parameters.set("fact: ilut level-of-fill", ILUT_fill_level);
    ifpack_parameters.set("fact: absolute threshold", Absolute_threshold);
    ifpack_parameters.set("fact: relative threshold", Relative_threshold);

    // create the preconditioner
    Ifpack ifpack_factory;
    Ifpack_Preconditioner* tmp_pt =
      ifpack_factory.Create(Preconditioner_type, epetra_matrix_pt, Overlap);

    tmp_pt->SetParameters(ifpack_parameters);
    tmp_pt->Initialize();
    tmp_pt->Compute();

    // Now store the pointer to the newly created Ifpack_Preconditioner
    Epetra_preconditioner_pt = tmp_pt;
  }


  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////


} // namespace oomph
