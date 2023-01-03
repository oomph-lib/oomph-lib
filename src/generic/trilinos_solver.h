// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_TRILINOS_SOLVER_HEADER
#define OOMPH_TRILINOS_SOLVER_HEADER

#include "trilinos_preconditioners.h"
#include "iterative_linear_solver.h"

#include "AztecOO.h"

namespace oomph
{
  //=============================================================================
  /// An Epetra_Operator class for oomph-lib preconditioners.
  /// A helper class for TrilinosOomphLibPreconditioner to allow an oomph-lib
  /// preconditioner (i.e. one derived from Preconditioner) to be used with
  /// a trilinos solver (TrilinosAztecOOSolver)
  //=============================================================================
  class OomphLibPreconditionerEpetraOperator : public Epetra_Operator
  {
  public:
    /// Constructor - takes the pointer to the oomph-lib
    /// preconditioner and the distribution of the preconditioner
    /// Note: the oomph-lib preconditioner must be setup.
    /// If use_eptra_values is true then the epetra vector values is used
    /// within the vectors passed to the oomph-lib
    /// preconditioner. If this is true none of the vector rebuild methods can
    /// be called.
    OomphLibPreconditionerEpetraOperator(Preconditioner* preconditioner_pt,
                                         bool use_epetra_values = false)
#ifdef OOMPH_HAS_MPI
      : Operator_comm(
          preconditioner_pt->distribution_pt()->communicator_pt()->mpi_comm()),
        Use_epetra_values(use_epetra_values)
#else
      : Operator_comm(), Use_epetra_values(use_epetra_values)
#endif
    {
      // set the ooomph-lib preconditioner
      Oomph_lib_preconditioner_pt = preconditioner_pt;

      // set the preconditioner label
      Preconditioner_label = "oomph-lib Preconditioner";

      // setup the Epetra_map
      Operator_map_pt = TrilinosEpetraHelpers::create_epetra_map(
        Oomph_lib_preconditioner_pt->distribution_pt());
    }

    /// Destructor - deletes the Epetra_map and My_global_rows vector
    /// (if MPI)
    ~OomphLibPreconditionerEpetraOperator()
    {
      delete Operator_map_pt;
      Operator_map_pt = 0;
    }

    /// Broken copy constructor
    OomphLibPreconditionerEpetraOperator(
      const OomphLibPreconditionerEpetraOperator&) = delete;

    /// Broken assignment operator.
    void operator=(const OomphLibPreconditionerEpetraOperator&) = delete;

    /// Broken Epetra_Operator member - SetUseTranspose
    int SetUseTranspose(bool UseTranspose)
    {
      std::ostringstream error_message;
      error_message << "SetUseTranspose() is a pure virtual Epetra_Operator "
                    << "member that is not required for a Preconditioner"
                    << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    /// Broken Epetra_Operator member - Apply
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      std::ostringstream error_message;
      error_message << "Apply() is a pure virtual Epetra_Operator member"
                    << "that is not required for a Preconditioner" << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


    /// applies the oomph-lib preconditioner. Converts the Epetra vector
    /// applys the preconditioner by calling the oomph-lib preconditioner's
    /// preconditioner_solve functionality.
    /// NOTE : the oomph-lib preconditioner is setup prior to being passed to
    /// this class
    int ApplyInverse(const Epetra_MultiVector& epetra_r,
                     Epetra_MultiVector& epetra_z) const
    {
      // the oomph-lib vector for r
      DoubleVector oomph_r;

      // copy the Epetra_MultiVector r into an oomph-lib vector
      double** r_pt;
      epetra_r.ExtractView(&r_pt);
      if (Use_epetra_values)
      {
        oomph_r.set_external_values(
          Oomph_lib_preconditioner_pt->distribution_pt(), *r_pt, false);
      }
      else
      {
        oomph_r.build(Oomph_lib_preconditioner_pt->distribution_pt(), 0.0);
        unsigned nrow_local =
          Oomph_lib_preconditioner_pt->distribution_pt()->nrow_local();
        for (unsigned i = 0; i < nrow_local; i++)
        {
          oomph_r[i] = r_pt[0][i];
        }
      }

      // oomph-lib vector for Y
      DoubleVector oomph_z;
      if (Use_epetra_values)
      {
        double** z_pt;
        epetra_z.ExtractView(&z_pt);
        DoubleVector oomph_z;
        oomph_z.set_external_values(
          Oomph_lib_preconditioner_pt->distribution_pt(), *z_pt, false);
      }
      else
      {
        oomph_z.build(Oomph_lib_preconditioner_pt->distribution_pt(), 0.0);
      }

      // apply the preconditioner
      Oomph_lib_preconditioner_pt->preconditioner_solve(oomph_r, oomph_z);

      // if not using external data copy the oomph-lib vector oomph_Y back to Y
      if (!Use_epetra_values)
      {
        unsigned nrow_local =
          Oomph_lib_preconditioner_pt->distribution_pt()->nrow_local();
        for (unsigned i = 0; i < nrow_local; i++)
        {
          epetra_z.ReplaceMyValue(i, 0, oomph_z[i]);
        }
      }

      // return 0 to indicate success
      return 0;
    }

    /// Broken Epetra_Operator member - NormInf
    double NormInf() const
    {
      std::ostringstream error_message;
      error_message << "NormInf() is a pure virtual Epetra_Operator member"
                    << "that is not required for a Preconditioner" << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    /// Epetra_Operator::Label - returns a string describing the operator
    const char* Label() const
    {
      return Preconditioner_label.c_str();
    }

    /// Broken Epetra_Operator member - UseTranspose
    bool UseTranspose() const
    {
      std::ostringstream error_message;
      error_message
        << "UseTranspose() is a pure virtual Epetra_Operator member "
        << "that is not required for a Preconditioner" << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    /// Broken Epetra_Operator member - HasNormInf
    bool HasNormInf() const
    {
      std::ostringstream error_message;
      error_message << "HasNormInf() is a pure virtual Epetra_Operator member "
                    << "that is not required for a Preconditioner" << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    /// Returns the Epetra MPI_Comm object
    const Epetra_Comm& Comm() const
    {
      return Operator_comm;
    }

    /// Epetra_Operator member - OperatorDomainMap
    const Epetra_Map& OperatorDomainMap() const
    {
      return *Operator_map_pt;
    }

    /// Epetra_Operator member - OperatorRangeMap
    const Epetra_Map& OperatorRangeMap() const
    {
      return *Operator_map_pt;
    }


  private:
    /// A pointer to the oomph-lib preconditioner
    Preconditioner* Oomph_lib_preconditioner_pt;

#ifdef OOMPH_HAS_MPI
    /// An Epetra MPI Comm object
    Epetra_MpiComm Operator_comm;
#else
    /// An Epetra Serial Comm object
    Epetra_SerialComm Operator_comm;
#endif

    /// Use the epetra data within the vectors passed to the oomph-lib
    /// preconditioner. If this is true none of the vector rebuild methods can
    /// be called.
    bool Use_epetra_values;

    /// A pointer to an Epetra_Map object - describes distribution of the
    /// preconditioner, in this instance it is primarily used to prescribe the
    /// distribution
    /// of the residual and solution vector
    Epetra_Map* Operator_map_pt;

    /// a label for the preconditioner ( for Epetra_Operator::Label() )
    std::string Preconditioner_label;
  };


  //=============================================================================
  /// An interface to the Trilinos AztecOO classes allowing it to be used
  /// as an Oomph-lib LinearSolver.
  /// The AztecOO solver is a Krylov Subspace solver; the solver type (either
  /// CG, GMRES or BiCGStab) can be set using solver_type(). This solver can be
  /// preconditioned with Trilinos Preconditioners (derived from
  /// TrilinosPreconditionerBase) or Oomph-lib preconditioners (derived from
  /// Preconditioner). Preconditioners are set using preconditioner_pt().
  //=============================================================================
  class TrilinosAztecOOSolver : public IterativeLinearSolver
  {
  public:
    /// Constructor.
    TrilinosAztecOOSolver()
    {
      // By default use workaround for creating of epetra matrix that
      // respects aztecoo's ordering requirements
      Use_aztecoo_workaround_for_epetra_matrix_setup = true;

      // set pointer to Null
      AztecOO_solver_pt = 0;

      // initially assume not problem based solve
      Using_problem_based_solve = false;

      // if a problem based solve is performed then it should generate a
      // serial matrix
      Assemble_serial_jacobian = false;

      // by default we do not use external values in the oomph-lib
      // preconditioner
      If_oomphlib_preconditioner_use_epetra_values = false;

      // null the pts
      Problem_pt = 0;
      Epetra_matrix_pt = 0;
      Epetra_preconditioner_pt = 0;

      // set solver defaults
      Solver_type = GMRES;
      Tolerance = 10e-10;

      // don't delete the matrix
      Delete_matrix = false;
    }

    /// Destructor - delete the solver and the matrices
    ~TrilinosAztecOOSolver()
    {
      // delete
      clean_up_memory();

      // if Problem_based solve then the oomph matrix was generated by this
      // class and should be deleted
      if (Using_problem_based_solve)
      {
        delete Oomph_matrix_pt;
        Oomph_matrix_pt = 0;
      }
    }

    /// Broken copy constructor.
    TrilinosAztecOOSolver(const TrilinosAztecOOSolver&) = delete;

    /// Broken assignment operator.
    void operator=(const TrilinosAztecOOSolver&) = delete;

    /// Enable workaround for creating of epetra matrix that respects
    /// aztecoo's ordering requirements
    void enable_aztecoo_workaround_for_epetra_matrix_setup()
    {
      Use_aztecoo_workaround_for_epetra_matrix_setup = true;
    }

    /// Disable workaround for creating of epetra matrix that respects
    /// aztecoo's ordering requirements
    void disable_aztecoo_workaround_for_epetra_matrix_setup()
    {
      Use_aztecoo_workaround_for_epetra_matrix_setup = false;
    }

    /// Is workaround for creating of epetra matrix that respects
    /// aztecoo's ordering requirements enabled?
    bool is_aztecoo_workaround_for_epetra_matrix_setup_enabled()
    {
      return Use_aztecoo_workaround_for_epetra_matrix_setup;
    }

    /// Clean up method - deletes the solver, the matrices and the
    /// preconditioner
    void clean_up_memory()
    {
      // delete the solver
      delete AztecOO_solver_pt;
      AztecOO_solver_pt = 0;

      // delete the Epetra_Operator preconditioner (only if it is a wrapper to
      // an oomphlib preconditioner in which case only the wrapper is deleted
      // and not the actual preconditioner). if the preconditioner is Trilinos
      // preconditioner then the Epetra_Operator is deleted when that
      // preconditioner is deleted
      if (Epetra_preconditioner_pt != 0)
      {
        if (dynamic_cast<OomphLibPreconditionerEpetraOperator*>(
              Epetra_preconditioner_pt) != 0)
        {
          delete Epetra_preconditioner_pt;
          Epetra_preconditioner_pt = 0;
        }
      }

      // don't delete the preconditioner but call its clean up memory method
      if (this->preconditioner_pt() != 0)
      {
        this->preconditioner_pt()->clean_up_memory();
      }


      // delete the matrices
      // This must now happen after the preconditioner delete because the
      // ML preconditioner (and maybe others) use the communicator from the
      // matrix in the destructor
      delete Epetra_matrix_pt;
      Epetra_matrix_pt = 0;
    }

    /// Function which uses problem_pt's get_jacobian(...) function to
    /// generate a linear system which is then solved. This function deletes
    /// any existing internal data and then generates a new AztecOO solver.
    void solve(Problem* const& problem_pt, DoubleVector& solution);

    /// Function to solve the linear system defined by matrix_pt and rhs.
    /// \b NOTE 1. The matrix has to be of type CRDoubleMatrix or
    /// DistributedCRDoubleMatrix.
    /// \b NOTE 2. This function will delete any existing internal data and
    /// generate a new AztecOO solver.
    void solve(DoubleMatrixBase* const& matrix_pt,
               const DoubleVector& rhs,
               DoubleVector& solution);

    /// Function to resolve a linear system using the existing solver
    /// data, allowing a solve with a new right hand side vector. This
    /// function must be used after a call to solve(...) with
    /// enable_resolve set to true.
    void resolve(const DoubleVector& rhs, DoubleVector& solution);

    /// Disable resolve function (overloads the LinearSolver
    /// disable_resolve function).
    void disable_resolve()
    {
      Enable_resolve = false;
      clean_up_memory();
    }

    /// Call if the matrix can be deleted
    void enable_delete_matrix()
    {
      Delete_matrix = true;
    }

    /// Call if the matrix can not be deleted (default)
    void disable_delete_matrix()
    {
      Delete_matrix = false;
    }

    /// Access function to Max_iter
    unsigned& max_iter()
    {
      return Max_iter;
    }

    /// Acess function to Iterations
    unsigned iterations() const
    {
      return Iterations;
    }

    /// Access function to Tolerance
    double& tolerance()
    {
      return Tolerance;
    }

    /// Access function to Solver_type
    unsigned& solver_type()
    {
      return Solver_type;
    }

    /// Function to return Jacobian_setup_time;
    double jacobian_setup_time()
    {
      return Jacobian_setup_time;
    }

    /// Function to return Linear_solver_solution_time
    double linear_solver_solution_time()
    {
      return Linear_solver_solution_time;
    }

    /// Set the assembly of the serial jacobian
    /// when performing a problem-based solve
    void enable_assemble_serial_jacobian()
    {
      Assemble_serial_jacobian = true;
    }

    /// Unset the assembly of the serial jacobian
    void disable_assemble_serial_jacobian()
    {
      Assemble_serial_jacobian = false;
    }

    /// if this solver is using an oomph-lib preconditioner then the vectors
    /// passed to preconditioner_solve(...) should be using the values in the
    /// epetra vector as external data. If the vectors are using external
    /// data then rebuild(...) methods cannot be used be used in the
    /// preconditioner.
    // bool& if_oomphlib_preconditioner_use_epetra_values()
    // {
    //  return If_oomphlib_preconditioner_use_epetra_values;
    // }

    /// Enumerated list to define which AztecOO solver is used
    enum AztecOO_solver_types
    {
      CG,
      GMRES,
      BiCGStab
    };

  protected:
    /// Helper function performs the actual solve once the AztecOO
    /// solver is set up
    void solve_using_AztecOO(Epetra_Vector*& rhs_pt, Epetra_Vector*& soln_pt);

    /// Helper function for setting up the solver. Converts the oomph-lib
    /// matrices to Epetra matrices, sets up the preconditioner, creates the
    /// Trilinos Aztec00 solver and passes in the matrix, preconditioner and
    /// parameters.
    void solver_setup(DoubleMatrixBase* const& matrix_pt);

    /// Use workaround for creating of epetra matrix that respects
    /// aztecoo's ordering requirements
    bool Use_aztecoo_workaround_for_epetra_matrix_setup;

    /// Stores number of iterations used
    unsigned Iterations;

    /// Pointer to the AztecOO solver
    AztecOO* AztecOO_solver_pt;

    /// Stores set up time for Jacobian
    double Jacobian_setup_time;

    /// Stores time for the solution (excludes time to set up preconditioner)
    double Linear_solver_solution_time;

    /// Trilinos copies matrix data from oomph-lib's own CRDoubleMatrix
    /// or DistributedCRDoubleMatrix to Trilinos's Epetra format - the Trilinos
    /// solver no longer requires the oomph-lib matrices and therefore they
    /// could be deleted to save memory. This must be requested explicitly by
    /// setting this flag to TRUE. \b NOTE: The matrix is deleted after the
    /// preconditioner is setup.
    bool Delete_matrix;

    /// If true, when performing a problem based solve a serial matrix
    /// will be requested from Problem::get_jacobian(...). Defaults to true
    bool Assemble_serial_jacobian;

    /// Defines which solver is set up - available types are
    /// defined in AztecOO_solver_types
    unsigned Solver_type;

    /// Helper flag keeping track of whether we called the
    /// linear algebra or problem-based solve function.
    bool Using_problem_based_solve;

    /// A pointer for the linear system matrix in Epetra_CrsMatrix format
    Epetra_CrsMatrix* Epetra_matrix_pt;

    /// A pointer to the Epetra_Operator for the preconditioner. This is
    /// only used if the preconditioner NOT a Trilinos preconditioner.
    Epetra_Operator* Epetra_preconditioner_pt;

    /// Oomph lib matrix pointer
    DoubleMatrixBase* Oomph_matrix_pt;

    /// A pointer to the underlying problem (NULL if MATRIX based solve)
    /// The problem_pt is stored here in a problem based solve for the
    /// preconditioner
    Problem* Problem_pt;

    /// if this solver is using an oomph-lib preconditioner then the vectors
    /// passed to preconditioner_solve(...) should be using the values in the
    /// epetra vector as external data. If the vectors are using external
    /// data then rebuild(...) methods cannot be used.
    bool If_oomphlib_preconditioner_use_epetra_values;
  };

} // namespace oomph
#endif
