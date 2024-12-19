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
#ifndef OOMPH_TRILINOS_OPERATORS_HEADER
#define OOMPH_TRILINOS_OPERATORS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

// Needed in trilinos (as of gcc 4.6.* or so...)
#include <cstddef>

// trilinos headers
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Ifpack.h"

// oomph-lib headers
#include "trilinos_helpers.h"
#include "preconditioner.h"

namespace oomph
{
  //=============================================================================
  /// Base class for Trilinos preconditioners as oomph-lib preconditioner.
  //=============================================================================
  class TrilinosPreconditionerBase : public Preconditioner
  {
  public:
    /// Constructor.
    TrilinosPreconditionerBase()
    {
      // Initialise pointers
      Epetra_preconditioner_pt = 0;
      Epetra_matrix_pt = 0;
    }

    /// Destructor.
    virtual ~TrilinosPreconditionerBase()
    {
      clean_up_memory();
    }

    /// Static double that accumulates the preconditioner
    /// solve time of all instantiations of this class. Reset
    /// it manually, e.g. after every Newton solve.
    static double Cumulative_preconditioner_solve_time;

    /// deletes the preconditioner, matrices and maps
    void clean_up_memory()
    {
      // delete the Epetra preconditioner
      delete Epetra_preconditioner_pt;
      Epetra_preconditioner_pt = 0;

      // delete the epetra matrix
      delete Epetra_matrix_pt;
      Epetra_matrix_pt = 0;
    }

    /// Broken copy constructor.
    TrilinosPreconditionerBase(const TrilinosPreconditionerBase&) = delete;

    /// Broken assignment operator.
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const TrilinosPreconditionerBase&) = delete;*/

    /// Function to set up a preconditioner for the linear system
    /// defined by matrix_pt. This function must be called before using
    /// preconditioner_solve.
    /// \b NOTE 1. matrix_pt must point to an object of class CRDoubleMatrix or
    /// DistributedCRDoubleMatrix
    /// This method should be called by oomph-lib solvers and preconditioners
    void setup();

    /// Function to setup a preconditioner for the linear system defined
    /// by the oomph-lib oomph_matrix_pt and Epetra epetra_matrix_pt matrices.
    /// This method is called by Trilinos solvers.
    void setup(Epetra_CrsMatrix* epetra_matrix_pt);

    /// applies the preconditioner
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    /// Access function to Epetra_preconditioner_pt.
    /// For use with \c TrilinosAztecOOSolver
    Epetra_Operator*& epetra_operator_pt()
    {
      return Epetra_preconditioner_pt;
    }

    /// Access function to Epetra_preconditioner_pt (const version)
    /// For use with \c TrilinosAztecOOSolver
    Epetra_Operator* epetra_operator_pt() const
    {
      return Epetra_preconditioner_pt;
    }

  protected:
    /// Function to set up a specific Trilinos preconditioner.
    /// This is called by setup(...).
    virtual void setup_trilinos_preconditioner(
      Epetra_CrsMatrix* epetra_matrix_pt) = 0;

    /// The preconditioner which will be set up using function
    /// setup_trilinos_preconditioner(...)
    Epetra_Operator* Epetra_preconditioner_pt;

    /// Pointer used to store the epetra matrix - only used when this
    /// preconditioner is setup using the oomph-lib interface
    Epetra_CrsMatrix* Epetra_matrix_pt;
  };


  //============================================================================
  /// An interface to the Trilinos ML class - provides a function
  /// to construct a serial ML object, and functions to modify some
  /// of the ML paramaters.
  //============================================================================
  class TrilinosMLPreconditioner : public TrilinosPreconditionerBase
  {
  public:
    /// Constructor. Build with Smooth Aggretation (SA) default
    /// settings, but our own default number of V cycles (initialised
    /// to 1 to replicate TrilinosML's own behaviour).
    TrilinosMLPreconditioner()
    {
      // set default values
      ML_Epetra::SetDefaults("SA", ML_parameters);

      // Set number of MG cycles performed in preconditioner
      ML_parameters.set("cycle applications", Default_n_cycles);
    }

    /// Destructor empty -- clean up is done in base class
    virtual ~TrilinosMLPreconditioner() {}

    /// Broken copy constructor.
    TrilinosMLPreconditioner(const TrilinosMLPreconditioner&) = delete;

    /// Broken assignment operator.
    /*void operator=(const TrilinosMLPreconditioner&) = delete;*/

    /// Set control flags to values for Petrov-Galerkin
    /// preconditioning - for non symmetric systems
    void set_NSSA_default_values()
    {
      ML_Epetra::SetDefaults("NSSA", ML_parameters);
    }


    /// Set control flags to values for classical smoothed aggregation-
    /// based 2-level domain decomposition
    void set_DD_default_values()
    {
      ML_Epetra::SetDefaults("DD", ML_parameters);
    }


    /// Set control flags to values 3-level algebraic domain
    /// decomposition
    void set_DDML_default_values()
    {
      ML_Epetra::SetDefaults("DD-ML", ML_parameters);
    }

    /// Set control flags to values for classical smoothed
    /// aggregation preconditioning
    void set_SA_default_values()
    {
      ML_Epetra::SetDefaults("SA", ML_parameters);
    }

    /// Function to set maximum number of levels
    void set_max_levels(int max_levels)
    {
      ML_parameters.set("max levels", max_levels);
    }

    /// Function to set the number of cycles used
    void set_n_cycles(int n_cycles)
    {
      ML_parameters.set("cycle applications", n_cycles);
    }

    /// Function to set Smoother_damping
    void set_smoother_damping(double smoother_damping)
    {
      ML_parameters.set("smoother: damping factor", smoother_damping);
    }

    /// Function to set Smoother_sweeps
    void set_smoother_sweeps(int smoother_sweeps)
    {
      ML_parameters.set("smoother: sweeps", smoother_sweeps);
    }

    /// Function to set smoother type to "Jacobi"
    void set_smoother_jacobi()
    {
      ML_parameters.set("smoother: type", "Jacobi");
    }

    /// Function to set smoother type to "symmetric Gauss-Seidel"
    void set_smoother_gauss_seidel()
    {
      ML_parameters.set("smoother: type", "symmetric Gauss-Seidel");
    }

    /// Function to set output - controls level of information output by ML
    void set_output(int output)
    {
      ML_parameters.set("output", output);
    }


    /// Default number of V cycles (one to be consistent with
    /// previous default) (It's an int because Trilinos wants it to be!)
    static int Default_n_cycles;

  protected:
    /// Function to set up the ML preconditioner. It is assumed
    /// Trilinos_matrix_pt points to a suitable matrix.
    void setup_trilinos_preconditioner(Epetra_CrsMatrix* epetra_matrix_pt);

    // Parameter list of control flags for the preconditioner
    Teuchos::ParameterList ML_parameters;
  };


  //============================================================================
  /// An interface to the Trilinos IFPACK class- provides a function
  /// to construct an IFPACK object, and functions to modify some
  /// of the IFPACK paramaters.
  //============================================================================
  class TrilinosIFPACKPreconditioner : public TrilinosPreconditionerBase
  {
  public:
    /// Constructor.
    TrilinosIFPACKPreconditioner()
    {
      // set default values
      Preconditioner_type = "ILU";
      ILU_fill_level = 0;
      ILUT_fill_level = 1.0;
      Overlap = 0;
      Absolute_threshold = 0.0;
      Relative_threshold = 1.0;
    }

    /// Destructor -- empty, cleanup is done in base class
    virtual ~TrilinosIFPACKPreconditioner() {}

    /// Broken copy constructor.
    TrilinosIFPACKPreconditioner(const TrilinosIFPACKPreconditioner&) = delete;

    /// Broken assignment operator.
    /*void operator=(const TrilinosIFPACKPreconditioner&) = delete;*/

    /// Function to set Preconditioner_type to "ILU"
    void set_preconditioner_ILU()
    {
      Preconditioner_type = "ILU";
    }

    /// Function to set Preconditioner_type to "ILUT"
    void set_preconditioner_ILUT()
    {
      Preconditioner_type = "ILUT";
    }

    /// Access function for ILU_fill_level
    int& ilu_fill_level()
    {
      return ILU_fill_level;
    }

    /// Access function for ILUT_fill_level
    double& ilut_fill_level()
    {
      return ILUT_fill_level;
    }

    /// Access function for the absolute threshold
    double& absolute_threshold()
    {
      return Absolute_threshold;
    }

    /// Access function for the relative threshold
    double& relative_threshold()
    {
      return Relative_threshold;
    }

  protected:
    /// Function to set up an IFPACK preconditioner. It is assumed
    /// Trilinos_matrix_pt points to a suitable matrix.
    void setup_trilinos_preconditioner(Epetra_CrsMatrix* epetra_matrix_pt);

    /// Type of ILU preconditioner
    string Preconditioner_type;

    /// Level of fill for "ILU"
    int ILU_fill_level;

    /// Level of fill for "ILUT"
    double ILUT_fill_level;

    /// Value of overlap level - used in parallel ILU
    int Overlap;

    /// Value of absolute threshold, used to peturb diagonal
    double Absolute_threshold;

    /// Value of relative threshold, used to pertub diagonal
    double Relative_threshold;
  };
} // namespace oomph
#endif
