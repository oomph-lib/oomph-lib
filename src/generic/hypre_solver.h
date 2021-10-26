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
#ifndef OOMPH_HYPRE_SOLVER_HEADER
#define OOMPH_HYPRE_SOLVER_HEADER

// Headers required for hypre
#include "_hypre_utilities.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"

// OOMPH-LIB includes
#include "linear_algebra_distribution.h"
#include "linear_solver.h"
#include "preconditioner.h"

// hypre's global error flag
extern int hypre__global_error;


namespace oomph
{
  //==================================================================
  /// Helper functions for use with the Hypre library
  //==================================================================
  namespace HypreHelpers
  {
    /// Default for AMG strength (0.25 recommended for 2D problems;
    /// larger (0.5-0.75, say) for 3D
    extern double AMG_strength;

    /// Default AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    ///  3 = modified RS with 3rd pass to add C points on the boundaries
    ///  6 = Falgout (uses 1 then CLJP using interior coarse points as
    ///      first independent set)
    ///  8 = PMIS (parallel coarsening using independent sets - lower
    ///      complexities than 0, maybe also slower convergence)
    ///  10= HMIS (one pass RS on each processor then PMIS on interior
    ///      coarse points as first independent set)
    ///  11= One pass RS on each processor (not recommended)
    extern unsigned AMG_coarsening;

    /// AMG interpolation truncation factor
    extern double AMG_truncation; // 0.0;

    /// Helper function to check the Hypre error flag, return the
    /// message associated with any error, and reset the error flag to zero.
    int check_HYPRE_error_flag(std::ostringstream& message);

    /// Helper function to create a HYPRE_IJVector and HYPRE_ParVector.
    /// + If no MPI then serial vectors are created
    /// + If MPI and serial input vector then distributed hypre vectors are
    ///   created
    /// + If MPI and distributed input vector the distributed output vectors
    ///   are created.
    extern void create_HYPRE_Vector(const DoubleVector& oomph_vec,
                                    const LinearAlgebraDistribution* dist_pt,
                                    HYPRE_IJVector& hypre_ij_vector,
                                    HYPRE_ParVector& hypre_par_vector);

    /// Helper function to create an empty HYPRE_IJVector and
    /// HYPRE_ParVector.
    /// + If no MPI then serial vectors are created
    /// + If MPI and serial distribution then distributed hypre vectors are
    ///   created
    /// + If MPI and distributed input distribution the distributed output
    ///   vectors are created.
    void create_HYPRE_Vector(const LinearAlgebraDistribution* oomph_vec,
                             HYPRE_IJVector& hypre_ij_vector,
                             HYPRE_ParVector& hypre_par_vector);

    /// Helper function to create a serial HYPRE_IJMatrix and
    /// HYPRE_ParCSRMatrix from a CRDoubleMatrix
    void create_HYPRE_Matrix(CRDoubleMatrix* oomph_matrix,
                             HYPRE_IJMatrix& hypre_ij_matrix,
                             HYPRE_ParCSRMatrix& hypre_par_matrix,
                             LinearAlgebraDistribution* dist_pt);

    /// Helper function to set Euclid options using a command line
    /// like array.
    void euclid_settings_helper(const bool& use_block_jacobi,
                                const bool& use_row_scaling,
                                const bool& use_ilut,
                                const int& level,
                                const double& drop_tol,
                                const int& print_level,
                                HYPRE_Solver& euclid_object);

  } // namespace HypreHelpers


  //=====================================================================
  /// An interface class to the suite of Hypre solvers and preconditioners
  /// to allow use of:
  ///
  ///   BoomerAMG (AMG), CG, GMRES or BiCGStab, Euclid (ILU) or
  ///   ParaSails (Approximate inverse)
  ///
  /// Hypre's Krylov subspace solvers (CG, GMRES and BiCGStab) may be
  /// preconditioned using:
  ///
  ///    	BoomerAMG, Euclid or ParaSails
  ///
  //====================================================================
  class HypreInterface
  {
  public:
    ///  Constructor
    HypreInterface()
    {
#ifdef PARANOID
#ifndef HYPRE_SEQUENTIAL
      // For the MPI version of Hypre, check MPI_Helpers::setup has been called
      if (MPI_Helpers::mpi_has_been_initialised() == false)
      {
        std::ostringstream error_message;
        error_message << "When using the MPI version of Hypre please first\n"
                      << "call function MPI_Helpers::setup()\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
#endif

      // setup the distribution
      Hypre_distribution_pt = new LinearAlgebraDistribution();

      // These keep track of which solver and preconditioner
      // (if any) currently exist
      Existing_solver = None;
      Existing_preconditioner = None;

      // Do we want to output info and results of timings?
      Output_info = true;

      // General control paramaters
      Tolerance = 1e-10;
      Max_iter = 100;
      Hypre_method = GMRES;
      Internal_preconditioner = None;

      // Default AMG control parameters -- these seem OK;
      // see hypre documenation for details.
      AMG_using_simple_smoothing = true;
      if (MPI_Helpers::communicator_pt()->nproc() > 1)
      {
        // Jacobi in parallel
        AMG_simple_smoother = 0;
      }
      else
      {
        // Gauss Seidel in serial
        AMG_simple_smoother = 1;
      }
      AMG_complex_smoother = 6;
      AMG_smoother_iterations = 2;
      AMG_damping = 1.0;
      AMG_strength = HypreHelpers::AMG_strength; // 0.25;
      AMG_truncation = HypreHelpers::AMG_truncation; // 0.0;
      AMG_coarsening = HypreHelpers::AMG_coarsening; // 0
      AMG_max_levels = 100;
      AMG_max_row_sum = 1.0;
      AMG_print_level = 0;

      // Parameters for using Euclicd as an AMG smoother (defaults copied
      // from the normal defaults listed in the manual).
      AMGEuclidSmoother_use_block_jacobi = false;
      AMGEuclidSmoother_use_row_scaling = false;
      AMGEuclidSmoother_use_ilut = false;
      AMGEuclidSmoother_level = 1;
      AMGEuclidSmoother_drop_tol = 0; // No dropping
      AMGEuclidSmoother_print_level = 0;

      // Print level for CG, GMRES and BiCGStab
      Krylov_print_level = 0;

      // Set ParaSails default values
      ParaSails_symmetry = 2;
      ParaSails_nlevel = 1;
      ParaSails_thresh = 0.1;
      ParaSails_filter = 0.1;

      // Set Euclid default values
      Euclid_droptol = 0.0;
      Euclid_rowScale = false;
      Euclid_using_ILUT = false;
      Euclid_using_BJ = false;
      Euclid_level = 1;
      Euclid_print_level = 0;

      // Set to true to periodically check the hypre error flag
      // and output any messages
      Hypre_error_messages = false;
    }

    /// Destructor.
    ~HypreInterface()
    {
      // call function to delete solver data
      hypre_clean_up_memory();

      // delete teh oomph-lib distribution
      delete Hypre_distribution_pt;
    }

    /// Broken copy constructor.
    HypreInterface(const HypreInterface&) = delete;

    /// Broken assignment operator.
    void operator=(const HypreInterface&) = delete;

    /// Turn on  the hypre error messages
    void enable_hypre_error_messages()
    {
      Hypre_error_messages = true;
    }

    /// Turn off hypre error messages
    void disable_hypre_error_messages()
    {
      Hypre_error_messages = false;
    }

    /// Enumerated flag to define which Hypre methods are used
    /// CAREFUL: DON'T CHANGE THE ORDER OF THESE!
    enum Hypre_method_types
    {
      CG,
      GMRES,
      BiCGStab,
      BoomerAMG,
      Euclid,
      ParaSails,
      None
    };

    /// Function to return value of which solver (if any) is currently stored.
    unsigned existing_solver()
    {
      return Existing_solver;
    }

    /// Function return value of which preconditioner (if any) is stored.
    unsigned existing_preconditioner()
    {
      return Existing_preconditioner;
    }

    //??ds comment, write access functions and move to protected?
    bool AMGEuclidSmoother_use_block_jacobi;
    bool AMGEuclidSmoother_use_row_scaling;
    bool AMGEuclidSmoother_use_ilut;
    unsigned AMGEuclidSmoother_level;
    double AMGEuclidSmoother_drop_tol;
    unsigned AMGEuclidSmoother_print_level;

  protected:
    /// Function deletes all solver data.
    void hypre_clean_up_memory();

    /// Function which sets values of First_global_row,
    /// Last_global_row and other partitioning data and creates the distributed
    /// Hypre matrix (stored in Matrix_ij/Matrix_par) from the CRDoubleMatrix.
    void hypre_matrix_setup(CRDoubleMatrix* matrix_pt);

    /// Sets up the data required for to use as an oomph-lib
    /// LinearSolver or Preconditioner. This must be called after
    /// the Hypre matrix has been generated using hypre_matrix_setup(...).
    void hypre_solver_setup();

    /// Helper function performs a solve if any solver
    /// exists.
    void hypre_solve(const DoubleVector& rhs, DoubleVector& solution);

    /// Flag is true to output info and results of timings
    bool Output_info;

    /// Maximum number of iterations used in solver.
    unsigned Max_iter;

    /// Tolerance used to terminate solver.
    double Tolerance;

    /// Hypre method flag. Valid values are specified in enumeration
    unsigned Hypre_method;

    /// Preconditioner method flag used with Hypre's PCG,
    /// GMRES and BiCGStab in solve(...) or resolve(...). Valid values
    /// are BoomerAMG, Euclid, ParaSails or None (all enumerated above),
    /// for any other value no preconditioner is set.
    unsigned Internal_preconditioner;

    /// Used to set the Hypre printing level for AMG
    /// 0: no printout
    /// 1: print setup information
    /// 2: print solve information
    /// 3: print setup and solve information
    unsigned AMG_print_level;

    /// Maximum number of levels used in AMG
    unsigned AMG_max_levels;

    /// Parameter to identify diagonally dominant parts of the matrix in AMG
    double AMG_max_row_sum;

    /// Flag to determine whether simple smoothers (determined by the
    ///  AMG_simple_smoother flag) or complex smoothers (determined by the
    ///  AMG_complex_smoother flag are used in AMG
    bool AMG_using_simple_smoothing;

    /// Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    ///  2 = Gauss-Seidel, interior points in parallel, boundary sequential
    ///      (slow in parallel!)
    ///  3 = hybrid Gauss-Seidel or SOR, forward solve
    ///  4 = hybrid Gauss-Seidel or SOR, backward solve
    ///  6 = hybrid symmetric Gauss-Seidel or SSOR
    /// To use these methods set AMG_using_simple_smoothing to true
    unsigned AMG_simple_smoother;

    /// Complex smoothing methods used in BoomerAMG. Relaxation types
    /// are:
    ///  6 = Schwarz
    ///  7 = Pilut
    ///  8 = ParaSails
    ///  9 = Euclid
    /// To use these methods set AMG_using_simple_smoothing to false
    unsigned AMG_complex_smoother;

    /// The number of smoother iterations to apply
    unsigned AMG_smoother_iterations;

    /// Damping factor for BoomerAMG smoothed Jacobi or hybrid SOR
    double AMG_damping;

    /// Connection strength threshold parameter for BoomerAMG
    double AMG_strength;

    /// Interpolation truncation factor for BoomerAMG
    double AMG_truncation;

    /// AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    ///  3 = modified RS with 3rd pass to add C points on the boundaries
    ///  6 = Falgout (uses 1 then CLJP using interior coarse points as
    ///      first independent set)
    ///  8 = PMIS (parallel coarsening using independent sets - lower
    ///      complexities than 0, maybe also slower convergence)
    ///  10= HMIS (one pass RS on each processor then PMIS on interior
    ///      coarse points as first independent set)
    ///  11= One pass RS on each processor (not recommended)
    unsigned AMG_coarsening;

    /// ParaSails symmetry flag, used to inform ParaSails of
    /// Symmetry of definitenss of problem and type of ParaSails
    /// preconditioner:
    /// 0 = nonsymmetric and/or indefinite problem, nonsymmetric preconditioner
    /// 1 = SPD problem, and SPD (factored preconditioner)
    /// 2 = nonsymmetric, definite problem and SDP (factored preconditoner)
    int ParaSails_symmetry;

    /// ParaSails nlevel parameter
    int ParaSails_nlevel;

    /// ParaSails thresh parameter
    double ParaSails_thresh;

    /// ParaSails filter parameter
    double ParaSails_filter;

    /// Euclid drop tolerance for ILU(k) and ILUT factorization
    double Euclid_droptol;

    /// Flag to switch on Euclid row scaling
    bool Euclid_rowScale;

    /// Flag to determine if ILUT (if true) or ILU(k) is used in Euclid
    bool Euclid_using_ILUT;

    /// Flag to determine if Block Jacobi is used instead of PILU
    bool Euclid_using_BJ;

    /// Euclid level parameter for ILU(k) factorization
    int Euclid_level;

    /// Flag to set the level of printing from Euclid
    /// when the Euclid destructor is called
    /// 0: no printing (default)
    /// 1: prints summary of runtime settings and timings
    /// 2: as 1 plus prints memory usage
    unsigned Euclid_print_level;

    /// Used to set the Hypre printing level for the Krylov
    /// subspace solvers
    unsigned Krylov_print_level;

    /// Flag to determine if non-zero values of the Hypre error flag
    /// plus Hypre error messages are output to screen at various points
    /// in the solve function, i.e. after:
    /// 1. setting up the Hypre matrix
    /// 2. setting up the preconditioner
    /// 3. setting up the solver
    /// 4. setting up the Hypre vectors
    /// 5. solving
    /// 6. getting the solution vector
    /// 7. deallocation of solver data
    bool Hypre_error_messages;

    /// Internal flag which is true when hypre_setup or hypre_solve
    /// can delete input matrix.
    bool Delete_input_data;

#ifdef OOMPH_HAS_MPI
    /// Internal flag which tell the solver if the rhs Vector is
    /// distributed or not
    bool Using_distributed_rhs;

    /// Internal flag which tell the solver if the solution Vector to
    /// be returned is distributed or not
    bool Returning_distributed_solution;
#endif

  private:
    /// Hypre copies matrix data from oomph-lib's CRDoubleMatrix
    /// into its own data structures, doubling the memory requirements.
    /// As far as the Hypre solvers are concerned the oomph-lib matrix
    /// is no longer required and could be deleted to save memory.
    /// However, this will not be the expected behaviour and therefore
    /// needs to be requested explicitly by the user by changing this
    /// flag from false (its default) to true.
    bool Delete_matrix;

    /// The Hypre_IJMatrix version of the matrix used in solve(...),
    /// resolve(...) or preconditioner_solve(...).
    HYPRE_IJMatrix Matrix_ij;

    /// The Hypre_ParCSRMatrix version of the matrix used in solve(...),
    /// resolve(...) or preconditioner_solve(...).
    HYPRE_ParCSRMatrix Matrix_par;

    /// The Hypre solver used in solve(...), resolve(...) or
    /// preconditioner_solve(...). [This is a C structure!]
    HYPRE_Solver Solver;

    /// The internal Hypre preconditioner used in conjunction with
    /// Solver. [This is a C structure!]
    HYPRE_Solver Preconditioner;

    /// Used to keep track of which solver (if any) is currently stored.
    unsigned Existing_solver;

    /// Used to keep track of which preconditioner (if any) is currently stored.
    unsigned Existing_preconditioner;

    /// the distribution for this helpers-
    LinearAlgebraDistribution* Hypre_distribution_pt;
  };


  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////

  //=====================================================================
  /// An LinearSolver class using the suite of Hypre solvers
  /// to allow
  ///
  ///   BoomerAMG (AMG), CG, GMRES or BiCGStab
  ///
  /// to be used for LinearSolver::solve(...) or LinearSolver::resolve(...).
  ///
  /// The Krylov subspace solvers (CG, GMRES and BiCGStab) may be
  /// preconditioned using:
  ///
  ///    	BoomerAMG, Euclid or ParaSails
  ///
  /// By default GMRES without preconditioning is used.
  //====================================================================
  class HypreSolver : public LinearSolver, public HypreInterface
  {
  public:
    ///  Constructor
    HypreSolver()
    {
      // Hypre copies matrix data from oomph-lib's CRDoubleMatrix
      // and Distributed CRDoubleMatrix into its own data structures,
      // doubling the memory requirements for the matrix.
      // As far as the Hypre solvers are concerned the oomph-lib matrix
      // is no longer required and could be deleted to save memory.
      // However, this will not be the expected behaviour and therefore
      // needs to be requested explicitly by the user by changing this
      // flag from false (its default) to true.
      Delete_matrix = false;

      // Do we want to output results of timings?
      Doc_time = true;
    }

    /// Empty destructor.
    ~HypreSolver() {}

    /// Broken copy constructor.
    HypreSolver(const HypreSolver&) = delete;

    /// Broken assignment operator.
    void operator=(const HypreSolver&) = delete;

    /// Disable resolve function (overloads the LinearSolver
    /// disable_resolve function).
    void disable_resolve()
    {
      Enable_resolve = false;
      clean_up_memory();
    }

    /// Access function to Max_iter
    unsigned& max_iter()
    {
      return Max_iter;
    }

    /// Access function to Tolerance
    double& tolerance()
    {
      return Tolerance;
    }

    /// Access function to Hypre_method flag -- specified via enumeration.
    unsigned& hypre_method()
    {
      return Hypre_method;
    }

    /// Access function to Internal_preconditioner flag -- specified
    /// via enumeration.
    unsigned& internal_preconditioner()
    {
      return Internal_preconditioner;
    }

    /// Function to select use of 'simple' AMG smoothers as controlled
    /// by AMG_simple_smoother flag
    void amg_using_simple_smoothing()
    {
      AMG_using_simple_smoothing = true;
    }

    /// Access function to AMG_simple_smoother flag
    unsigned& amg_simple_smoother()
    {
      return AMG_simple_smoother;
    }

    /// Function to select use of 'complex' AMG smoothers as controlled
    /// by AMG_complex_smoother flag
    void amg_using_complex_smoothing()
    {
      AMG_using_simple_smoothing = false;
    }

    /// Access function to AMG_complex_smoother flag
    unsigned& amg_complex_smoother()
    {
      return AMG_complex_smoother;
    }

    /// Access function to AMG_print_level
    unsigned& amg_print_level()
    {
      return AMG_print_level;
    }

    /// Access function to AMG_smoother_iterations
    unsigned& amg_smoother_iterations()
    {
      return AMG_smoother_iterations;
    }

    /// Access function to AMG_coarsening flag
    unsigned& amg_coarsening()
    {
      return AMG_coarsening;
    }

    /// Access function to AMG_max_levels
    unsigned& amg_max_levels()
    {
      return AMG_max_levels;
    }

    /// Access function to AMG_damping parameter
    double& amg_damping()
    {
      return AMG_damping;
    }

    /// Access function to AMG_strength
    double& amg_strength()
    {
      return AMG_strength;
    }

    /// Access function to AMG_max_row_sum
    double& amg_max_row_sum()
    {
      return AMG_max_row_sum;
    }

    /// Access function to AMG_truncation
    double& amg_truncation()
    {
      return AMG_truncation;
    }

    /// Access function to ParaSails symmetry flag
    int& parasails_symmetry()
    {
      return ParaSails_symmetry;
    }

    /// Access function to ParaSails nlevel parameter
    int& parasails_nlevel()
    {
      return ParaSails_nlevel;
    }

    /// Access function to ParaSails thresh parameter
    double& parasails_thresh()
    {
      return ParaSails_thresh;
    }

    /// Access function to ParaSails filter parameter
    double& parasails_filter()
    {
      return ParaSails_filter;
    }

    /// Access function to Euclid drop tolerance parameter
    double& euclid_droptol()
    {
      return Euclid_droptol;
    }

    /// Access function to Euclid level parameter
    /// for ILU(k) factorization
    int& euclid_level()
    {
      return Euclid_level;
    }

    /// Enable euclid rowScaling
    void enable_euclid_rowScale()
    {
      Euclid_rowScale = true;
    }

    /// Disable euclid row scaling
    void disable_euclid_rowScale()
    {
      Euclid_rowScale = false;
    }

    /// Enable use of Block Jacobi
    /// as opposed to PILU
    void enable_euclid_using_BJ()
    {
      Euclid_using_BJ = true;
    }

    /// Disable use of Block Jacobi,
    /// so PILU will be used
    void disable_euclid_using_BJ()
    {
      Euclid_using_BJ = false;
    }

    /// Function to switch on ILU(k) factorization for Euclid
    /// (default is ILU(k) factorization)
    void euclid_using_ILUK()
    {
      Euclid_using_ILUT = false;
    }

    /// Function to switch on ILUT factorization for Euclid
    /// (default is ILU(k) factorization)
    void euclid_using_ILUT()
    {
      Euclid_using_ILUT = true;
    }

    /// Function to set the level of printing from Euclid
    /// when the Euclid destructor is called
    /// 0: no printing (default)
    /// 1: prints summary of runtime settings and timings
    /// 2: as 1 plus prints memory usage
    unsigned& euclid_print_level()
    {
      return Euclid_print_level;
    }

    /// Access function to Krylov_print_level
    unsigned& krylov_print_level()
    {
      return Krylov_print_level;
    }

    /// Hypre copies matrix data from oomph-lib's CRDoubleMatrix
    /// or DistributedCRDoubleMatrix into its own data structures,
    /// doubling the memory requirements for the matrix.
    /// As far as the Hypre solvers are concerned the oomph-lib matrix
    /// is no longer required and could be deleted to save memory.
    /// However, this will not be the expected behaviour and therefore
    /// needs to be requested explicitly by the user by calling this function
    /// which changes the
    /// flag from false (its default) to true.
    void enable_delete_matrix()
    {
      Delete_matrix = true;
    }

    /// Hypre copies matrix data from oomph-lib's CRDoubleMatrix
    /// or DistributedCRDoubleMatrix into its own data structures,
    /// doubling the memory requirements for the matrix.
    /// Calling this function ensures that the matrix is not deleted.
    void disable_delete_matrix()
    {
      Delete_matrix = false;
    }


    /// Function which uses problem_pt's get_jacobian(...) function to
    /// generate a linear system which is then solved. This function deletes
    /// any existing internal data and then generates a new Hypre solver.
    void solve(Problem* const& problem_pt, DoubleVector& solution);

    /// Function to solve the linear system defined by matrix_pt
    /// and rhs. This function will delete any existing internal data
    /// and generate a new Hypre solver.
    /// \b Note: The matrix has to be of type CRDoubleMatrix or
    /// Distributed CRDoubleMatrix.
    /// \b Note: Hypre copies matrix data from oomph-lib's CRDoubleMatrix
    /// or DistributedCRDoubleMatrix into its own data structures,
    /// doubling the memory requirements for the matrix.
    /// As far as the Hypre solvers are concerned, the oomph-lib matrix
    /// is no longer required and could be deleted to save memory.
    /// However, this will not be the expected behaviour and therefore
    /// needs to be requested explicitly by the user by calling the
    /// enable_delete_matrix() function.
    void solve(DoubleMatrixBase* const& matrix_pt,
               const DoubleVector& rhs,
               DoubleVector& solution);

    /// Function to resolve a linear system using the existing solver
    /// data, allowing a solve with a new right hand side vector. This
    /// function must be used after a call to solve(...) with
    /// enable_resolve set to true.
    void resolve(const DoubleVector& rhs, DoubleVector& solution);

    /// Function deletes all solver data.
    void clean_up_memory();

  private:
    /// Hypre copies matrix data from oomph-lib's CRDoubleMatrix
    /// or DistributedCRDoubleMatrix into its own data structures,
    /// doubling the memory requirements for the matrix.
    /// As far as the Hypre solvers are concerned the oomph-lib matrix
    /// is no longer required and could be deleted to save memory.
    /// However, this will not be the expected behaviour and therefore
    /// needs to be requested explicitly by the user by changing this
    /// flag from false (its default) to true.
    bool Delete_matrix;
  };


  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////

  //=====================================================================
  /// An Preconditioner class using the suite of Hypre preconditioners
  /// to allow
  ///
  ///      BoomerAMG (Algebraic Multi Grid),
  ///      Euclid (ILU) or
  ///      ParaSails (Approximate inverse)
  ///
  /// to be used for Preconditioner::preconditioner_solve(...).
  /// By default BoomerAMG is used.
  //====================================================================
  class HyprePreconditioner : public Preconditioner, public HypreInterface
  {
  public:
    ///  Constructor. Provide optional string that is used
    /// in annotation of performance
    HyprePreconditioner(const std::string& context_string = "")
    {
      Context_string = context_string;

      // Hypre copies matrix data from oomph-lib's CRDoubleMatrix
      // or DistributedCRDoubleMatrix into its own data structures,
      // doubling the memory requirements for the matrix.
      // As far as the Hypre solvers are concerned the oomph-lib matrix
      // is no longer required and could be deleted to save memory.
      // However, this will not be the expected behaviour and therefore
      // needs to be requested explicitly by the user by changing this
      // flag from false (its default) to true.
      Delete_matrix = false;

      // Do we want to output results of timings?
      Doc_time = true;

      // General control paramaters for using HYPRE as preconditioner
      Tolerance = 0.0;
      Max_iter = 1;
      Hypre_method = BoomerAMG;
      Internal_preconditioner = None;

      // Initialise private double that accumulates the preconditioner
      // solve time of thi instantiation of this class. Is reported
      // in destructor if Report_my_cumulative_preconditioner_solve_time
      // is set to true
      My_cumulative_preconditioner_solve_time = 0.0;
      Report_my_cumulative_preconditioner_solve_time = false;
    }

    /// Destructor.
    ~HyprePreconditioner()
    {
      if (Report_my_cumulative_preconditioner_solve_time)
      {
        oomph_info << "~HyprePreconditioner() in context \" " << Context_string
                   << "\": My_cumulative_preconditioner_solve_time = "
                   << My_cumulative_preconditioner_solve_time << std::endl;
      }
    }

    /// Broken copy constructor.
    HyprePreconditioner(const HyprePreconditioner&) = delete;

    /// Broken assignment operator.
    void operator=(const HyprePreconditioner&) = delete;

    /// Static double that accumulates the preconditioner
    /// solve time of all instantiations of this class. Reset
    /// it manually, e.g. after every Newton solve, using
    /// reset_cumulative_solve_times().
    static double Cumulative_preconditioner_solve_time;

    /// Static double that accumulates the preconditioner
    /// solve time of all instantiations of this class, labeled by
    /// context string. Reset
    /// it manually, e.g. after every Newton solve, using
    /// reset_cumulative_solve_times().
    static std::map<std::string, double> Context_based_cumulative_solve_time;


    /// Static unsigned that accumulates the number of preconditioner
    /// solves of all instantiations of this class. Reset
    /// it manually, e.g. after every Newton solve, using
    /// reset_cumulative_solve_times().
    static unsigned Cumulative_npreconditioner_solve;

    /// Static unsigned that accumulates the number of preconditioner
    /// solves of all instantiations of this class, labeled by
    /// context string. Reset
    /// it manually, e.g. after every Newton solve, using
    /// reset_cumulative_solve_times().
    static std::map<std::string, unsigned>
      Context_based_cumulative_npreconditioner_solve;

    /// Static unsigned that stores nrow for the most recent
    /// instantiations of this class, labeled by
    /// context string. Reset
    /// it manually, e.g. after every Newton solve, using
    /// reset_cumulative_solve_times().
    static std::map<std::string, unsigned> Context_based_nrow;

    /// Report cumulative solve times of all instantiations of this
    /// class
    static void report_cumulative_solve_times();

    /// Reset cumulative solve times
    static void reset_cumulative_solve_times();

    /// Enable reporting of cumulative solve time in destructor
    void enable_report_my_cumulative_preconditioner_solve_time()
    {
      Report_my_cumulative_preconditioner_solve_time = true;
    }

    /// Disable reporting of cumulative solve time in destructor
    void disable_report_my_cumulative_preconditioner_solve_time()
    {
      Report_my_cumulative_preconditioner_solve_time = false;
    }

    /// Enable documentation of preconditioner timings
    void enable_doc_time()
    {
      Doc_time = true;
    }

    /// Disable documentation of preconditioner timings
    void disable_doc_time()
    {
      Doc_time = false;
    }

    /// Access function to Hypre_method flag -- specified via enumeration.
    unsigned& hypre_method()
    {
      return Hypre_method;
    }

    /// Access function to Internal_preconditioner flag -- specified
    /// via enumeration.
    unsigned& internal_preconditioner()
    {
      return Internal_preconditioner;
    }

    /// Function to select BoomerAMG as the preconditioner
    void use_BoomerAMG()
    {
      Hypre_method = BoomerAMG;
    }

    /// Function to set the number of times to apply BoomerAMG
    void set_amg_iterations(const unsigned& amg_iterations)
    {
      Max_iter = amg_iterations;
    }

    /// Return function for Max_iter
    unsigned& amg_iterations()
    {
      return Max_iter;
    }

    /// Function to select use of 'simple' AMG smoothers as controlled
    /// by the flag AMG_simple_smoother
    void amg_using_simple_smoothing()
    {
      AMG_using_simple_smoothing = true;
    }

    /// Access function to AMG_simple_smoother flag
    unsigned& amg_simple_smoother()
    {
      return AMG_simple_smoother;
    }

    /// Function to select use of 'complex' AMG smoothers as controlled
    /// by the flag AMG_complex_smoother
    void amg_using_complex_smoothing()
    {
      AMG_using_simple_smoothing = false;
    }

    /// Access function to AMG_complex_smoother flag
    unsigned& amg_complex_smoother()
    {
      return AMG_complex_smoother;
    }

    /// Return function for the AMG_using_simple_smoothing_flag
    bool& amg_using_simple_smoothing_flag()
    {
      return AMG_using_simple_smoothing;
    }

    /// Access function to AMG_print_level
    unsigned& amg_print_level()
    {
      return AMG_print_level;
    }

    /// Access function to AMG_smoother_iterations
    unsigned& amg_smoother_iterations()
    {
      return AMG_smoother_iterations;
    }

    /// Access function to AMG_coarsening flag
    unsigned& amg_coarsening()
    {
      return AMG_coarsening;
    }

    /// Access function to AMG_max_levels
    unsigned& amg_max_levels()
    {
      return AMG_max_levels;
    }

    /// Access function to AMG_damping parameter
    double& amg_damping()
    {
      return AMG_damping;
    }

    /// Access function to AMG_strength
    double& amg_strength()
    {
      return AMG_strength;
    }

    /// Access function to AMG_max_row_sum
    double& amg_max_row_sum()
    {
      return AMG_max_row_sum;
    }

    /// Access function to AMG_truncation
    double& amg_truncation()
    {
      return AMG_truncation;
    }

    /// Function to select ParaSails as the preconditioner
    void use_ParaSails()
    {
      Hypre_method = ParaSails;
    }

    /// Access function to ParaSails symmetry flag
    int& parasails_symmetry()
    {
      return ParaSails_symmetry;
    }

    /// Access function to ParaSails nlevel parameter
    int& parasails_nlevel()
    {
      return ParaSails_nlevel;
    }

    /// Access function to ParaSails thresh parameter
    double& parasails_thresh()
    {
      return ParaSails_thresh;
    }

    /// Access function to ParaSails filter parameter
    double& parasails_filter()
    {
      return ParaSails_filter;
    }

    /// Function to select use Euclid as the preconditioner
    void use_Euclid()
    {
      Hypre_method = Euclid;
    }

    /// Access function to Euclid drop tolerance parameter
    double& euclid_droptol()
    {
      return Euclid_droptol;
    }

    /// Access function to Euclid level parameter for ILU(k) factorization
    int& euclid_level()
    {
      return Euclid_level;
    }

    /// Enable euclid rowScaling
    void enable_euclid_rowScale()
    {
      Euclid_rowScale = true;
    }

    /// Disable euclid row scaling
    void disable_euclid_rowScale()
    {
      Euclid_rowScale = false;
    }

    /// Enable use of Block Jacobi
    /// as opposed to PILU
    void enable_euclid_using_BJ()
    {
      Euclid_using_BJ = true;
    }

    /// Disable use of Block Jacobi,
    /// so PILU will be used
    void disable_euclid_using_BJ()
    {
      Euclid_using_BJ = false;
    }

    /// Function to switch on ILU(k) factorization for Euclid
    /// (default is ILU(k) factorization)
    void euclid_using_ILUK()
    {
      Euclid_using_ILUT = false;
    }

    /// Function to switch on ILUT factorization for Euclid
    /// (default is ILU(k) factorization)
    void euclid_using_ILUT()
    {
      Euclid_using_ILUT = true;
    }

    /// Function to set the level of printing from Euclid
    /// when the Euclid destructor is called
    /// 0: no printing (default)
    /// 1: prints summary of runtime settings and timings
    /// 2: as 1 plus prints memory usage
    unsigned& euclid_print_level()
    {
      return Euclid_print_level;
    }

    /// Hypre copies matrix data from oomph-lib's CRDoubleMatrix
    /// or DistributedCRDoubleMatrix into its own data structures,
    /// doubling the memory requirements for the matrix.
    /// As far as the Hypre solvers are concerned the oomph-lib matrix
    /// is no longer required and could be deleted to save memory.
    /// However, this will not be the expected behaviour and therefore
    /// needs to be requested explicitly by the user by calling this function
    /// which changes the
    /// flag from false (its default) to true.
    void enable_delete_matrix()
    {
      Delete_matrix = true;
    }

    /// Hypre copies matrix data from oomph-lib's CRDoubleMatrix
    /// or DistributedCRDoubleMatrix into its own data structures,
    /// doubling the memory requirements for the matrix.
    /// Calling this function ensures that the matrix is not deleted.
    void disable_delete_matrix()
    {
      Delete_matrix = false;
    }

    /// Function to set up a preconditioner for the linear
    /// system defined by matrix_pt. This function is required when
    /// preconditioning and must be called before using the
    /// preconditioner_solve(...) function. This interface allows
    /// HyprePreconditioner to be used as a Preconditioner object
    /// for oomph-lib's own IterativeLinearSolver class.
    /// \b Note: matrix_pt must point to an object of type
    /// CRDoubleMatrix or DistributedCRDoubleMatrix.
    /// \b Note: Hypre copies matrix data from oomph-lib's CRDoubleMatrix
    /// and DistributedCRDoubleMatrix into its own data structures,
    /// doubling the memory requirements for the matrix.
    /// As far as the Hypre solvers are concerned, the oomph-lib matrix
    /// is no longer required and could be deleted to save memory.
    /// However, this will not be the expected behaviour and therefore
    /// needs to be requested explicitly by the user by calling the
    /// enable_delete_matrix() function.
    void setup();

    /// Function applies solver to vector r for preconditioning.
    /// This requires a call to setup(...) first.
    /// \b Note: Hypre copies matrix data from oomph-lib's CRDoubleMatrix
    /// or DistributedCRDoubleMatrix into its own data structures,
    /// doubling the memory requirements for the matrix.
    /// As far as the Hypre solvers are concerned, the oomph-lib matrix
    /// is no longer required and could be deleted to save memory.
    /// However, this will not be the expected behaviour and therefore
    /// needs to be requested explicitly by the user with the
    /// delete_matrix() function.
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    /// Function deletes all solver data.
    void clean_up_memory();

  private:
    /// Hypre copies matrix data from oomph-lib's CRDoubleMatrix
    /// or DistributedCRDoubleMatrix into its own data structures,
    /// doubling the memory requirements for the matrix.
    /// As far as the Hypre solvers are concerned the oomph-lib matrix
    /// is no longer required and could be deleted to save memory.
    /// However, this will not be the expected behaviour and therefore
    /// needs to be requested explicitly by the user by changing this
    /// flag from false (its default) to true.
    bool Delete_matrix;

    // Flag is true to output results of timings
    bool Doc_time;

    /// Private double that accumulates the preconditioner
    /// solve time of thi instantiation of this class. Is reported
    /// in destructor if Report_my_cumulative_preconditioner_solve_time
    /// is set to true
    double My_cumulative_preconditioner_solve_time;

    /// Bool to request report of My_cumulative_preconditioner_solve_time
    /// in destructor
    bool Report_my_cumulative_preconditioner_solve_time;

    /// String can be used to provide context for annotation
    std::string Context_string;
  };


  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////

  //==================================================================
  /// Default settings for various uses of the Hypre solver
  //==================================================================
  namespace Hypre_default_settings
  {
    /// Set default parameters for use as preconditioner in
    /// for momentum block in Navier-Stokes problem
    extern void set_defaults_for_navier_stokes_momentum_block(
      HyprePreconditioner* hypre_preconditioner_pt);

    /// Set default parameters for use as preconditioner in
    /// 2D Poisson-type problem.
    extern void set_defaults_for_2D_poisson_problem(
      HyprePreconditioner* hypre_preconditioner_pt);

    /// Set default parameters for use as preconditioner in
    /// 3D Poisson-type problem.
    extern void set_defaults_for_3D_poisson_problem(
      HyprePreconditioner* hypre_preconditioner_pt);

  } // namespace Hypre_default_settings
} // namespace oomph

#endif
