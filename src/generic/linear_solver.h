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

// This header defines a class for linear solvers

// Include guards
#ifndef OOMPH_LINEAR_SOLVER_HEADER
#define OOMPH_LINEAR_SOLVER_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

// oomph-lib headers
#include "Vector.h"
#include "double_vector.h"
#include "matrices.h"

namespace oomph
{
  // Forward declaration of problem class
  class Problem;

  //====================================================================
  /// Base class for all linear solvers. This merely defines standard
  /// interfaces for linear solvers, so that different solvers can be
  /// used in a clean and transparent manner. Note that LinearSolvers
  /// are primarily used to solve the linear systems arising in
  /// oomph-lib's Newton iteration. Their primary solve function
  /// therefore takes a pointer to the associated problem, construct its
  /// Jacobian matrix and residual vector, and return the solution
  /// of the linear system formed by the Jacobian and the residual vector.
  ///  We also provide broken virtual interfaces
  /// to a linear-algebra-type solve function in which the matrix
  /// and the rhs can be specified, but this are not guaranteed to
  /// implemented for all linear solvers (e.g. for frontal solvers).
  //====================================================================
  class LinearSolver : public DistributableLinearAlgebraObject
  {
  protected:
    /// Boolean that indicates whether the matrix (or its factors, in
    /// the case of direct solver) should be stored so that the resolve function
    /// can be used.
    bool Enable_resolve;

    /// Boolean flag that indicates whether the time taken
    // for the solves should be documented
    bool Doc_time;

    /// flag that indicates whether the gradient required for the
    /// globally convergent Newton method should be computed or not
    bool Compute_gradient;

    /// flag that indicates whether the gradient was computed or not
    bool Gradient_has_been_computed;

    /// DoubleVector storing the gradient for the globally convergent
    /// Newton method
    DoubleVector Gradient_for_glob_conv_newton_solve;

  public:
    /// Empty constructor, initialise the member data
    LinearSolver()
      : Enable_resolve(false),
        Doc_time(true),
        Compute_gradient(false),
        Gradient_has_been_computed(false)
    {
    }

    /// Broken copy constructor
    LinearSolver(const LinearSolver& dummy) = delete;

    /// Broken assignment operator
    void operator=(const LinearSolver&) = delete;

    /// Empty virtual destructor
    virtual ~LinearSolver() {}

    /// Enable documentation of solve times
    void enable_doc_time()
    {
      Doc_time = true;
    }

    /// Disable documentation of solve times
    void disable_doc_time()
    {
      Doc_time = false;
    }

    /// Is documentation of solve times enabled?
    bool is_doc_time_enabled() const
    {
      return Doc_time;
    }

    /// Boolean flag indicating if resolves are enabled
    bool is_resolve_enabled() const
    {
      return Enable_resolve;
    }

    /// Enable resolve (i.e. store matrix and/or LU decomposition, say)
    /// Virtual so it can be overloaded to perform additional tasks
    virtual void enable_resolve()
    {
      Enable_resolve = true;
    }

    /// Disable resolve (i.e. store matrix and/or LU decomposition, say)
    /// This function simply resets an internal flag. It's virtual so
    /// it can be overloaded to perform additional tasks such as
    /// cleaning up memory that is only required for the resolve.
    virtual void disable_resolve()
    {
      Enable_resolve = false;
    }

    /// Solver: Takes pointer to problem and returns the results vector
    /// which contains the solution of the linear system defined by
    /// the problem's fully assembled Jacobian and residual vector.
    virtual void solve(Problem* const& problem_pt, DoubleVector& result) = 0;

    /// Linear-algebra-type solver: Takes pointer to a matrix and rhs
    /// vector and returns the solution of the linear system.
    virtual void solve(DoubleMatrixBase* const& matrix_pt,
                       const DoubleVector& rhs,
                       DoubleVector& result)
    {
      throw OomphLibError(
        "DoubleVector based solve function not implemented for this solver",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Linear-algebra-type solver: Takes pointer to a matrix and rhs
    /// vector and returns the solution of the linear system.
    virtual void solve(DoubleMatrixBase* const& matrix_pt,
                       const Vector<double>& rhs,
                       Vector<double>& result)
    {
      throw OomphLibError(
        "Vector<double> based solve function not implemented for this solver",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Solver: Takes pointer to problem and returns the results vector
    /// which contains the solution of the linear system defined by the
    /// problem's fully assembled Jacobian and residual vector (broken virtual).
    virtual void solve_transpose(Problem* const& problem_pt,
                                 DoubleVector& result)
    {
      // Create an output stream
      std::ostringstream error_message_stream;

      // Create the error message
      error_message_stream << "The function to solve the transposed system has "
                           << "not yet been\nimplemented for this solver."
                           << std::endl;

      // Throw the error message
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    } // End of solve_transpose

    /// Linear-algebra-type solver: Takes pointer to a matrix and rhs
    /// vector and returns the solution of the linear system.
    virtual void solve_transpose(DoubleMatrixBase* const& matrix_pt,
                                 const DoubleVector& rhs,
                                 DoubleVector& result)
    {
      throw OomphLibError(
        "DoubleVector based solve function not implemented for this solver",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Linear-algebra-type solver: Takes pointer to a matrix and rhs
    /// vector and returns the solution of the linear system.
    virtual void solve_transpose(DoubleMatrixBase* const& matrix_pt,
                                 const Vector<double>& rhs,
                                 Vector<double>& result)
    {
      throw OomphLibError(
        "Vector<double> based solve function not implemented for this solver",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Resolve the system defined by the last assembled jacobian
    /// and the rhs vector. Solution is returned in the vector result.
    /// (broken virtual)
    virtual void resolve(const DoubleVector& rhs, DoubleVector& result)
    {
      throw OomphLibError(
        "Resolve function not implemented for this linear solver",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Solver: Resolve the system defined by the last assembled jacobian
    /// and the rhs vector. Solution is returned in the vector result.
    /// (broken virtual)
    virtual void resolve_transpose(const DoubleVector& rhs,
                                   DoubleVector& result)
    {
      // Create an output stream
      std::ostringstream error_message_stream;

      // Create the error message
      error_message_stream
        << "The function to resolve the transposed system has "
        << "not yet been\nimplemented for this solver." << std::endl;

      // Throw the error message
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    } // End of resolve_transpose

    /// Empty virtual function that can be overloaded in specific
    /// linear solvers to clean up any memory that may have been
    /// allocated (e.g. when preparing for a re-solve).
    virtual void clean_up_memory() {}

    ///  returns the time taken to assemble the Jacobian matrix and
    /// residual vector (needs to be overloaded for each solver)
    virtual double jacobian_setup_time() const
    {
      throw OomphLibError(
        "jacobian_setup_time has not been implemented for this linear solver",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
      return 0;
    }

    /// return the time taken to solve the linear system (needs to be
    /// overloaded for each linear solver)
    virtual double linear_solver_solution_time() const
    {
      throw OomphLibError(
        "linear_solver_solution_time() not implemented for this linear solver",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
      return 0;
    }

    /// function to enable the computation of the gradient required
    /// for the globally convergent Newton method
    virtual void enable_computation_of_gradient()
    {
      throw OomphLibError(
        "enable_computation_of_gradient() not implemented for "
        "this linear solver",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// function to disable the computation of the gradient required
    /// for the globally convergent Newton method
    void disable_computation_of_gradient()
    {
      Compute_gradient = false;
    }

    /// function to reset the size of the gradient before each Newton
    /// solve
    void reset_gradient()
    {
      Gradient_for_glob_conv_newton_solve.clear();
    }

    /// function to access the gradient, provided it has been computed
    void get_gradient(DoubleVector& gradient)
    {
#ifdef PARANOID
      if (Gradient_has_been_computed)
      {
#endif
        gradient = Gradient_for_glob_conv_newton_solve;
#ifdef PARANOID
      }
      else
      {
        throw OomphLibError(
          "The gradient has not been computed for this linear solver!",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }
  };

  //=============================================================================
  /// Dense LU decomposition-based solve of full assembled linear system.
  /// VERY inefficient but useful to illustrate the principle.
  /// Only suitable for use with Serial matrices and vectors.
  /// This solver will only work with non-distributed matrices and vectors
  /// (note: DenseDoubleMatrix is not distributable)
  //============================================================================
  class DenseLU : public LinearSolver
  {
    /// The DenseDoubleMatrix class is a friend
    friend class DenseDoubleMatrix;

  public:
    /// Constructor, initialise storage
    DenseLU()
      : Jacobian_setup_time(0.0),
        Solution_time(0.0),
        Sign_of_determinant_of_matrix(0),
        Index(0),
        LU_factors(0)
    {
      // Shut up!
      Doc_time = false;
    }

    /// Broken copy constructor
    DenseLU(const DenseLU& dummy) = delete;

    /// Broken assignment operator
    void operator=(const DenseLU&) = delete;

    /// Destructor, clean up the stored LU factors
    ~DenseLU()
    {
      clean_up_memory();
    }

    /// Solver: Takes pointer to problem and returns the results Vector
    /// which contains the solution of the linear system defined by
    /// the problem's fully assembled Jacobian and residual Vector.
    void solve(Problem* const& problem_pt, DoubleVector& result);

    /// Linear-algebra-type solver: Takes pointer to a matrix and rhs
    /// vector and returns the solution of the linear system.
    void solve(DoubleMatrixBase* const& matrix_pt,
               const DoubleVector& rhs,
               DoubleVector& result);

    /// Linear-algebra-type solver: Takes pointer to a matrix and rhs
    /// vector and returns the solution of the linear system.
    void solve(DoubleMatrixBase* const& matrix_pt,
               const Vector<double>& rhs,
               Vector<double>& result);

    ///  returns the time taken to assemble the jacobian matrix and
    /// residual vector
    double jacobian_setup_time() const
    {
      return Jacobian_setup_time;
    }

    /// return the time taken to solve the linear system (needs to be
    /// overloaded for each linear solver)
    virtual double linear_solver_solution_time() const
    {
      return Solution_time;
    }

  protected:
    /// Perform the LU decomposition of the matrix
    void factorise(DoubleMatrixBase* const& matrix_pt);

    /// Do the backsubstitution step to solve the system LU result = rhs
    void backsub(const DoubleVector& rhs, DoubleVector& result);

    /// perform back substitution using Vector<double>
    void backsub(const Vector<double>& rhs, Vector<double>& result);

    /// Clean up the stored LU factors
    void clean_up_memory();

    /// Jacobian setup time
    double Jacobian_setup_time;

    /// Solution time
    double Solution_time;

    /// Sign of the determinant of the matrix
    /// (obtained during the LU decomposition)
    int Sign_of_determinant_of_matrix;

  private:
    /// Pointer to storage for the index of permutations in the LU solve
    long* Index;

    /// Pointer to storage for the LU decomposition
    double* LU_factors;
  };


  //====================================================================
  /// Dense LU decomposition-based solve of linear system
  /// assembled via finite differencing of the residuals Vector.
  /// Even more inefficient than DenseLU but excellent sanity check!
  //====================================================================
  class FD_LU : public DenseLU
  {
  public:
    /// Constructor: empty
    FD_LU() : DenseLU() {}

    /// Broken copy constructor
    FD_LU(const FD_LU& dummy) = delete;

    /// Broken assignment operator
    void operator=(const FD_LU&) = delete;

    /// Solver: Takes pointer to problem and returns the results Vector
    /// which contains the solution of the linear system defined by
    /// the problem's residual Vector (Jacobian computed by FD approx.)
    void solve(Problem* const& problem_pt, DoubleVector& result);

    /// Linear-algebra-type solver: Takes pointer to a matrix and rhs
    /// vector and returns the solution of the linear system.
    void solve(DoubleMatrixBase* const& matrix_pt,
               const DoubleVector& rhs,
               DoubleVector& result)
    {
      DenseLU::solve(matrix_pt, rhs, result);
    }

    /// Linear-algebra-type solver: Takes pointer to a matrix
    /// and rhs vector and returns the solution of the linear system
    /// Call the broken base-class version. If you want this, please
    /// implement it
    void solve(DoubleMatrixBase* const& matrix_pt,
               const Vector<double>& rhs,
               Vector<double>& result)
    {
      LinearSolver::solve(matrix_pt, rhs, result);
    }
  };


  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////


  //=============================================================================
  /// SuperLU Project Solver class. This is a combined wrapper for both
  /// SuperLU and SuperLU Dist.
  /// See http://crd.lbl.gov/~xiaoye/SuperLU/
  /// Default Behaviour: If this solver is distributed over more than one
  /// processor then SuperLU Dist is used.
  /// Member data naming convention: member data associated with the SuperLU
  /// Dist solver begins Dist_... and member data associated with the serial
  /// SuperLU solver begins Serial_... .
  //=============================================================================
  class SuperLUSolver : public LinearSolver
  {
  public:
    /// enum type to specify the solver behaviour.
    /// Default - will employ superlu dist if more than 1 processor.
    /// Serial - will always try use superlu (serial).
    /// Distributed - will always try to use superlu dist.
    enum Type
    {
      Default,
      Serial,
      Distributed
    };

    /// Constructor. Set the defaults.
    SuperLUSolver() : Serial_f_factors(0)
    {
      // Set solver wide default values and null pointers
      Doc_stats = false;
      Suppress_solve = false;
      Using_dist = false;
      Solver_type = Default;

#ifdef OOMPH_HAS_MPI
      // Set default values and nullify pointers for SuperLU Dist
      Dist_use_global_solver = false;
      Dist_delete_matrix_data = false;
      Dist_allow_row_and_col_permutations = true;
      Dist_solver_data_pt = 0;
      Dist_global_solve_data_allocated = false;
      Dist_distributed_solve_data_allocated = false;
      Dist_value_pt = 0;
      Dist_index_pt = 0;
      Dist_start_pt = 0;
#endif

      // Set default values and null pointers for SuperLU (serial)
      Serial_compressed_row_flag = true;
      Serial_sign_of_determinant_of_matrix = 0;
      Serial_n_dof = 0;
    }

    /// Broken copy constructor
    SuperLUSolver(const SuperLUSolver& dummy) = delete;

    /// Broken assignment operator
    void operator=(const SuperLUSolver&) = delete;

    /// Destructor, clean up the stored matrices
    ~SuperLUSolver()
    {
      clean_up_memory();
    }

    /// function to enable the computation of the gradient
    void enable_computation_of_gradient()
    {
      Compute_gradient = true;
    }

    // SuperLUSolver methods
    /// /////////////////////

    /// Overload disable resolve so that it cleans up memory too
    void disable_resolve()
    {
      LinearSolver::disable_resolve();
      clean_up_memory();
    }

    /// Solver: Takes pointer to problem and returns the results Vector
    /// which contains the solution of the linear system defined by
    /// the problem's fully assembled Jacobian and residual Vector.
    void solve(Problem* const& problem_pt, DoubleVector& result);

    /// Linear-algebra-type solver: Takes pointer to a matrix and rhs
    /// vector and returns the solution of the linear system.
    /// The function returns the global result Vector.
    /// Note: if Delete_matrix_data is true the function
    /// matrix_pt->clean_up_memory() will be used to wipe the matrix data.
    void solve(DoubleMatrixBase* const& matrix_pt,
               const DoubleVector& rhs,
               DoubleVector& result);


    /*  /// Linear-algebra-type solver: Takes pointer to a matrix */
    /*  /// and rhs vector and returns the solution of the linear system */
    /*  /// Call the broken base-class version. If you want this, please  */
    /*  /// implement it */
    /*  void solve(DoubleMatrixBase* const &matrix_pt, */
    /*                     const Vector<double> &rhs, */
    /*                     Vector<double> &result) */
    /*   {LinearSolver::solve(matrix_pt,rhs,result);} */

    /// Solver: Takes pointer to problem and returns the results Vector
    /// which contains the solution of the linear system defined by
    /// the problem's fully assembled Jacobian and residual Vector.
    void solve_transpose(Problem* const& problem_pt, DoubleVector& result);

    /// Linear-algebra-type solver: Takes pointer to a matrix and rhs
    /// vector and returns the solution of the linear system.
    /// The function returns the global result Vector.
    /// Note: if Delete_matrix_data is true the function
    /// matrix_pt->clean_up_memory() will be used to wipe the matrix data.
    void solve_transpose(DoubleMatrixBase* const& matrix_pt,
                         const DoubleVector& rhs,
                         DoubleVector& result);

    /// Resolve the system defined by the last assembled jacobian
    /// and the specified rhs vector if resolve has been enabled.
    /// Note: returns the global result Vector.
    void resolve(const DoubleVector& rhs, DoubleVector& result);

    /// Resolve the (transposed) system defined by the last assembled
    /// Jacobian and the specified rhs vector if resolve has been enabled.
    void resolve_transpose(const DoubleVector& rhs, DoubleVector& result);

    /// Enable documentation of solver statistics
    void enable_doc_stats()
    {
      Doc_stats = true;
    }

    /// Disable documentation of solver statistics
    void disable_doc_stats()
    {
      Doc_stats = false;
    }

    /// returns the time taken to assemble the jacobian matrix and
    /// residual vector
    double jacobian_setup_time() const
    {
      return Jacobian_setup_time;
    }

    /// return the time taken to solve the linear system (needs to be
    /// overloaded for each linear solver)
    virtual double linear_solver_solution_time() const
    {
      return Solution_time;
    }

    /// Do the factorisation stage
    /// Note: if Delete_matrix_data is true the function
    /// matrix_pt->clean_up_memory() will be used to wipe the matrix data.
    void factorise(DoubleMatrixBase* const& matrix_pt);

    /// Do the backsubstitution for SuperLU solver
    /// Note: returns the global result Vector.
    void backsub(const DoubleVector& rhs, DoubleVector& result);

    /// Do the back-substitution for transposed system of the SuperLU
    /// solver Note: Returns the global result Vector.
    void backsub_transpose(const DoubleVector& rhs, DoubleVector& result);

    /// Clean up the memory allocated by the solver
    void clean_up_memory();

    /// Specify the solve type. Either default, serial or distributed.
    /// See enum SuperLU_solver_type for more details.
    void set_solver_type(const Type& t)
    {
      this->clean_up_memory();
      Solver_type = t;
    }

    // SuperLU (serial) methods
    /// ////////////////////////

    /// Use the compressed row format in superlu serial
    void use_compressed_row_for_superlu_serial()
    {
      Serial_compressed_row_flag = true;
    }

    /// Use the compressed column format in superlu serial
    void use_compressed_column_for_superlu_serial()
    {
      Serial_compressed_row_flag = false;
    }

#ifdef OOMPH_HAS_MPI

    // SuperLU Dist methods
    /// ////////////////////

    /// Set Delete_matrix_data flag. SuperLU_dist needs its own copy
    /// of the input matrix, therefore a copy must be made if any matrix
    /// used with this solver is to be preserved. If the input matrix can be
    /// deleted the flag can be set to true to reduce the amount of memory
    /// required, and the matrix data will be wiped using its clean_up_memory()
    /// function.  Default value is false.
    void enable_delete_matrix_data_in_superlu_dist()
    {
      Dist_delete_matrix_data = true;
    }

    /// Unset Delete_matrix_data flag. SuperLU_dist needs its own copy
    /// of the input matrix, therefore a copy must be made if any matrix
    /// used with this solver is to be preserved. If the input matrix can be
    void disable_delete_matrix_data_in_superlu_dist()
    {
      Dist_delete_matrix_data = false;
    }

    /// Set flag so that SuperLU_DIST is allowed to permute matrix rows
    /// and columns during factorisation. This is the default for SuperLU_DIST,
    /// and can lead to significantly faster solves.
    void enable_row_and_col_permutations_in_superlu_dist()
    {
      Dist_allow_row_and_col_permutations = true;
    }

    /// Set flag so that SuperLU_DIST is not allowed to permute matrix
    /// rows and columns during factorisation.
    void disable_row_and_col_permutations_in_superlu_dist()
    {
      Dist_allow_row_and_col_permutations = false;
    }

    /// Calling this method will ensure that when the problem based
    /// solve interface is used, a global (serial) jacobian will be
    /// assembled.
    /// Note: calling this function will delete any distributed solve data.
    void use_global_solve_in_superlu_dist()
    {
      if (!Dist_use_global_solver)
      {
        clean_up_memory();
        Dist_use_global_solver = true;
      }
    }

    /// Calling this method will ensure that when the problem based
    /// solve interface is used, a distributed jacobian will be
    /// assembled.
    /// Note: calling this function will delete any global solve data.
    void use_distributed_solve_in_superlu_dist()
    {
      if (Dist_use_global_solver)
      {
        clean_up_memory();
        Dist_use_global_solver = false;
      }
    }

#endif

  private:
    /// factorise method for SuperLU (serial)
    void factorise_serial(DoubleMatrixBase* const& matrix_pt);

    /// backsub method for SuperLU (serial)
    void backsub_serial(const DoubleVector& rhs, DoubleVector& result);

    /// backsub method for SuperLU (serial)
    void backsub_transpose_serial(const DoubleVector& rhs,
                                  DoubleVector& result);

#ifdef OOMPH_HAS_MPI
    /// factorise method for SuperLU Dist
    void factorise_distributed(DoubleMatrixBase* const& matrix_pt);

    /// backsub method for SuperLU Dist
    void backsub_distributed(const DoubleVector& rhs, DoubleVector& result);

    /// backsub method for SuperLU Dist
    void backsub_transpose_distributed(const DoubleVector& rhs,
                                       DoubleVector& result);
#endif

    // SuperLUSolver member data
    /// /////////////////////////

    /// Jacobian setup time
    double Jacobian_setup_time;

    /// Solution time
    double Solution_time;

    /// Suppress solve?
    bool Suppress_solve;

    /// Set to true to output statistics (false by default).
    bool Doc_stats;

    /// the solver type. see SuperLU_solver_type for details.
    Type Solver_type;

    /// boolean flag indicating whether superlu dist is being used
    bool Using_dist;

    // SuperLU (serial) member data
    /// ////////////////////////////

    /// Storage for the LU factors as required by SuperLU
    void* Serial_f_factors;

    /// Info flag for the SuperLU solver
    int Serial_info;

    /// The number of unknowns in the linear system
    unsigned long Serial_n_dof;

    /// Sign of the determinant of the matrix
    int Serial_sign_of_determinant_of_matrix;

    /// Use compressed row version?
    bool Serial_compressed_row_flag;

  public:
    /// How much memory do the LU factors take up? In bytes
    double get_memory_usage_for_lu_factors();

    /// How much memory was allocated by SuperLU? In bytes
    double get_total_needed_memory();

  private:
#ifdef OOMPH_HAS_MPI

    // SuperLU Dist member data
    /// ////////////////////////
  public:
    /// Static flag that determines whether the warning about
    /// incorrect distribution of RHSs will be printed or not
    static bool Suppress_incorrect_rhs_distribution_warning_in_resolve;

  private:
    /// Flag that determines whether the MPIProblem based solve function
    /// uses the global or distributed version of SuperLU_DIST
    /// (default value is false).
    bool Dist_use_global_solver;

    ///  Flag is true if solve data has been generated for a global matrix
    bool Dist_global_solve_data_allocated;

    ///  Flag is true if solve data has been generated for distributed matrix
    bool Dist_distributed_solve_data_allocated;

    /// Storage for the LU factors and other data required by SuperLU
    void* Dist_solver_data_pt;

    /// Number of rows for the process grid
    int Dist_nprow;

    /// Number of columns for the process grid
    int Dist_npcol;

    /// Info flag for the SuperLU solver
    int Dist_info;

    /// If true then SuperLU_DIST is allowed to permute matrix rows
    /// and columns during factorisation. This is the default for SuperLU_DIST,
    /// and can lead to significantly faster solves, but has been known to
    /// fail, hence the default value is 0.
    bool Dist_allow_row_and_col_permutations;

    /// Delete_matrix_data flag. SuperLU_dist needs its own copy
    /// of the input matrix, therefore a copy must be made if any matrix
    /// used with this solver is to be preserved. If the input matrix can be
    /// deleted the flag can be set to true to reduce the amount of memory
    /// required, and the matrix data will be wiped using its clean_up_memory()
    /// function. Default value is false.
    bool Dist_delete_matrix_data;

    /// Pointer for storage of the matrix values required by SuperLU_DIST
    double* Dist_value_pt;

    /// Pointer for storage of matrix rows or column indices required
    /// by SuperLU_DIST
    int* Dist_index_pt;

    /// Pointers for storage of matrix column or row starts
    // required by SuperLU_DIST
    int* Dist_start_pt;

#endif
  }; // end of SuperLUSolver
} // namespace oomph
#endif
