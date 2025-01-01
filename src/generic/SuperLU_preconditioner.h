// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_SuperLU_Preconditioner_HEADER
#define OOMPH_SuperLU_Preconditioner_HEADER

// oomph-lib headers

#include "linear_solver.h"
#include "preconditioner.h"

namespace oomph
{
  //====================================================================
  /// An interface to allow SuperLU to be used as an (exact) Preconditioner
  //====================================================================
  class SuperLUPreconditioner : public Preconditioner
  {
  public:
    /// Constructor.
    SuperLUPreconditioner()
    {
      Solver.disable_doc_stats();
      Solver.disable_doc_time();
    }

    /// Destructor.
    ~SuperLUPreconditioner() {}

    /// Broken copy constructor.
    SuperLUPreconditioner(const SuperLUPreconditioner&) = delete;

    /// Broken assignment operator.
    void operator=(const SuperLUPreconditioner&) = delete;

    /// Function to set up a preconditioner for the linear
    /// system defined by matrix_pt. This function must be called
    /// before using preconditioner_solve.
    /// Note: matrix_pt must point to an object of class
    /// CRDoubleMatrix or CCDoubleMatrix
    void setup()
    {
      oomph_info << "Setting up SuperLU (exact) preconditioner" << std::endl;
      if (dynamic_cast<DistributableLinearAlgebraObject*>(matrix_pt()) != 0)
      {
        LinearAlgebraDistribution dist(
          (dynamic_cast<DistributableLinearAlgebraObject*>(matrix_pt()))
            ->distribution_pt());
        this->build_distribution(dist);
        Solver.factorise(matrix_pt());
      }
      else
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "SuperLUPreconditioner can only be applied to matrices derived \n"
          << "DistributableLinearAlgebraObject.\n"
          << "You are most likely to be here because you are using the\n "
          << "soon to be obsolete CCDoubleMatrix\n";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }

    /// Function applies SuperLU to vector r for (exact) preconditioning,
    /// this requires a call to setup(...) first.
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z)
    {
      Solver.resolve(r, z);
    }

    /// Function applies SuperLU to vector r for (exact) preconditioning
    /// (of the transposed matrix system) this requires a call to setup(...)
    /// first.
    void preconditioner_solve_transpose(const DoubleVector& r, DoubleVector& z)
    {
      Solver.resolve_transpose(r, z);
    }


    /// Clean up memory -- forward the call to the version in
    /// SuperLU in its LinearSolver incarnation.
    virtual void clean_up_memory()
    {
      Solver.clean_up_memory();
    }


    /// Get the amount of memory used to store the LU factors inside SuperLU
    double get_memory_usage_for_lu_factors()
    {
      // Return the appropriate result
      return Solver.get_memory_usage_for_lu_factors();
    } // End of get_memory_usage_for_lu_factors


    /// Get the total memory needed by SuperLU to store AND calculate
    /// the LU factors
    double get_total_memory_needed_for_superlu()
    {
      // Return the appropriate result
      return Solver.get_total_needed_memory();
    } // End of get_memory_usage_for_superlu


    /// Get the amount of memory taken up by SuperLU. The first entry
    /// of the returned result contains the memory used to store the LU
    /// factors and the second entry contains the total memory used to
    /// store AND calculate the LU factors
    Vector<double> get_memory_usage_for_superlu()
    {
      // Allocate storage for the memory statistics
      Vector<double> memory_usage(2, 0.0);

      // The first entry contains the memory used to store the LU factors
      memory_usage[0] = Solver.get_memory_usage_for_lu_factors();

      // The second entry contains the total memory used to both calculate
      // and store the LU factors
      memory_usage[1] = Solver.get_total_needed_memory();

      // Now return the calculated result
      return memory_usage;
    } // End of get_memory_usage_for_superlu


    /// Enable documentation of solver statistics
    void enable_doc_stats()
    {
      // Enable the documentation of statistics inside SuperLU
      Solver.enable_doc_stats();
    } // End of enable_doc_stats

    /// Enable documentation of solver statistics
    void disable_doc_stats()
    {
      // Disable the documentation of statistics inside SuperLU
      Solver.disable_doc_stats();
    } // End of disable_doc_stats


  private:
    /// the SuperLU solver emplyed by this preconditioner
    SuperLUSolver Solver;
  };

} // namespace oomph
#endif
