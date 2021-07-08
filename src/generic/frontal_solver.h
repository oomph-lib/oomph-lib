// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision$
// LIC//
// LIC// $LastChangedDate$
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
// This is the header file for the C-wrapper functions for the HSL MA42
// frontal solver

// Include guards to prevent multiple inclusions of the header
#ifndef HSL_MA42OOMPH__FRONTAL_SOLVER_HEADER
#define HSL_MA42OOMPH__FRONTAL_SOLVER_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

// oomph-lib headers
#include "Vector.h"
#include "linear_solver.h"

namespace oomph
{
  //====================================================================
  /// Linear solver class that provides a wrapper to the frontal
  /// solver MA42 from the <a href="http://www.hsl.rl.ac.uk/">HSL
  /// library;</a> see  <a href="http://www.hsl.rl.ac.uk/">
  /// http://www.hsl.rl.ac.uk/.</A>
  //====================================================================
  class HSL_MA42 : public LinearSolver
  {
  private:
    /// \short Special solver for problems with 1 dof (MA42 can't handle this
    /// case so solve() forwards the "solve" to this function.
    void solve_for_one_dof(Problem *const &problem_pt, DoubleVector &result);

    /// Doc the solver stats or stay quiet?
    bool Doc_stats;

    /// Reorder elements with Sloan's algorithm?
    bool Reorder_flag;

    /// Use direct access files?
    bool Use_direct_access_files;

    /// \short Factor to increase storage for lenbuf[0]; see MA42 documentation
    /// for details.
    double Lenbuf_factor0;

    /// \short Factor to increase storage for lenbuf[1]; see MA42 documentation
    /// for details
    double Lenbuf_factor1;

    /// \short Factor to increase storage for lenbuf[2]; see MA42 documentation
    /// for details.
    double Lenbuf_factor2;

    /// \short Factor to increase storage for front size; see MA42 documentation
    /// for details.
    double Front_factor;

    /// \short Factor to increase size of direct access files; see
    /// MA42 documentation for details.
    double Lenfle_factor;

    /// Control flag for MA42; see MA42 documentation for details
    int Icntl[8];

    /// Control flag for MA42; see MA42 documentation for details
    int Isave[45];

    /// Control flag for MA42; see MA42 documentation for details
    int Info[23];

    /// Workspace storage for MA42
    double *W;

    /// Size of the workspace array, W
    int Lw;

    /// Integer workspace storage for MA42
    int *IW;

    /// Size of the integer workspace array
    int Liw;

    /// Size of the linear system
    unsigned long N_dof;

  public:
    /// \short Constructor: By default suppress verbose output (stats), don't
    /// reorder elements and don't use direct access files
    HSL_MA42() : W(0), Lw(0), IW(0), Liw(0), N_dof(0)
    {
      Doc_stats = false;
      Reorder_flag = false;
      Use_direct_access_files = false;

      // Default values for memory allocation
      Lenbuf_factor0 = 1.2;
      Lenbuf_factor1 = 1.2;
      Lenbuf_factor2 = 1.5;
      Front_factor = 1.1;
      Lenfle_factor = 1.5;
    }

    /// Destructor, clean up the allocated memory
    ~HSL_MA42()
    {
      clean_up_memory();
    }

    /// Broken copy constructor
    HSL_MA42(const HSL_MA42 &)
    {
      BrokenCopy::broken_copy("HSL_MA42");
    }

    /// Broken assignment operator
    void operator=(const HSL_MA42 &)
    {
      BrokenCopy::broken_assign("HSL_MA42");
    }

    /// Clean up memory
    void clean_up_memory()
    {
      if (IW)
      {
        delete[] IW;
        IW = 0;
        Liw = 0;
      }
      if (W)
      {
        delete[] W;
        W = 0;
        Lw = 0;
      }
    }

    /// Overload disable resolve so that it cleans up memory too
    void disable_resolve()
    {
      LinearSolver::disable_resolve();
      clean_up_memory();
    }

    /// \short Solver: Takes pointer to problem and returns the results Vector
    /// which contains the solution of the linear system defined by
    /// the problem's fully assembled Jacobian and residual Vector.
    void solve(Problem *const &problem_pt, DoubleVector &result);

    /// \short Linear-algebra-type solver: Takes pointer to a matrix and rhs
    /// vector and returns the solution of the linear system.
    /// Call the broken base-class version. If you want this, please implement
    /// it
    void solve(DoubleMatrixBase *const &matrix_pt,
               const DoubleVector &rhs,
               DoubleVector &result)
    {
      LinearSolver::solve(matrix_pt, rhs, result);
    }

    /// \short Linear-algebra-type solver: Takes pointer to a matrix
    /// and rhs vector and returns the solution of the linear system
    /// Call the broken base-class version. If you want this, please
    /// implement it
    void solve(DoubleMatrixBase *const &matrix_pt,
               const Vector<double> &rhs,
               Vector<double> &result)
    {
      LinearSolver::solve(matrix_pt, rhs, result);
    }

    /// \short Return the solution to the linear system Ax = result, where
    /// A is the most recently factorised jacobian matrix of the problem
    /// problem_pt. The solution is returned in the result vector.
    void resolve(const DoubleVector &rhs, DoubleVector &result);

    /// \short Function to reorder the elements based on Sloan's algorithm
    void reorder_elements(Problem *const &problem_pt);

    /// \short Enable documentation of statistics
    void enable_doc_stats()
    {
      Doc_stats = true;
    }

    /// \short Disable documentation of statistics
    void disable_doc_stats()
    {
      Doc_stats = false;
    }

    /// \short Enable reordering using Sloan's algorithm
    void enable_reordering()
    {
      Reorder_flag = true;
    }

    /// \short Disable reordering
    void disable_reordering()
    {
      Reorder_flag = false;
    }

    /// \short Enable use of direct access files
    void enable_direct_access_files()
    {
      Use_direct_access_files = true;
    }

    /// \short Disable use of direct access files
    void disable_direct_access_files()
    {
      Use_direct_access_files = false;
    }

    /// \short Factor to increase storage for lenbuf[0]; see MA42 documentation
    /// for details.
    double &lenbuf_factor0()
    {
      return Lenbuf_factor0;
    }

    /// \short Factor to increase storage for lenbuf[1]; see MA42 documentation
    /// for details.
    double &lenbuf_factor1()
    {
      return Lenbuf_factor1;
    }

    /// \short Factor to increase storage for lenbuf[2]; see MA42 documentation
    /// for details.
    double &lenbuf_factor2()
    {
      return Lenbuf_factor2;
    }

    /// \short Factor to increase storage for front size; see MA42 documentation
    /// for details.
    double &front_factor()
    {
      return Front_factor;
    }

    /// \short Factor to increase the size of the direct access files;
    /// see MA42 documentation for details.
    double &lenfle_factor()
    {
      return Lenfle_factor;
    }
  };

} // namespace oomph

#endif
