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
// Interface to HSL frontal solver (fortran)


#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

// oomph-lib headers
#include "cfortran.h"
#include "frontal.h"
#include "frontal_solver.h"
#include "problem.h"
#include "mesh.h"
#include "matrices.h"
#include "assembly_handler.h"

namespace oomph
{
  //====================================================================
  /// Special solver for problems with one DOF -- HSL_MA42 can't handle
  /// that!
  //====================================================================
  void HSL_MA42::solve_for_one_dof(Problem* const& problem_pt,
                                   DoubleVector& result)
  {
    // Find the number of elements
    unsigned long n_el = problem_pt->mesh_pt()->nelement();

#ifdef PARANOID
    // Find the number of dofs (variables)
    int n_dof = problem_pt->ndof();
    // Check
    if (n_dof != 1)
    {
      throw OomphLibError(
        "You can only call this if the problem has just one dof!",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Global jac and residuals are scalars!
    double global_jac = 0.0;
    double global_res = 0.0;

    // Locally cache pointer to assembly handler
    AssemblyHandler* const assembly_handler_pt =
      problem_pt->assembly_handler_pt();


    // Do the assembly
    for (unsigned long e = 0; e < n_el; e++)
    {
      // Get pointer to the element
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);

      // Get number of variables in element
      int nvar = assembly_handler_pt->ndof(elem_pt);

      // Only assemble if there's something to be assembled
      if (nvar > 0)
      {
        // Set up matrices and Vector for call to get_residuals
        Vector<double> residuals(nvar);
        DenseMatrix<double> jacobian(nvar);

        // Get the residuals and jacobian
        assembly_handler_pt->get_jacobian(elem_pt, residuals, jacobian);

        // Add contribution
        global_jac += jacobian(0, 0);
        global_res += residuals[0];
      }
    }

    // "Solve"
    result[0] = global_res / global_jac;

    // If we are storing the matrix, allocate the storage
    if (Enable_resolve)
    {
      // If storage has been allocated delete it
      if (W != 0)
      {
        delete[] W;
      }
      // Now allocate new storage
      W = new double[1];
      // Store the global jacobian
      W[0] = global_jac;
      // Set the number of degrees of freedom
      N_dof = 1;
    }
    // Set the sign of the jacobian
    problem_pt->sign_of_jacobian() =
      static_cast<int>(std::fabs(global_jac) / global_jac);
  }

  //====================================================================
  /// Wrapper for HSL MA42 frontal solver
  //====================================================================
  void HSL_MA42::solve(Problem* const& problem_pt, DoubleVector& result)
  {
#ifdef OOMPH_HAS_MPI
    if (problem_pt->communicator_pt()->nproc() > 1)
    {
      std::ostringstream error_message;
      error_message
        << "HSL_MA42 solver cannot be used in parallel.\n"
        << "Please change to another linear solver.\n"
        << "If you want to use a frontal solver then try MumpsSolver\n";

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // Find the number of elements
    unsigned long n_el = problem_pt->mesh_pt()->nelement();

    // Find the number of dofs (variables) and store for resolves
    N_dof = problem_pt->ndof();

    // Build the distribution .. this is a serial solver
    LinearAlgebraDistribution dist(problem_pt->communicator_pt(), N_dof, false);
    this->build_distribution(dist);

    // Cast the number of dofs into an integer for the HSL solver
    int n_dof = N_dof;

    // Special case: Just one variable! MA42 can't handle that
    if (n_dof == 1)
    {
      solve_for_one_dof(problem_pt, result);
      return;
    }

    // Reorder?
    if (this->Reorder_flag)
    {
      reorder_elements(problem_pt);
    }

    // Set the control flags
    int ifsize[5];
    double cntl[2], rinfo[2];

    // Set up memory for last and for ndf
    int ndf, nmaxe = 2, nrhs = 1, lrhs = 1;

    // Memory for last and dx should be allocated dynamically from the heap
    // rather than the stack, otherwise one gets a nasty segmentation fault
    int* last = new int[n_dof];
    // Set up memory for dx
    double** dx = new double*;
    *dx = new double[n_dof];

    // Provide storage for exact values required for lenbuf array
    int exact_lenbuf[3] = {0, 0, 0};
    bool exact_lenbuf_available = false;

    // Storage for recommendations for lenbuf so we can check how
    // close our factors are
    int lenbuf0_recommendation = 0;
    int lenbuf1_recommendation = 0;
    int lenbuf2_recommendation = 0;

    // Provide storage
    int lenbuf[3];


    // Provide storage for exact values required for lenfle array
    int exact_lenfle[3] = {0, 0, 0};
    bool exact_lenfle_available = false;

    // Storage for recommendation for lenfle so we can check how
    // close our factor is
    int lenfle0_recommendation = 0;
    int lenfle1_recommendation = 0;
    int lenfle2_recommendation = 0;

    // Provide storage
    int lenfle[3];


    // Storage for recommendations for front size so we can check how
    // close our factors are
    double front0_recommendation = 0;
    double front1_recommendation = 0;


    // Flag to indicate if exact, required values for nfront are available
    // from previous unsucessful solve
    bool exact_nfront_available = false;

    // Storage for front sizes
    int nfront[2];

    // Locally cache pointer to assembly handler
    AssemblyHandler* const assembly_handler_pt =
      problem_pt->assembly_handler_pt();


    // Loop over size allocation/solver until we've made all the arrays
    //=================================================================
    // big enough...
    //==============
    bool not_done = true;
    while (not_done)
    {
      // Initialise frontal solver
      //=========================
      MA42ID(Icntl, cntl, Isave);


      // Loop over the elements to setup variables associated with elements
      //==================================================================
      for (unsigned long e = 0; e < n_el; e++)
      {
        // Pointer to the element
        GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);

        // MH: changed array to allocateable
        int* ivar;

        // Get number of variables in element
        int nvar = assembly_handler_pt->ndof(elem_pt);


        // Multiple equations
        //-------------------
        if (nvar != 1)
        {
          // Copy global equation numbers into local array
          ivar = new int[nvar];

          for (int i = 0; i < nvar; i++)
          {
            // Need to add one to equation numbers
            ivar[i] = 1 + assembly_handler_pt->eqn_number(elem_pt, i);
          }
        }

        // Just one equation
        //------------------
        else
        {
          // Copy global equation numbers into local array - minimum length: 2
          ivar = new int[2];

          // Here's the number of the only equation
          int only_eqn = assembly_handler_pt->eqn_number(elem_pt, 0);

          // Need to add one to equation numbers
          ivar[0] = 1 + only_eqn;

          // Now add a dummy eqn either before or after the current one
          int dummy_eqn;
          if (only_eqn != (n_dof - 1))
          {
            dummy_eqn = only_eqn + 1;
          }
          else
          {
            dummy_eqn = only_eqn - 1;
          }
          ivar[1] = 1 + dummy_eqn;

          nvar = 2;
        }


        // Now call the frontal routine
        // GB: added check to exclude case with no dofs
        if (nvar)
        {
          MA42AD(nvar, ivar, ndf, last, n_dof, Icntl, Isave, Info);
          // ALH : throw an exception if the frontal_solver fails
          if (Info[0] < 0)
          {
            throw NewtonSolverError(true);
          }
        }

        // Cleanup
        delete[] ivar;
        ivar = 0;
      }


      // Loop over the elements again for symbolic factorisation
      //=======================================================
      for (unsigned long e = 0; e < n_el; e++)
      {
        // Pointer to the element
        GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);

        // MH: changed array to allocateable
        int* ivar;

        // Get number of variables in element
        int nvar = assembly_handler_pt->ndof(elem_pt);

        // Multiple equations
        //-------------------
        if (nvar != 1)
        {
          // Copy global equation numbers into local array
          ivar = new int[nvar];

          for (int i = 0; i < nvar; i++)
          {
            // Need to add one to equation numbers
            ivar[i] = 1 + assembly_handler_pt->eqn_number(elem_pt, i);
          }
        }

        // Just one equation
        //------------------
        else
        {
          // Copy global equation numbers into local array - minimum length: 2
          ivar = new int[2];

          // Here's the number of the only equation
          int only_eqn = assembly_handler_pt->eqn_number(elem_pt, 0);

          // Need to add one to equation numbers
          ivar[0] = 1 + only_eqn;

          // Now add a dummy eqn either before or after the current one
          int dummy_eqn;
          if (only_eqn != (n_dof - 1))
          {
            dummy_eqn = only_eqn + 1;
          }
          else
          {
            dummy_eqn = only_eqn - 1;
          }
          ivar[1] = 1 + dummy_eqn;

          nvar = 2;
        }

        // GB: added check to exclude case with no dofs
        if (nvar)
        {
          MA42JD(nvar, ivar, ndf, last, nmaxe, ifsize, Icntl, Isave, Info);

          // ALH : throw an exception if the frontal_solver fails
          if (Info[0] < 0)
          {
            throw NewtonSolverError(true);
          }
        }

        // Cleanup
        delete[] ivar;
        ivar = 0;

      } // end of symbolic factorisation


      // Front sizes: "Somewhat larger than..."
      //---------------------------------------
      front0_recommendation = ifsize[0];
      front1_recommendation = ifsize[1];
      if (!exact_nfront_available)
      {
        nfront[0] = int(ceil(Front_factor * double(front0_recommendation)));
        nfront[1] = int(ceil(Front_factor * double(front1_recommendation)));

        if (n_dof < nfront[0]) nfront[0] = n_dof;
        if (n_dof < nfront[1]) nfront[1] = n_dof;
      }
      // We have a problem if the front size is cocked


      // Sizes for lenbuf[]: "Somewhat larger than..."
      //----------------------------------------------
      lenbuf0_recommendation = ifsize[2] + ndf;
      // If we are storing the matrix,
      // we need to allocate storage for lenbuf_1
      if (Enable_resolve)
      {
        lenbuf1_recommendation = ifsize[3];
      }
      // Otherwise set it to zero
      else
      {
        lenbuf1_recommendation = 0;
      }
      lenbuf2_recommendation = ifsize[4];


      // Enable use of direct access files?
      //----------------------------------
      if (Use_direct_access_files)
      {
        if (Doc_stats)
        {
          oomph_info << "Using direct access files" << std::endl;
        }

        // Unit numbers for direct access files (middle one only used for
        // re-solve; set to zero as this is not used here).
        int istrm[3];
        istrm[0] = 17;
        // If we are storing the matrix, set the stream
        if (Enable_resolve)
        {
          istrm[1] = 19;
        }
        // otherwise, set it to zero
        else
        {
          istrm[1] = 0;
        }
        istrm[2] = 18;

        // Set up the memory: Need memory "somewhat larger" than ...
        double factor = 1.1;
        lenbuf[0] = int(ceil(factor * double(10 * (nfront[0] + 1))));
        // If we are storing the matrix, set the value
        if (Enable_resolve)
        {
          lenbuf[1] = int(ceil(factor * double(10 * (nfront[0]))));
        }
        // Otherwise, set to zero
        else
        {
          lenbuf[1] = 0;
        }
        lenbuf[2] = int(ceil(factor * double(10 * (2 * nfront[0] + 5))));


        // Initial size allocation: Need memory "somewhat larger" than ...
        if (exact_lenfle_available)
        {
          lenfle[0] = exact_lenfle[0];
          lenfle[1] = exact_lenfle[1];
          lenfle[2] = exact_lenfle[2];
        }
        else
        {
          lenfle0_recommendation = ifsize[2] + ndf;
          lenfle1_recommendation = ifsize[3];
          lenfle2_recommendation = ifsize[4];
          lenfle[0] = int(ceil(Lenfle_factor * double(lenfle0_recommendation)));
          lenfle[1] = int(ceil(Lenfle_factor * double(lenfle1_recommendation)));
          lenfle[2] = int(ceil(Lenfle_factor * double(lenfle2_recommendation)));
        }

        // Setup direct access files
        MA42PD(istrm, lenbuf, lenfle, Icntl, Isave, Info);
      }
      else
      {
        if (Doc_stats)
        {
          oomph_info << "Not using direct access files" << std::endl;
        }

        // Initial size allocation: Need memory "somewhat larger" than ...
        if (exact_lenbuf_available)
        {
          lenbuf[0] = exact_lenbuf[0];
          lenbuf[1] = exact_lenbuf[1];
          lenbuf[2] = exact_lenbuf[2];
        }
        else
        {
          lenbuf[0] =
            int(ceil(Lenbuf_factor0 * double(lenbuf0_recommendation)));
          // If we are storing the matrix, set the buffer size
          if (Enable_resolve)
          {
            lenbuf[1] =
              int(ceil(Lenbuf_factor1 * double(lenbuf1_recommendation)));
          }
          // Otherwise its zero
          else
          {
            lenbuf[1] = 0;
          }
          lenbuf[2] =
            int(ceil(Lenbuf_factor2 * double(lenbuf2_recommendation)));
        }
      }


      if (Doc_stats)
      {
        oomph_info << "\n FRONT SIZE (min and actual): " << ifsize[0] << " "
                   << nfront[0] << std::endl;
      }


      // Initialise success
      bool success = true;

      // Workspace arrays: calculate amount in local variables initiall
      int lw = 1 + lenbuf[0] + lenbuf[1] + nfront[0] * nfront[1];
      if (lrhs * nfront[0] > nrhs * nfront[1])
      {
        lw += lrhs * nfront[0];
      }
      else
      {
        lw += nrhs * nfront[1];
      }

      int liw = lenbuf[2] + 2 * nfront[0] + 4 * nfront[1];

      // Set up memory for w and iw
      // Again allocate dynamically from the heap, rather than the stack
      // If the required amount of storage has changed, reallocate
      // and store the values in the object member data
      if (lw != Lw)
      {
        // Set the new value of Lw (member data)
        Lw = lw;
        // Delete the previously allocated storage
        if (W) delete[] W;
        // Now allocate new storage
        W = new (std::nothrow) double[Lw];
        if (!W)
        {
          throw OomphLibError("Out of memory in allocation of w\n",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // If the required amount of storage has changed, reallocate
      if (liw != Liw)
      {
        // Set the new value of Liw (member data)
        Liw = liw;
        // Delete the previously allocated storgae
        if (IW != 0)
        {
          delete[] IW;
        }
        // Now allocate new storage
        IW = new (std::nothrow) int[Liw];
        if (!IW)
        {
          throw OomphLibError("Out of memory in allocation of iw\n",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Doc
      if (Doc_stats)
      {
        double temp = (lenbuf[0] * sizeof(double) + lenbuf[2] * sizeof(int)) /
                      (1024.0 * 1024.0);
        oomph_info << "\n ESTIMATED MEMORY USAGE " << temp << "Mb" << std::endl;
      }


      // Now do the actual frontal assembly and solve loop
      //=================================================
      for (unsigned long e = 0; e < n_el; e++)
      {
        // Get pointer to the element
        GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);

        // Get number of variables in element
        int nvar = assembly_handler_pt->ndof(elem_pt);

        // MH: changed array to allocatable
        int* ivar;

        // Multiple equations
        //-------------------
        if (nvar != 1)
        {
          // Copy global equation numbers into local array
          ivar = new int[nvar];

          for (int i = 0; i < nvar; i++)
          {
            // Need to add one to equation numbers
            ivar[i] = 1 + assembly_handler_pt->eqn_number(elem_pt, i);
          }
          nmaxe = nvar;
        }

        // Just one equation
        //------------------
        else
        {
          // Copy global equation numbers into local array - minimum length: 2
          ivar = new int[2];

          // Here's the number of the only equation
          int only_eqn = assembly_handler_pt->eqn_number(elem_pt, 0);

          // Need to add one to equation numbers
          ivar[0] = 1 + only_eqn;

          // Now add a dummy eqn either before or after the current one
          int dummy_eqn;
          if (only_eqn != (n_dof - 1))
          {
            dummy_eqn = only_eqn + 1;
          }
          else
          {
            dummy_eqn = only_eqn - 1;
          }
          ivar[1] = 1 + dummy_eqn;

          nmaxe = 2;
        }

        // Set up matrices and Vector for call to get_residuals
        Vector<double> residuals(nvar);
        DenseMatrix<double> jacobian(nvar);

        // Get the residuals and jacobian
        assembly_handler_pt->get_jacobian(elem_pt, residuals, jacobian);

        // Add dummy rows and columns if required
        if (nvar == 1)
        {
          double onlyjac = jacobian(0, 0);
          jacobian.resize(2, 2);
          jacobian(0, 0) = onlyjac;
          jacobian(1, 0) = 0.0;
          jacobian(0, 1) = 0.0;
          jacobian(1, 1) = 0.0;

          double onlyres = residuals[0];
          residuals.resize(2);
          residuals[0] = onlyres;
          residuals[1] = 0.0;

          nvar = 2;
        }


        // GB: added check to exclude case with no dofs
        if (nvar)
        {
          // Now translate the residuals and jacobian into form to be
          // taken by stupid FORTRAN frontal solver
          // double avar[nvar][nmaxe], rhs[1][nmaxe];
          // Assign memory on the heap rather than the stack because the ndofs
          // could get too large
          double** avar = new double*[nvar];
          double* tmp = new double[nvar * nmaxe];
          for (int i = 0; i < nvar; i++)
          {
            avar[i] = &tmp[i * nmaxe];
          }
          double** rhs = new double*[1];
          rhs[0] = new double[nmaxe];


          for (int i = 0; i < nmaxe; i++)
          {
            rhs[0][i] = residuals[i];
            for (int j = 0; j < nvar; j++)
            {
              // Note need to transpose here
              avar[j][i] = jacobian(i, j);
            }
          }

          // Call the frontal solver
          MA42BD(nvar,
                 ivar,
                 ndf,
                 last,
                 nmaxe,
                 avar,
                 nrhs,
                 rhs,
                 lrhs,
                 n_dof,
                 dx,
                 nfront,
                 lenbuf,
                 Lw,
                 W,
                 Liw,
                 IW,
                 Icntl,
                 cntl,
                 Isave,
                 Info,
                 rinfo);


          // Error check:
          if (Info[0] < 0)
          {
            if (Doc_stats)
              oomph_info << "Error: Info[0]=" << Info[0] << std::endl;

            // Solve isn't going to be successful
            success = false;

            // Error: Entries in lenbuf too small -- this is covered
            if (Info[0] == -16)
            {
              if (Doc_stats)
                oomph_info << "lenbuf[] too small -- can recover..."
                           << std::endl;
            }
            else if (Info[0] == -12)
            {
              if (Doc_stats)
                oomph_info << "nfront[] too small -- can recover..."
                           << std::endl;
            }
            else if (Info[0] == -17)
            {
              if (Doc_stats)
                oomph_info << "lenfle[] too small -- can recover..."
                           << std::endl;
            }
            else
            {
              if (Doc_stats)
              {
                oomph_info << "Can't recover from this error" << std::endl;
              }
              throw NewtonSolverError(true);
            }
          }

          // Clean up the memory
          delete[] rhs[0];
          rhs[0] = 0;
          delete[] rhs;
          rhs = 0;
          delete[] avar[0];
          avar[0] = 0;
          delete[] avar;
        }

        // Cleanup
        delete[] ivar;
        ivar = 0;

      } // end of loop over elements for assemble/solve

      // Set the sign of the jacobian matrix
      problem_pt->sign_of_jacobian() = Info[1];

      // If we are not storing the matrix, then delete the assigned workspace
      if (!Enable_resolve)
      {
        // Free the workspace memory assigned and set stored values to zero
        delete[] IW;
        IW = 0;
        Liw = 0;
        delete[] W;
        W = 0;
        Lw = 0;
        // Reset the number of dofs
        N_dof = 0;
      }

      // We've done the elements -- did we encounter an error
      // that forces us to re-allocate workspace?
      if (success)
      {
        // We're done!
        not_done = false;
      }
      else
      {
        // Reset sizes for lenbuf
        if (!Use_direct_access_files)
        {
          exact_lenbuf[0] = Info[4];
          exact_lenbuf[1] = Info[5];
          exact_lenbuf[2] = Info[6];
          exact_lenbuf_available = true;
        }

        // nfront[] has already been overwritten with required values -- don't
        // need to update anything. Just indicate that new values are
        // available so they don't have to be re-assigned.
        exact_nfront_available = true;

        // Reset sizes for lenfle
        exact_lenfle[0] = lenbuf[0] * Info[9];
        exact_lenfle[1] = lenbuf[1] * Info[10];
        exact_lenfle[2] = lenbuf[2] * Info[11];
        exact_lenfle_available = true;
      }

    } // end of loop over shuffling of workspace array sizes until it all
      // fits...

    result.build(&dist);
    // Now load the values back into result
    for (int i = 0; i < n_dof; i++) result[i] = dx[0][i];

    // Print the actual memory used
    if (Doc_stats)
    {
      if (!Use_direct_access_files)
      {
        double temp = (Info[4] * sizeof(double) + Info[6] * sizeof(int)) /
                      (1024.0 * 1024.0);
        oomph_info << " TOTAL MEMORY USED " << temp << "Mb" << std::endl;
      }
      oomph_info << std::endl;
      oomph_info << "lenbuf[]: " << lenbuf[0] << " " << lenbuf[1] << " "
                 << lenbuf[2] << " " << std::endl;
      oomph_info << "lenbuf[] factors required and initially specified:"
                 << std::endl;
      oomph_info << "lenbuf[0]: " << Info[4] / double(lenbuf0_recommendation)
                 << " " << Lenbuf_factor0 << std::endl;
      oomph_info << "lenbuf[1]: " << Info[5] / double(lenbuf1_recommendation)
                 << " " << Lenbuf_factor1 << std::endl;

      oomph_info << "lenbuf[2]: " << Info[6] / double(lenbuf2_recommendation)
                 << " " << Lenbuf_factor2 << std::endl;
      oomph_info << "nfront[] factors required and initially specified:"
                 << std::endl;
      oomph_info << "nfront[0]: " << nfront[0] / double(front0_recommendation)
                 << " " << Front_factor << std::endl;
      oomph_info << "nfront[1]: " << nfront[1] / double(front1_recommendation)
                 << " " << Front_factor << std::endl;
      if (Use_direct_access_files)
      {
        oomph_info << "lenfle[]: " << lenfle[0] << " " << lenfle[1] << " "
                   << lenfle[2] << " " << std::endl;
        oomph_info << "lenfle[] factors required and initially specified:"
                   << std::endl;
        oomph_info << "lenfle[0]: "
                   << Info[9] * lenbuf[0] / double(lenfle0_recommendation)
                   << " " << Lenfle_factor << std::endl;
        oomph_info << "lenfle[1]: "
                   << Info[10] * lenbuf[1] / double(lenfle1_recommendation)
                   << " " << Lenfle_factor << std::endl;
        oomph_info << "lenfle[2]: "
                   << Info[11] * lenbuf[2] / double(lenfle2_recommendation)
                   << " " << Lenfle_factor << std::endl;
      }
    }

    // Free the memory assigned
    delete[] * dx;
    delete dx;
    delete[] last;
  }

  //====================================================================
  /// Wrapper for HSL MA42 frontal solver
  //====================================================================
  void HSL_MA42::resolve(const DoubleVector& rhs, DoubleVector& result)
  {
    // If we haven't stored the factors complain
    if (W == 0)
    {
      throw OomphLibError("Resolve called without first solving",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Find the number of dofs (variables)
    int n_dof = N_dof;

    // Check that the number of DOFs is equal to the size of the RHS vector
#ifdef PARANOID
    if (n_dof != static_cast<int>(rhs.nrow()))
    {
      throw OomphLibError(
        "RHS does not have the same dimension as the linear system",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Special case: just one variable! MA42 can't handle it
    // Solve using the stored jacobian (single double)
    if (n_dof == 1)
    {
      result[0] = rhs[0] / W[0];
      return;
    }

    // Otherwise actually call the MA42 routine
    // We are solving the original matrix (not the transpose)
    int trans = 0;
    // There is only one rhs
    int nrhs = 1;
    // The size of the vectors is the number of degrees of freedom in the
    // problem
    int lx = n_dof;

    // Allocate storage for the RHS
    double** b = new double*;
    *b = new double[n_dof];
    // Load in the RHS vector
    for (int i = 0; i < n_dof; i++)
    {
      b[0][i] = rhs[i];
    }

    // Storage for the results
    double** x = new double*;
    *x = new double[n_dof];

    // Now call the resolver
    MA42CD(trans, nrhs, lx, b, x, Lw, W, Liw, IW, Icntl, Isave, Info);

    // If there has been an error throw it
    if (Info[0] < 0)
    {
      throw NewtonSolverError(true);
    }

    result.build(this->distribution_pt());
    // Put the answer into the result array
    for (int i = 0; i < n_dof; ++i)
    {
      result[i] = x[0][i];
    }

    // Delete the allocated storage
    delete[] * x;
    delete x;
    delete[] * b;
    delete b;
  }


  //=========================================================================
  /// Function to reorder the elements according to Sloan's algorithm
  //=========================================================================
  void HSL_MA42::reorder_elements(Problem* const& problem_pt)
  {
    // Find the number of elements
    unsigned long n_el = problem_pt->mesh_pt()->nelement();

    // Find the number of dofs
    int n_dof = problem_pt->ndof();

    Vector<int> order(n_el);

    // Control parameters
    int icntl[10], info[15], direct, nsup, vars[n_dof], svar[n_dof];
    int liw, lw, *perm = 0;
    double wt[3], rinfo[6];


    // Direct ordering: 1
    direct = 1;

    // Initial guess for workspaces (deliberately too small so the
    // routine fails and returns the required sizes
    liw = 1;
    lw = 1;

    // Initialise things here
    MC63ID(icntl);


    // Suppress printing of error and warning messages?
    if (!Doc_stats)
    {
      icntl[0] = -1;
      icntl[1] = -1;
    }

    // Set up iw and w
    int* iw = new int[liw];
    double* w = new double[lw];

    // Set up the vectors that hold the element connectivity information
    // Can initialise the first entry of eltptr to 1
    Vector<int> eltvar, eltptr(1, 1);

    // Locally cache pointer to assembly handler
    AssemblyHandler* const assembly_handler_pt =
      problem_pt->assembly_handler_pt();

    // Now loop over all the elements and assemble eltvar and eltptr
    for (unsigned long e = 0; e < n_el; e++)
    {
      // Set up the default order
      order[e] = e;

      // Get pointer to the element
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);

      // Get the number of variables in the element via the assembly handler
      int nvar = assembly_handler_pt->ndof(elem_pt);

      // Multiple equations
      //-------------------
      if (nvar != 1)
      {
        // Copy the equation numbers into the global array
        for (int i = 0; i < nvar; i++)
        {
          // Need to add one to equation numbers
          eltvar.push_back(1 + assembly_handler_pt->eqn_number(elem_pt, i));
        }
        eltptr.push_back(eltptr.back() + nvar);
      }
      // Just one equation
      //------------------
      else
      {
        // Here's the number of the only equation
        int only_eqn = assembly_handler_pt->eqn_number(elem_pt, 0);

        // Need to add one to equation numbers
        eltvar.push_back(1 + only_eqn);

        // Now add a dummy eqn either before or after the current one
        int dummy_eqn;
        if (only_eqn != (n_dof - 1))
        {
          dummy_eqn = only_eqn + 1;
        }
        else
        {
          dummy_eqn = only_eqn - 1;
        }

        eltvar.push_back(1 + dummy_eqn);

        eltptr.push_back(eltptr.back() + 2);
      }

    } // End of loop over the elements

    // Set the value of n_e (number of entries in the element variable list
    int n_e = eltvar.size();

    do
    {
      // Call the reordering routine
      MC63AD(direct,
             n_dof,
             n_el,
             n_e,
             &eltvar[0],
             &eltptr[0],
             &order[0],
             perm,
             nsup,
             vars,
             svar,
             wt,
             liw,
             iw,
             lw,
             w,
             icntl,
             info,
             rinfo);

      // If not enough space in iw or w, allocate enough and try again
      if (info[0] == -4)
      {
        if (Doc_stats)
          oomph_info << " Reallocating liw to " << info[4] << std::endl;
        delete[] iw;
        liw = info[4];
        iw = new int[liw];
      }

      // If not enough space in w, allocate enough and try again
      if (info[0] == -2)
      {
        if (Doc_stats)
          oomph_info << " Reallocating lw to " << info[5] << std::endl;
        delete[] w;
        lw = info[5];
        w = new double[lw];
      }

    } while (info[0] < 0);


    if (Doc_stats)
    {
      oomph_info << "\n Initial wavefront details  :\n      max " << rinfo[0]
                 << "\tmean " << rinfo[1] << "\tprofile " << rinfo[2];
      oomph_info << "\n Reordered wavefront details:\n      max " << rinfo[3]
                 << "\tmean " << rinfo[4] << "\tprofile " << rinfo[5];
      oomph_info << std::endl;
    }

    // Free the memory allocated
    delete[] w;
    delete[] iw;


    // Store re-ordered elements in temporary vector
    Vector<GeneralisedElement*> tmp_element_pt(n_el);
    for (unsigned e = 0; e < n_el; e++)
    {
      tmp_element_pt[e] = problem_pt->mesh_pt()->element_pt(order[e] - 1);
    }
    for (unsigned e = 0; e < n_el; e++)
    {
      problem_pt->mesh_pt()->element_pt(e) = tmp_element_pt[e];
    }
  }

} // namespace oomph
