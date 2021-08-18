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

#include <algorithm>

#include "segregated_fsi_solver.h"
#include <algorithm>

namespace oomph
{
  //============================================================================
  /// Extrapolate solid data and update fluid mesh. Extrapolation is based
  /// on previous history values and is therefore only called in
  /// time-dependent runs.
  //============================================================================
  void SegregatableFSIProblem::extrapolate_solid_data()
  {
    // Number of solid Data items:
    unsigned n_data = Solid_data_pt.size();

    // Loop over all solid Data items:
    for (unsigned i = 0; i < n_data; i++)
    {
      // Number of values stored at this Data item:
      unsigned n_value = Solid_data_pt[i]->nvalue();

      // Number of values stored at this Data item:
      unsigned n_prev = Solid_data_pt[i]->time_stepper_pt()->nprev_values();

      // Can't use this extrapolation if we don't have at least two
      // history values; may be able to do better for higher order
      // timesteppers, though.
      if (n_prev >= 2)
      {
        // Extrapolate all values
        for (unsigned k = 0; k < n_value; k++)
        {
          // Linear extrapolation based on previous two values,
          // assuming constant timestep.
          double new_value =
            2.0 * Solid_data_pt[i]->value(1, k) - Solid_data_pt[i]->value(2, k);
          *(Solid_data_pt[i]->value_pt(0, k)) = new_value;
        }
      }
    }

    // Update fluid node position in response to changes in wall shape
    // (If fluid mesh was specified via non-NULL Fluid_mesh_pt
    // this node update will only act on the fluid nodes).
    if (Fluid_mesh_pt != 0)
    {
      Fluid_mesh_pt->node_update();
    }
    else
    {
      mesh_pt()->node_update();
    }
  }


  //============================================================================
  /// Rebuild global mesh for monolithic problem
  //============================================================================
  void SegregatableFSIProblem::rebuild_monolithic_mesh()
  {
    // Get rid of the previous submeshes
    flush_sub_meshes();

    // Add original submeshes
    unsigned orig_n_sub_mesh = Orig_sub_mesh_pt.size();
    for (unsigned i = 0; i < orig_n_sub_mesh; i++)
    {
      add_sub_mesh(Orig_sub_mesh_pt[i]);
    }

    // Rebuild global mesh
    rebuild_global_mesh();
  }


  //============================================================================
  /// Only include solid elements in the Problem's mesh. This is
  /// called before the segregated solid solve. The solid elements are
  /// identified by the user via the solid_mesh_pt argument
  /// in the pure virtual function identify_fluid_and_solid_dofs(...).
  /// If a NULL pointer is returned by this function (i.e. if the user
  /// hasn't bothered to identify the solids elements in a submesh, then
  /// no stripping of non-solid elements is performed. The result
  /// of the computation will be correct but
  /// it won't be very efficient.
  //============================================================================
  void SegregatableFSIProblem::use_only_solid_elements()
  {
    // Wipe global mesh and rebuild it with just the solid elements in it
    if (Solid_mesh_pt != 0)
    {
      flush_sub_meshes();
      add_sub_mesh(Solid_mesh_pt);
      rebuild_global_mesh();
    }
  }


  //============================================================================
  /// Only include fluid elements in the Problem's mesh. This is
  /// called before the segregated fluid solve. The fluid elements are
  /// identified by the user via the fluid_mesh_pt argument
  /// in the pure virtual function identify_fluid_and_solid_dofs(...).
  /// If a NULL pointer is returned by this function (i.e. if the user
  /// hasn't bothered to identify the fluids elements in a submesh, then
  /// no stripping of non-fluid elements is performed. The result
  /// of the computation will be correct but
  /// it won't be very efficient.
  //============================================================================
  void SegregatableFSIProblem::use_only_fluid_elements()
  {
    // Wipe global mesh and rebuild it with just the fluid elements in it
    if (Fluid_mesh_pt != 0)
    {
      flush_sub_meshes();
      add_sub_mesh(Fluid_mesh_pt);
      rebuild_global_mesh();
    }
  }


  //============================================================================
  /// Pin all fluid dofs
  //============================================================================
  void SegregatableFSIProblem::pin_fluid_dofs()
  {
    // Number of fluid Data items:
    unsigned n_data = Fluid_data_pt.size();

    // Loop over all fluid Data items:
    for (unsigned i = 0; i < n_data; i++)
    {
      // Number of values stored at this Data item:
      unsigned n_value = Fluid_data_pt[i]->nvalue();

      // Pin all values
      for (unsigned k = 0; k < n_value; k++)
      {
        Fluid_data_pt[i]->pin(k);
      }
    }
  }


  //============================================================================
  /// Restore the pinned status of all fluid dofs
  //============================================================================
  void SegregatableFSIProblem::restore_fluid_dofs()
  {
    // Number of fluid Data items:
    unsigned n_data = Fluid_data_pt.size();

    // Loop over all fluid Data items:
    for (unsigned i = 0; i < n_data; i++)
    {
      // Number of values stored at this Data item:
      unsigned n_value = Fluid_data_pt[i]->nvalue();

      // Pin all values
      for (unsigned k = 0; k < n_value; k++)
      {
        if (Fluid_value_is_pinned[i][k])
        {
          Fluid_data_pt[i]->pin(k);
        }
        else
        {
          Fluid_data_pt[i]->unpin(k);
        }
      }
    }
  }


  //============================================================================
  /// Pin all solid dofs
  //============================================================================
  void SegregatableFSIProblem::pin_solid_dofs()
  {
    // Number of solid Data items:
    unsigned n_data = Solid_data_pt.size();

    // Loop over all solid Data items:
    for (unsigned i = 0; i < n_data; i++)
    {
      // Number of values stored at this Data item:
      unsigned n_value = Solid_data_pt[i]->nvalue();

      // Pin all values
      for (unsigned k = 0; k < n_value; k++)
      {
        Solid_data_pt[i]->pin(k);
      }
    }
  }


  //============================================================================
  /// Restore the pinned status of all solid dofs
  //============================================================================
  void SegregatableFSIProblem::restore_solid_dofs()
  {
    // Number of solid Data items:
    unsigned n_data = Solid_data_pt.size();

    // Loop over all solid Data items:
    for (unsigned i = 0; i < n_data; i++)
    {
      // Number of values stored at this Data item:
      unsigned n_value = Solid_data_pt[i]->nvalue();

      // Pin all values
      for (unsigned k = 0; k < n_value; k++)
      {
        if (Solid_value_is_pinned[i][k])
        {
          Solid_data_pt[i]->pin(k);
        }
        else
        {
          Solid_data_pt[i]->unpin(k);
        }
      }
    }
  }


  //============================================================================
  /// Store the current solid values as reference values for future
  /// convergence check. Also add another entry to pointwise Aitken history.
  //============================================================================
  void SegregatableFSIProblem::store_solid_dofs()
  {
    // Counter for the number of values:
    unsigned value_count = 0;

    // Number of solid Data items:
    unsigned n_data = Solid_data_pt.size();

    // Loop over all solid Data items:
    for (unsigned i = 0; i < n_data; i++)
    {
      // Number of values stored at this Data item:
      unsigned n_value = Solid_data_pt[i]->nvalue();

      // Loop over all values
      for (unsigned k = 0; k < n_value; k++)
      {
        // Make backup
        Previous_solid_value[value_count] = Solid_data_pt[i]->value(k);

        // Store in pointwise Aitken history
        if (Use_pointwise_aitken && (Pointwise_aitken_counter >= 0))
        {
          Pointwise_aitken_solid_value[value_count][Pointwise_aitken_counter] =
            Solid_data_pt[i]->value(k);
        }

        // Increment counter
        value_count++;
      }
    }

    // We stored another level of Aitken history values
    Pointwise_aitken_counter++;
  }


  //============================================================================
  /// Get rms of change in the solid dofs; the max. change of the
  /// solid dofs and the rms norm of the solid dofs themselves.
  /// Change is computed relative to the reference values stored when
  /// store_solid_dofs() was last called.
  //============================================================================
  void SegregatableFSIProblem::get_solid_change(double& rms_change,
                                                double& max_change,
                                                double& rms_norm)
  {
    // Initialise
    rms_change = 0.0;
    max_change = 0.0;
    rms_norm = 0.0;

    // Counter for the number of values:
    unsigned value_count = 0;

    // Number of solid Data items:
    unsigned n_data = Solid_data_pt.size();

    // Loop over all solid Data items:
    for (unsigned i = 0; i < n_data; i++)
    {
      // Number of values stored at this Data item:
      unsigned n_value = Solid_data_pt[i]->nvalue();

      // Loop over all values
      for (unsigned k = 0; k < n_value; k++)
      {
        // Change
        double change =
          Previous_solid_value[value_count] - Solid_data_pt[i]->value(k);

        // Max change?
        if (std::fabs(change) > max_change) max_change = std::fabs(change);

        // Add square of change relative to previous value
        rms_change += pow(change, 2);

        // Add square of previous value
        rms_norm += pow(Previous_solid_value[value_count], 2);

        // Increment counter
        value_count++;
      }
    }

    // Turn into rms:
    rms_change = sqrt(rms_change / double(value_count));
    rms_norm = sqrt(rms_norm / double(value_count));
  }


  //============================================================================
  /// Under-relax solid, either by classical relaxation or by Irons & Tuck.
  //============================================================================
  void SegregatableFSIProblem::under_relax_solid()
  {
    // Irons and Tuck extrapolation/relaxation; an extension of Aitken's method
    //-------------------------------------------------------------------------
    if (Use_irons_and_tuck_extrapolation)
    {
      double top = 0.0;
      double den = 0.0;
      double crit = 0.0;

      // Counter for the number of values:
      unsigned value_count = 0;

      // Number of solid Data items:
      unsigned n_data = Solid_data_pt.size();

      // Loop over all solid Data items:
      for (unsigned i = 0; i < n_data; i++)
      {
        // Number of values stored at this Data item:
        unsigned n_value = Solid_data_pt[i]->nvalue();

        // Loop over all values
        for (unsigned k = 0; k < n_value; k++)
        {
          // Prediction from solid solver
          double new_solid_value = Solid_data_pt[i]->value(k);

          // Previus value
          double old_solid_value = Previous_solid_value[value_count];

          // Change
          double change = old_solid_value - new_solid_value;

          // Change of change
          double del2 = Del_irons_and_tuck[value_count] - change;

          // Update change
          Del_irons_and_tuck[value_count] = change;

          // Update top
          top += del2 * change;

          // Update denominator
          den += del2 * del2;

          // Update convergence criterion
          crit += std::fabs(change);

          // Increment counter
          value_count++;
        }
      }

      // Update relaxation factor. The if buffers the case in which
      // we haven't realised that we've converged (so that den=0).
      // This can happen, e.g. if the convergence assessment is based on the
      // global residual or during validation. In that case we
      // obviously don't want any changes to the iterates.
      if (den != 0.0)
      {
        double new_r = R_irons_and_tuck + (R_irons_and_tuck - 1.0) * top / den;
        R_irons_and_tuck = new_r;
      }
      else
      {
        R_irons_and_tuck = 0.0;
      }

      // Loop over all solid Data items for update
      value_count = 0;
      for (unsigned i = 0; i < n_data; i++)
      {
        // Number of values stored at this Data item:
        unsigned n_value = Solid_data_pt[i]->nvalue();

        // Loop over all values
        for (unsigned k = 0; k < n_value; k++)
        {
          // Compute relaxed/extrapolated value
          double new_value = Solid_data_pt[i]->value(k) +
                             R_irons_and_tuck * Del_irons_and_tuck[value_count];

          // Assign
          Solid_data_pt[i]->set_value(k, new_value);

          // Increment counter
          value_count++;
        }
      }
      return;
    }

    // Standard relaxation
    //--------------------
    else
    {
      // No relaxation: Can return immediately
      if (Omega_relax == 1.0) return;

      // Counter for the number of values:
      unsigned value_count = 0;

      // Number of solid Data items:
      unsigned n_data = Solid_data_pt.size();

      // Loop over all solid Data items:
      for (unsigned i = 0; i < n_data; i++)
      {
        // Number of values stored at this Data item:
        unsigned n_value = Solid_data_pt[i]->nvalue();

        // Loop over all values
        for (unsigned k = 0; k < n_value; k++)
        {
          // Prediction from solid solver
          double new_solid_value = Solid_data_pt[i]->value(k);

          // Previus value
          double old_solid_value = Previous_solid_value[value_count];

          // Relax
          Solid_data_pt[i]->set_value(k,
                                      Omega_relax * new_solid_value +
                                        (1.0 - Omega_relax) * old_solid_value);

          // Increment counter
          value_count++;
        }
      }
    }
  }

  //============================================================================
  /// Pointwise Aitken extrapolation for solid variables
  //============================================================================
  void SegregatableFSIProblem::pointwise_aitken_extrapolate()
  {
    // Counter for the number of values:
    unsigned value_count = 0;

    // Number of solid Data items:
    unsigned n_data = Solid_data_pt.size();

    // Loop over all solid Data items:
    for (unsigned i = 0; i < n_data; i++)
    {
      // Number of values stored at this Data item:
      unsigned n_value = Solid_data_pt[i]->nvalue();

      // Loop over all values
      for (unsigned k = 0; k < n_value; k++)
      {
        // Shorthand
        double v0 = Pointwise_aitken_solid_value[value_count][0];
        double v1 = Pointwise_aitken_solid_value[value_count][1];
        double v2 = Pointwise_aitken_solid_value[value_count][2];

        double new_value = v2;

        double max_diff = std::max(std::fabs(v1 - v0), std::fabs(v2 - v1));
        if (max_diff > 1.0e-10)
        {
          new_value = v2 - std::pow((v2 - v1), int(2)) / (v2 - 2.0 * v1 + v0);
        }
        Solid_data_pt[i]->set_value(k, new_value);

        // Increment counter
        value_count++;
      }
    }

    // Reset the counter for the Aitken convergence check
    // (setting counter to -1 forces three new genuine
    // iterates to be computed).
    Pointwise_aitken_counter = -1;
  }


  //============================================================================
  /// Segregated fixed point solver with optional pointwise Aitken acceleration
  /// on selected solid dofs. Returns PicardConvergenceData object
  /// that contains the vital stats of the iteration.
  //============================================================================
  PicardConvergenceData SegregatableFSIProblem::segregated_solve()
  {
    // Initialise timer for essential bits of code
    reset_timer();
    // Start timer for overall solve
    clock_t t_total_start = clock();

    // If convergence is based on max. global residual we may as well
    // document it...
    bool get_max_global_residual = Doc_max_global_residual;
    if (Convergence_criterion ==
        Assess_convergence_based_on_max_global_residual)
    {
      get_max_global_residual = true;
    }

    // Create object to doc convergence stats
    PicardConvergenceData conv_data;

    // Initialise total time for computation of global residuals
    double cpu_for_global_residual = 0.0;

    // Update anything that needs updating
    actions_before_segregated_solve();

    // Set flags to values that are appropriate if Picard iteration
    // does not converge with Max_picard iterations
    bool converged = false;
    unsigned iter_taken = Max_picard;

    // This one will be overwritten during the convergence checks
    double tol_achieved = 0.0;

    // Store the current values of the solid dofs as reference values
    // and for the pointwise Aitken extrapolation
    store_solid_dofs();

    // Loop over Picard iterations
    //----------------------------
    for (unsigned iter = 1; iter <= Max_picard; iter++)
    {
      // Calculate the initial residuals
      if (iter == 1)
      {
        // Problem is always non-linear?
        // Perform any actions before the convergence check
        actions_before_segregated_convergence_check();

        double max_res = 0.0;
        if (get_max_global_residual)
        {
          clock_t t_start = clock();
          restore_solid_dofs();
          restore_fluid_dofs();
          rebuild_monolithic_mesh();
          assign_eqn_numbers();
          // This is now a full problem
          Solve_type = Full_solve;
          DoubleVector residual;
          get_residuals(residual);
          // Get maximum residuals using our own abscmp function
          max_res = residual.max();
          clock_t t_end = clock();
          cpu_for_global_residual += double(t_end - t_start) / CLOCKS_PER_SEC;
        }

        oomph_info << "==================================================\n";
        oomph_info << "Initial iteration             : " << 0 << std::endl;
        oomph_info << "RMS  change           :       " << 0 << std::endl;
        oomph_info << "Max. change           :       " << 0 << std::endl;
        oomph_info << "RMS norm              :       " << 0 << std::endl;
        if (get_max_global_residual)
        {
          oomph_info << "Max. global residual  :       " << max_res
                     << std::endl;
        }
        oomph_info << "==================================================\n\n";

        // Check for convergence, but this only makes sense
        // for the maximum (rather than relative case)
        if ((Convergence_criterion ==
             Assess_convergence_based_on_max_global_residual) &&
            (max_res < Convergence_tolerance))
        {
          oomph_info
            << "\n\n\n////////////////////////////////////////////////////////"
            << "\nPicard iteration converged after " << 0 << " steps!"
            << std::endl
            << "Convergence was based on max. residual of coupled eqns \n"
            << "being less than " << Convergence_tolerance << "." << std::endl
            << "////////////////////////////////////////////////////////\n\n\n"
            << std::endl;

          // Converged, so bail out
          // Number of iterations taken
          iter_taken = 0;

          // We have converged (overwrites default of false)
          conv_data.set_solver_converged();

          // Break loop using a GO TO! This is the first (and hopefully
          // the last) one in the whole of oomph-lib. Here it's
          // it's actually the cleanest way to exit from these
          // nested control structures
          goto jump_out_of_picard;
        }
      }

      // Stage 1: Freeze wall and solve single-physics fluid problem
      //------------------------------------------------------------
      // Pin the solid dofs, restore the pinned status of the fluid dofs
      // and re-assign the equation numbers
      restore_fluid_dofs();
      pin_solid_dofs();
      use_only_fluid_elements();
      assign_eqn_numbers();

      // Solve the fluid problem for the current wall geometry
      oomph_info << "\n\nDOING FLUID SOLVE\n\n";
      // This is a fluid solve at the moment done by a newton solve
      Solve_type = Fluid_solve;
      newton_solve();

      // Stage 2: Freeze fluid variables and update wall shape
      //------------------------------------------------------

      // Now restore the pinned status of the wall and pin the fluid
      // dofs in anticipation of the wall update:
      restore_solid_dofs();
      pin_fluid_dofs();
      use_only_solid_elements();
      assign_eqn_numbers();

      // Solve the solid problem for the wall solution
      oomph_info << "\n\nDOING SOLID SOLVE\n\n";
      // This is a solid solve, at the moment done by a newton method
      Solve_type = Solid_solve;
      newton_solve();

      // Under-relax, either by classical relaxation or by Irons & Tuck
      // Note that we are under-relaxing before the convergence check
      // because under-relaxtion may be required to acheieve any
      // kind of convergence. If the convergence check is on the RELATIVE
      // change of the residuals, however, then a small under-relaxation
      // parameter will cause a false convergence. You have been warned!
      under_relax_solid();

      // Stage 3: Convergence check (possibly again after pointwise Aitken
      //------------------------------------------------------------------
      // extrapolation)
      //---------------
      // Assume that we are not doing Aitken acceleration
      Recheck_convergence_after_pointwise_aitken = false;
      do
      {
        // Perform any actions before the convergence check
        actions_before_segregated_convergence_check();

        // Get the change in the solid variables
        double rms_change;
        double max_change;
        double rms_norm;
        double max_res = 0.0;
        get_solid_change(rms_change, max_change, rms_norm);


        // If we are computing the maximum global residual, do so
        if (get_max_global_residual)
        {
          clock_t t_start = clock();

          // Free all dofs so the residuals of the global equations are computed
          restore_solid_dofs();
          restore_fluid_dofs();
          rebuild_monolithic_mesh();
          assign_eqn_numbers();

          // We are now treating this as a full solve
          Solve_type = Full_solve;

          // Get the residuals
          DoubleVector residual;
          get_residuals(residual);

          // Get maximum residuals, using our own abscmp function
          max_res = residual.max();

          clock_t t_end = clock();
          cpu_for_global_residual += double(t_end - t_start) / CLOCKS_PER_SEC;
        }

        oomph_info << "==================================================\n";
        oomph_info << "Iteration             : " << iter << std::endl;
        oomph_info << "RMS  change           :       " << rms_change
                   << std::endl;
        oomph_info << "Max. change           :       " << max_change
                   << std::endl;
        oomph_info << "RMS norm              :       " << rms_norm << std::endl;
        if (get_max_global_residual)
        {
          oomph_info << "Max. global residual  :       " << max_res
                     << std::endl;
        }
        oomph_info << "==================================================\n\n";

        // Check for convergence
        switch (Convergence_criterion)
        {
          case Assess_convergence_based_on_absolute_solid_change:
            tol_achieved = max_change;
            if (tol_achieved < Convergence_tolerance)
            {
              oomph_info
                << "\n\n\n/////////////////////////////////////////////////////"
                   "///"
                << "\nPicard iteration converged after " << iter << " steps!"
                << std::endl
                << "Convergence was based on absolute change in solid dofs \n"
                << "being less than " << Convergence_tolerance << "."
                << std::endl
                << "////////////////////////////////////////////////////////"
                   "\n\n\n"
                << std::endl;
              converged = true;
            }
            break;


          case Assess_convergence_based_on_relative_solid_change:
            tol_achieved = std::fabs(rms_change / rms_norm);
            if (tol_achieved < Convergence_tolerance)
            {
              oomph_info
                << "\n\n\n/////////////////////////////////////////////////////"
                   "//"
                << "\nPicard iteration converged after " << iter << " steps!"
                << std::endl
                << "Convergence was based on relative change in solid dofs \n"
                << "being less than " << Convergence_tolerance << "."
                << std::endl
                << "////////////////////////////////////////////////////////"
                   "\n\n\n"
                << std::endl;
              converged = true;
            }
            break;

          case Assess_convergence_based_on_max_global_residual:
            tol_achieved = max_res;
            if (tol_achieved < Convergence_tolerance)
            {
              oomph_info
                << "\n\n\n/////////////////////////////////////////////////////"
                   "///"
                << "\nPicard iteration converged after " << iter << " steps!"
                << std::endl
                << "Convergence was based on max. residual of coupled eqns \n"
                << "being less than " << Convergence_tolerance << "."
                << std::endl
                << "////////////////////////////////////////////////////////"
                   "\n\n\n"
                << std::endl;
              converged = true;
            }
            break;
        }

        // If converged bail out
        if (converged)
        {
          // Number of iterations taken
          iter_taken = iter;

          // We have converged (overwrites default of false)
          conv_data.set_solver_converged();

          // Break loop using a GO TO! This is the first (and hopefully
          // the last) one in the whole of oomph-lib. Here it's
          // it's actually the cleanest way to exit from these
          // nested control structures
          goto jump_out_of_picard;
        }

        // Store the current values of the solid dofs as reference values
        // and for the pointwise Aitken extrapolation
        store_solid_dofs();

        // Correct wall shape by pointwise Aitken extrapolation if possible
        //-----------------------------------------------------------------
        // and desired
        //------------
        // This is an acceleration method for the (possibly under-relaxed)
        // sequence of iterates.
        if ((3 == Pointwise_aitken_counter) && (Use_pointwise_aitken) &&
            (iter > Pointwise_aitken_start))
        {
          pointwise_aitken_extrapolate();

          // Repeat the convergence check
          Recheck_convergence_after_pointwise_aitken = true;
        }
        else
        {
          // Didn't change anything: Don't repeat the convergence check
          Recheck_convergence_after_pointwise_aitken = false;
        }
      }
      // Repeat convergence while we are doing aitken extrapolation
      while (Recheck_convergence_after_pointwise_aitken);

    } // End of loop over iterations


    // Jump address for converged or diverged iteration
  jump_out_of_picard:

    // Reset everything
    restore_solid_dofs();
    restore_fluid_dofs();
    rebuild_monolithic_mesh();
    assign_eqn_numbers();
    // This is a full solve again
    Solve_type = Full_solve;

    // Do any updates that are required
    actions_after_segregated_solve();

    // Number of iterations (either this is still Max_iter from
    // the initialisation or it's been overwritten on convergence)
    conv_data.niter() = iter_taken;

    // Total cpu time
    clock_t t_total_end = clock();
    conv_data.cpu_total() =
      double(t_total_end - t_total_start) / CLOCKS_PER_SEC;

    // Total essential cpu time (exluding doc etc)
    conv_data.essential_cpu_total() = t_spent_on_actual_solve();

    // cpu time for check/doc of global residual
    conv_data.cpu_for_global_residual() = cpu_for_global_residual;

    // Final tolerance achieved by the iteration
    conv_data.tol_achieved() = tol_achieved;

    // Doc non-convergence
    if (!converged)
    {
      switch (Convergence_criterion)
      {
        case Assess_convergence_based_on_absolute_solid_change:
          oomph_info
            << "\n\n\n////////////////////////////////////////////////////////"
            << "\nPicard iteration did not converge after " << iter_taken
            << " steps!" << std::endl
            << "Convergence was based on absolute change in solid dofs \n"
            << "being less than " << Convergence_tolerance << " " << std::endl
            << "but we achieved only " << tol_achieved << "." << std::endl
            << "////////////////////////////////////////////////////////\n\n\n"
            << std::endl;
          // Throw an error indicating if we ran out of iterations
          if (iter_taken == Max_picard)
          {
            throw SegregatedSolverError(true);
          }
          else
          {
            throw SegregatedSolverError(false);
          }
          break;


        case Assess_convergence_based_on_relative_solid_change:
          oomph_info
            << "\n\n\n///////////////////////////////////////////////////////"
            << "\nPicard iteration did not converge after " << iter_taken
            << " steps!" << std::endl
            << "Convergence was based on relative change in solid dofs \n"
            << "being less than " << Convergence_tolerance << " " << std::endl
            << "but we achieved only " << tol_achieved << "." << std::endl
            << "////////////////////////////////////////////////////////\n\n\n"
            << std::endl;
          // Throw an error indicating if we ran out of iterations
          if (iter_taken == Max_picard)
          {
            throw SegregatedSolverError(true);
          }
          else
          {
            throw SegregatedSolverError(false);
          }
          break;

        case Assess_convergence_based_on_max_global_residual:
          oomph_info
            << "\n\n\n////////////////////////////////////////////////////////"
            << "\nPicard iteration did not converge after " << iter_taken
            << " steps!" << std::endl
            << "Convergence was based on max. residual of coupled eqns \n"
            << "being less than " << Convergence_tolerance << " " << std::endl
            << "but we achieved only " << tol_achieved << "." << std::endl

            << "////////////////////////////////////////////////////////\n\n\n"
            << std::endl;
          // Throw an error indicating if we ran out of iterations
          if (iter_taken == Max_picard)
          {
            throw SegregatedSolverError(true);
          }
          else
          {
            throw SegregatedSolverError(false);
          }
          break;
      }
    }

    return conv_data;
  }


  //============================================================================
  /// Segregated fixed point solver with optional pointwise Aitken acceleration
  /// on selected solid dofs. Returns PicardConvergenceData object
  /// that contains the vital stats of the iteration.
  //============================================================================
  PicardConvergenceData SegregatableFSIProblem::steady_segregated_solve()
  {
    // Find out how many timesteppers there are
    unsigned n_time_steppers = ntime_stepper();

    // Vector of bools to store the is_steady status of the various
    // timesteppers when we came in here
    std::vector<bool> was_steady(n_time_steppers);

    // Loop over them all and make them (temporarily) static
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      was_steady[i] = time_stepper_pt(i)->is_steady();
      time_stepper_pt(i)->make_steady();
    }

    // Create object to doc convergence stats
    PicardConvergenceData conv_data;

    // Solve the non-linear problem by the segregated solver
    try
    {
      conv_data = segregated_solve();
    }
    // Catch any exceptions thrown in the segregated solver
    catch (SegregatedSolverError& error)
    {
      if (!error.Ran_out_of_iterations)
      {
        std::ostringstream error_stream;
        error_stream << "Error occured in Segregated solver. " << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        oomph_info << "Note: Ran out of iterations but continuing anyway"
                   << std::endl;
      }
    }

    // Reset the is_steady status of all timesteppers that
    // weren't already steady when we came in here and reset their
    // weights
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      if (!was_steady[i])
      {
        time_stepper_pt(i)->undo_make_steady();
      }
    }

    // Since we performed a steady solve, the history values
    // now have to be set as if we had performed an impulsive start from
    // the current solution. This ensures that the time-derivatives
    // evaluate to zero even now that the timesteppers have been
    // reactivated.
    assign_initial_values_impulsive();

    // Return the convergence data
    return conv_data;
  }


  //============================================================================
  /// Segregated fixed point solver with optional pointwise Aitken acceleration
  /// on selected solid dofs. Returns PicardConvergenceData object
  /// that contains the vital stats of the iteration.
  //============================================================================
  PicardConvergenceData SegregatableFSIProblem::unsteady_segregated_solve(
    const double& dt)
  {
    // We shift the values, so shift_values is true
    return unsteady_segregated_solve(dt, true);
  }

  //============================================================================
  /// Segregated fixed point solver with optional pointwise Aitken acceleration
  /// on selected solid dofs. Returns PicardConvergenceData object
  /// that contains the vital stats of the iteration.
  //============================================================================
  PicardConvergenceData SegregatableFSIProblem::unsteady_segregated_solve(
    const double& dt, const bool& shift_values)
  {
    // Shift the time values and the dts according to the control flag
    if (shift_values)
    {
      shift_time_values();
    }

    // Advance global time and set current value of dt
    time_pt()->time() += dt;
    time_pt()->dt() = dt;

    // Find out how many timesteppers there are
    unsigned n_time_steppers = ntime_stepper();

    // Loop over them all and set the weights
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      time_stepper_pt(i)->set_weights();
    }

    // Now update anything that needs updating before the timestep
    // This could be time-dependent boundary conditions, for example.
    actions_before_implicit_timestep();

    // Extrapolate the solid data and then update fluid mesh during unsteady run
    extrapolate_solid_data();

    // Create object to doc convergence stats
    PicardConvergenceData conv_data;

    try
    {
      // Solve the non-linear problem for this timestep with Newton's method
      conv_data = segregated_solve();
    }
    // Catch any exceptions thrown in the segregated solver
    catch (SegregatedSolverError& error)
    {
      if (!error.Ran_out_of_iterations)
      {
        std::ostringstream error_stream;
        error_stream << "Error occured in Segregated solver. " << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        oomph_info << "Note: Ran out of iterations but continuing anyway"
                   << std::endl;
      }
    }

    // Now update anything that needs updating after the timestep
    actions_after_implicit_timestep();

    return conv_data;
  }


  //============================================================================
  /// Setup segregated solver, using the information provided by the user
  /// in his/her implementation of the pure virtual function
  /// identify_fluid_and_solid_dofs(...).
  //============================================================================
  void SegregatableFSIProblem::setup_segregated_solver(
    const bool& full_setup_of_fluid_and_solid_dofs)
  {
    // If we are doing a full setup
    if (full_setup_of_fluid_and_solid_dofs)
    {
      // Identify the fluid and solid Data
      Vector<Data*> fluid_data_pt;
      Vector<Data*> solid_data_pt;
      identify_fluid_and_solid_dofs(
        fluid_data_pt, solid_data_pt, Fluid_mesh_pt, Solid_mesh_pt);

      if (Fluid_mesh_pt == 0)
      {
        oomph_info
          << std::endl
          << std::endl
          << "Warning: Your implementation of the pure virtual\n"
          << "         function identify_fluid_and_solid_dofs(...)\n"
          << "         returned a NULL pointer for Fluid_mesh_pt.\n"
          << "         --> The fluid elements will remain activated\n"
          << "         during the solid solve. This is inefficient!\n"
          << "         You should combine all fluid elements into a combined\n"
          << "         mesh and specify this mesh in your\n"
          << "         implementation of \n\n"
          << "         "
             "SegregatableFSIProblem::identify_fluid_and_solid_dofs(...)"
          << std::endl
          << std::endl;
      }


      if (Solid_mesh_pt == 0)
      {
        oomph_info
          << std::endl
          << std::endl
          << "Warning: Your implementation of the pure virtual\n"
          << "         function identify_fluid_and_solid_dofs(...)\n"
          << "         returned a NULL pointer for Solid_mesh_pt.\n"
          << "         --> The solid elements will remain activated\n"
          << "         during the fluid solve. This is inefficient!\n"
          << "         You should combine all solid elements into a combined\n"
          << "         mesh and specify this mesh in your\n"
          << "         implementation of \n\n"
          << "         "
             "SegregatableFSIProblem::identify_fluid_and_solid_dofs(...)"
          << std::endl
          << std::endl;
      }

      // Back up the pointers to the submeshes in the original problem
      // so we can restore the problem when we're done
      unsigned orig_n_sub_mesh = nsub_mesh();
      Orig_sub_mesh_pt.resize(orig_n_sub_mesh);
      for (unsigned i = 0; i < orig_n_sub_mesh; i++)
      {
        Orig_sub_mesh_pt[i] = mesh_pt(i);
      }

      // Fluid
      //------

      // Reset
      Fluid_data_pt.clear();
      Fluid_value_is_pinned.clear();

      // Make space for fluid data items
      unsigned n_fluid_data = fluid_data_pt.size();
      Fluid_data_pt.resize(n_fluid_data);
      Fluid_value_is_pinned.resize(n_fluid_data);

      // Counter for number of fluid values
      unsigned n_fluid_values = 0;

      // Loop over fluid data
      for (unsigned i = 0; i < n_fluid_data; i++)
      {
        // Copy data
        Fluid_data_pt[i] = fluid_data_pt[i];

        // Number of values stored at this Data item:
        unsigned n_value = fluid_data_pt[i]->nvalue();
        Fluid_value_is_pinned[i].resize(n_value);

        // Copy pinned status for all values
        for (unsigned k = 0; k < n_value; k++)
        {
          Fluid_value_is_pinned[i][k] = fluid_data_pt[i]->is_pinned(k);

          // Increment counter for number of fluid values
          n_fluid_values++;
        }
      }


      // Solid
      //------

      // Reset
      Solid_data_pt.clear();
      Solid_value_is_pinned.clear();
      Previous_solid_value.clear();
      Pointwise_aitken_solid_value.clear();
      Del_irons_and_tuck.clear();


      unsigned n_solid_data = solid_data_pt.size();

      // Make space for solid data items
      unsigned nsolid_data = solid_data_pt.size();
      Solid_data_pt.resize(nsolid_data);
      Solid_value_is_pinned.resize(nsolid_data);

      // Counter for number of solid values
      unsigned n_solid_values = 0;

      // Loop over solid data
      for (unsigned i = 0; i < n_solid_data; i++)
      {
        // Copy data
        Solid_data_pt[i] = solid_data_pt[i];

        // Number of values stored at this Data item:
        unsigned n_value = solid_data_pt[i]->nvalue();
        Solid_value_is_pinned[i].resize(n_value);

        // Copy pinned status for all values
        for (unsigned k = 0; k < n_value; k++)
        {
          Solid_value_is_pinned[i][k] = solid_data_pt[i]->is_pinned(k);

          // Increment counter for number of solid values
          n_solid_values++;
        }
      }

      // Make space for previous solid values
      Previous_solid_value.resize(n_solid_values);

      // Allocate storage and initialise Irons and Tuck extrapolation
      if (Use_irons_and_tuck_extrapolation)
      {
        Del_irons_and_tuck.resize(n_solid_values);
      }

      // Make space for pointwise Aitken extrapolation
      if (Use_pointwise_aitken)
      {
        Pointwise_aitken_solid_value.resize(n_solid_values);
        for (unsigned i = 0; i < n_solid_values; i++)
        {
          Pointwise_aitken_solid_value[i].resize(3);
        }
      }

    } // End of full setup


    // Initialise Irons and Tuck relaxation factor
    R_irons_and_tuck = 1.0 - Omega_relax;

    // Initialise Irons and Tuck delta values
    unsigned n_del = Del_irons_and_tuck.size();
    for (unsigned i = 0; i < n_del; i++)
    {
      Del_irons_and_tuck[i] = 1.0e20;
    }

    // Initialise counter for the number of pointwise Aitken values stored
    Pointwise_aitken_counter = 0;
  }
} // namespace oomph
