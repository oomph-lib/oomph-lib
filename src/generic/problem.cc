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

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

#include <list>
#include <algorithm>
#include <string>

#include "oomph_utilities.h"
#include "problem.h"
#include "timesteppers.h"
#include "explicit_timesteppers.h"
#include "generalised_timesteppers.h"
#include "refineable_mesh.h"
#include "triangle_mesh.h"
#include "linear_solver.h"
#include "eigen_solver.h"
#include "assembly_handler.h"
#include "dg_elements.h"
#include "partitioning.h"
#include "spines.h"

// Include to fill in additional_setup_shared_node_scheme() function
#include "refineable_mesh.template.cc"


namespace oomph
{
  //////////////////////////////////////////////////////////////////
  // Non-inline functions for the problem class
  //////////////////////////////////////////////////////////////////

  //=================================================================
  /// The continuation timestepper object
  //=================================================================
  ContinuationStorageScheme Problem::Continuation_time_stepper;

  //================================================================
  /// Constructor: Allocate space for one time stepper
  /// and set all pointers to NULL and set defaults for all
  /// parameters.
  //===============================================================
  Problem::Problem()
    : Mesh_pt(0),
      Time_pt(0),
      Explicit_time_stepper_pt(0),
      Saved_dof_pt(0),
      Default_set_initial_condition_called(false),
      Use_globally_convergent_newton_method(false),
      Empty_actions_before_read_unstructured_meshes_has_been_called(false),
      Empty_actions_after_read_unstructured_meshes_has_been_called(false),
      Store_local_dof_pt_in_elements(false),
      Calculate_hessian_products_analytic(false),
#ifdef OOMPH_HAS_MPI
      Doc_imbalance_in_parallel_assembly(false),
      Use_default_partition_in_load_balance(false),
      Must_recompute_load_balance_for_assembly(true),
      Halo_scheme_pt(0),
#endif
      Relaxation_factor(1.0),
      Newton_solver_tolerance(1.0e-8),
      Max_newton_iterations(10),
      Nnewton_iter_taken(0),
      Max_residuals(10.0),
      Time_adaptive_newton_crash_on_solve_fail(false),
      Jacobian_reuse_is_enabled(false),
      Jacobian_has_been_computed(false),
      Problem_is_nonlinear(true),
      Pause_at_end_of_sparse_assembly(false),
      Doc_time_in_distribute(false),
      Sparse_assembly_method(Perform_assembly_using_vectors_of_pairs),
      Sparse_assemble_with_arrays_initial_allocation(400),
      Sparse_assemble_with_arrays_allocation_increment(150),
      Numerical_zero_for_sparse_assembly(0.0),
      FD_step_used_in_get_hessian_vector_products(1.0e-8),
      Mass_matrix_reuse_is_enabled(false),
      Mass_matrix_has_been_computed(false),
      Discontinuous_element_formulation(false),
      Minimum_dt(1.0e-12),
      Maximum_dt(1.0e12),
      DTSF_max_increase(4.0),
      DTSF_min_decrease(0.8),
      Minimum_dt_but_still_proceed(-1.0),
      Scale_arc_length(true),
      Desired_proportion_of_arc_length(0.5),
      Theta_squared(1.0),
      Sign_of_jacobian(0),
      Continuation_direction(1.0),
      Parameter_derivative(1.0),
      Parameter_current(0.0),
      Use_continuation_timestepper(false),
      Dof_derivative_offset(1),
      Dof_current_offset(2),
      Ds_current(0.0),
      Desired_newton_iterations_ds(5),
      Minimum_ds(1.0e-10),
      Bifurcation_detection(false),
      Bisect_to_find_bifurcation(false),
      First_jacobian_sign_change(false),
      Arc_length_step_taken(false),
      Use_finite_differences_for_continuation_derivatives(false),
#ifdef OOMPH_HAS_MPI
      Dist_problem_matrix_distribution(Uniform_matrix_distribution),
      Parallel_sparse_assemble_previous_allocation(0),
      Problem_has_been_distributed(false),
      Bypass_increase_in_dof_check_during_pruning(false),
      Max_permitted_error_for_halo_check(1.0e-14),
#endif
      Shut_up_in_newton_solve(false),
      Always_take_one_newton_step(false),
      Timestep_reduction_factor_after_nonconvergence(0.5),
      Keep_temporal_error_below_tolerance(true)
  {
    Use_predictor_values_as_initial_guess = false;

    /// Setup terminate helper
    TerminateHelper::setup();

    // By default no submeshes:
    Sub_mesh_pt.resize(0);
    // No timesteppers
    Time_stepper_pt.resize(0);

    // Set the linear solvers, eigensolver and assembly handler
    Linear_solver_pt = Default_linear_solver_pt = new SuperLUSolver;
    Mass_matrix_solver_for_explicit_timestepper_pt = Linear_solver_pt;

    Eigen_solver_pt = Default_eigen_solver_pt = new ARPACK;

    Assembly_handler_pt = Default_assembly_handler_pt = new AssemblyHandler;

    // setup the communicator
#ifdef OOMPH_HAS_MPI
    if (MPI_Helpers::mpi_has_been_initialised())
    {
      Communicator_pt = new OomphCommunicator(MPI_Helpers::communicator_pt());
    }
    else
    {
      Communicator_pt = new OomphCommunicator();
    }
#else
    Communicator_pt = new OomphCommunicator();
#endif

    // just create an empty linear algebra distribution for the
    // DOFs
    // this is setup when assign_eqn_numbers(...) is called.
    Dof_distribution_pt = new LinearAlgebraDistribution;
  }

  //================================================================
  /// Destructor to clean up memory
  //================================================================
  Problem::~Problem()
  {
    // Delete the memory assigned for the global time
    // (it's created on the fly in Problem::add_time_stepper_pt()
    // so we are entitled to delete it.
    if (Time_pt != 0)
    {
      delete Time_pt;
      Time_pt = 0;
    }

    // We're not using the default linear solver,
    // somebody else must have built it, so that person
    // must be in charge of killing it.
    // We can safely delete the defaults, however
    delete Default_linear_solver_pt;

    delete Default_eigen_solver_pt;
    delete Default_assembly_handler_pt;
    delete Communicator_pt;
    delete Dof_distribution_pt;

    // Delete any copies of the problem that have been created for
    // use in adaptive bifurcation tracking.
    // ALH: This will eventually go
    unsigned n_copies = Copy_of_problem_pt.size();
    for (unsigned c = 0; c < n_copies; c++)
    {
      delete Copy_of_problem_pt[c];
    }

    // if this problem has sub meshes then we must delete the Mesh_pt
    if (Sub_mesh_pt.size() != 0)
    {
      Mesh_pt->flush_element_and_node_storage();
      delete Mesh_pt;
    }

    // Since we called the TerminateHelper setup function in the constructor,
    // we need to delete anything that was dynamically allocated (as it's
    // just a namespace and so doesn't have it's own destructor) in the function
    TerminateHelper::clean_up_memory();
  }

  //=================================================================
  /// Setup the count vector that records how many elements contribute
  /// to each degree of freedom. Returns the total number of elements
  /// in the problem
  //=================================================================
  unsigned Problem::setup_element_count_per_dof()
  {
    // Now set the element counter to have the current Dof distribution
    Element_count_per_dof.build(this->Dof_distribution_pt);
    // We need to use the halo scheme (assuming it has been setup)
#ifdef OOMPH_HAS_MPI
    Element_count_per_dof.build_halo_scheme(this->Halo_scheme_pt);
#endif

    // Loop over the elements and count the entries
    // and number of (non-halo) elements
    const unsigned n_element = this->mesh_pt()->nelement();
    unsigned n_non_halo_element_local = 0;
    for (unsigned e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = this->mesh_pt()->element_pt(e);
#ifdef OOMPH_HAS_MPI
      // Ignore halo elements
      if (!elem_pt->is_halo())
      {
#endif
        // Increment the number of non halo elements
        ++n_non_halo_element_local;
        // Now count the number of times the element contributes to a value
        // using the current assembly handler
        unsigned n_var = this->Assembly_handler_pt->ndof(elem_pt);
        for (unsigned n = 0; n < n_var; n++)
        {
          ++Element_count_per_dof.global_value(
            this->Assembly_handler_pt->eqn_number(elem_pt, n));
        }
#ifdef OOMPH_HAS_MPI
      }
#endif
    }

    // Storage for the total number of elements
    unsigned Nelement = 0;

    // Add together all the counts if we are in parallel
#ifdef OOMPH_HAS_MPI
    Element_count_per_dof.sum_all_halo_and_haloed_values();

    // If distributed, find the total number of elements in the problem
    if (this->Problem_has_been_distributed)
    {
      // Need to gather the total number of non halo elements
      MPI_Allreduce(&n_non_halo_element_local,
                    &Nelement,
                    1,
                    MPI_UNSIGNED,
                    MPI_SUM,
                    this->communicator_pt()->mpi_comm());
    }
    // Otherwise the total number is the same on each processor
    else
#endif
    {
      Nelement = n_non_halo_element_local;
    }

    return Nelement;
  }


#ifdef OOMPH_HAS_MPI

  //==================================================================
  /// Setup the halo scheme for the degrees of freedom
  //==================================================================
  void Problem::setup_dof_halo_scheme()
  {
    // Find the number of elements stored on this processor
    const unsigned n_element = this->mesh_pt()->nelement();

    // Work out the all global equations to which this processor
    // contributes
    Vector<unsigned> my_eqns;
    this->get_my_eqns(this->Assembly_handler_pt, 0, n_element - 1, my_eqns);

    // Build the halo scheme, based on the equations to which this
    // processor contributes
    Halo_scheme_pt =
      new DoubleVectorHaloScheme(this->Dof_distribution_pt, my_eqns);

    // Find pointers to all the halo dofs
    // There may be more of these than required by my_eqns
    //(but should not be less)
    std::map<unsigned, double*> halo_data_pt;
    this->get_all_halo_data(halo_data_pt);

    // Now setup the Halo_dofs
    Halo_scheme_pt->setup_halo_dofs(halo_data_pt, this->Halo_dof_pt);
  }

  //==================================================================
  /// Distribute the problem without doc; report stats if required.
  /// Returns actual partitioning used, e.g. for restart.
  //==================================================================
  Vector<unsigned> Problem::distribute(const bool& report_stats)
  {
    // Set dummy doc paramemters
    DocInfo doc_info;
    doc_info.disable_doc();

    // Set the sizes of the input and output vectors
    unsigned n_element = mesh_pt()->nelement();
    Vector<unsigned> element_partition(n_element, 0);

    // Distribute and return partitioning
    return distribute(element_partition, doc_info, report_stats);
  }

  //==================================================================
  /// Distribute the problem according to specified partition.
  /// If all entries in partitioning vector are zero we use METIS
  /// to do the partitioning after all.
  /// Returns actual partitioning used, e.g. for restart.
  //==================================================================
  Vector<unsigned> Problem::distribute(
    const Vector<unsigned>& element_partition, const bool& report_stats)
  {
#ifdef PARANOID
    bool has_non_zero_entry = false;
    unsigned n = element_partition.size();
    for (unsigned i = 0; i < n; i++)
    {
      if (element_partition[i] != 0)
      {
        has_non_zero_entry = true;
        break;
      }
    }
    if (!has_non_zero_entry)
    {
      std::ostringstream warn_message;
      warn_message << "WARNING: All entries in specified partitioning vector \n"
                   << "         are zero -- will ignore this and use METIS\n"
                   << "         to perform the partitioning\n";
      OomphLibWarning(
        warn_message.str(), "Problem::distribute()", OOMPH_EXCEPTION_LOCATION);
    }
#endif
    // Set dummy doc paramemters
    DocInfo doc_info;
    doc_info.disable_doc();

    // Distribute and return partitioning
    return distribute(element_partition, doc_info, report_stats);
  }

  //==================================================================
  /// Distribute the problem and doc to specified DocInfo.
  /// Returns actual partitioning used, e.g. for restart.
  //==================================================================
  Vector<unsigned> Problem::distribute(DocInfo& doc_info,
                                       const bool& report_stats)
  {
    // Set the sizes of the input and output vectors
    unsigned n_element = mesh_pt()->nelement();

    // Dummy input vector
    Vector<unsigned> element_partition(n_element, 0);

    // Distribute and return partitioning
    return distribute(element_partition, doc_info, report_stats);
  }

  //==================================================================
  /// Distribute the problem according to specified partition.
  /// (If all entries in partitioning vector are zero we use METIS
  /// to do the partitioning after all) and doc.
  /// Returns actual partitioning used, e.g. for restart.
  //==================================================================
  Vector<unsigned> Problem::distribute(
    const Vector<unsigned>& element_partition,
    DocInfo& doc_info,
    const bool& report_stats)
  {
    // Storage for number of processors and number of elements in global mesh
    int n_proc = this->communicator_pt()->nproc();
    int my_rank = this->communicator_pt()->my_rank();
    int n_element = mesh_pt()->nelement();

    // Vector to be returned
    Vector<unsigned> return_element_domain;

    // Buffer extreme cases
    if (n_proc == 1) // single-process job - don't do anything
    {
      if (report_stats)
      {
        std::ostringstream warn_message;
        warn_message << "WARNING: You've tried to distribute a problem over\n"
                     << "only one processor: this would make METIS crash.\n"
                     << "Ignoring your request for distribution.\n";
        OomphLibWarning(warn_message.str(),
                        "Problem::distribute()",
                        OOMPH_EXCEPTION_LOCATION);
      }
    }
    else if (n_proc > n_element) // more processors than elements
    {
      // Throw an error
      std::ostringstream error_stream;
      error_stream << "You have tried to distribute a problem\n"
                   << "but there are less elements than processors.\n"
                   << "Please re-run with more elements!\n"
                   << "Please also ensure that actions_before_distribute().\n"
                   << "and actions_after_distribute() are correctly set up.\n"
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      // We only distribute uniformly-refined meshes; buffer the case where
      // either mesh is not uniformly refined
      bool a_mesh_is_not_uniformly_refined = false;
      unsigned n_mesh = nsub_mesh();
      if (n_mesh == 0)
      {
        // Check refinement levels
        if (TreeBasedRefineableMeshBase* mmesh_pt =
              dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(0)))
        {
          unsigned min_ref_level = 0;
          unsigned max_ref_level = 0;
          mmesh_pt->get_refinement_levels(min_ref_level, max_ref_level);
          // If they are not the same
          if (max_ref_level != min_ref_level)
          {
            a_mesh_is_not_uniformly_refined = true;
          }
        }
      }
      else
      {
        for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
        {
          // Check refinement levels for each mesh individually
          // (one mesh is allowed to be "more uniformly refined" than another)
          if (TreeBasedRefineableMeshBase* mmesh_pt =
                dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(i_mesh)))
          {
            unsigned min_ref_level = 0;
            unsigned max_ref_level = 0;
            mmesh_pt->get_refinement_levels(min_ref_level, max_ref_level);
            // If they are not the same
            if (max_ref_level != min_ref_level)
            {
              a_mesh_is_not_uniformly_refined = true;
            }
          }
        }
      }

      // If any mesh is not uniformly refined
      if (a_mesh_is_not_uniformly_refined)
      {
        // Again it may make more sense to throw an error here as the user
        // will probably not be running a problem that is small enough to
        // fit the whole of on each processor
        std::ostringstream error_stream;
        error_stream << "You have tried to distribute a problem\n"
                     << "but at least one of your meshes is no longer\n"
                     << "uniformly refined.  In order to preserve the Tree\n"
                     << "and TreeForest structure, Problem::distribute() can\n"
                     << "only be called while meshes are uniformly refined.\n"
                     << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        // Is there any global data?  If so, distributing the problem won't work
        if (nglobal_data() > 0)
        {
          std::ostringstream error_stream;
          error_stream << "You have tried to distribute a problem\n"
                       << "and there is some global data.\n"
                       << "This is not likely to work...\n"
                       << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        double t_start = 0;
        if (Doc_time_in_distribute)
        {
          t_start = TimingHelpers::timer();
        }


#ifdef PARANOID
        unsigned old_ndof = ndof();
#endif

        // Need to partition the global mesh before distributing
        Mesh* global_mesh_pt = mesh_pt();

        // Vector listing the affiliation of each element
        unsigned nelem = global_mesh_pt->nelement();
        Vector<unsigned> element_domain(nelem);

        // Number of elements that I'm in charge of, based on any
        // incoming partitioning
        unsigned n_my_elements = 0;

        // Have we used the pre-set partitioning
        bool used_preset_partitioning = false;

        // Partition the mesh, unless the partition has already been passed in
        // If it hasn't then the sum of all the entries of the vector should be
        // 0
        unsigned sum_element_partition = 0;
        unsigned n_part = element_partition.size();
        for (unsigned e = 0; e < n_part; e++)
        {
          // ... another one for me.
          if (int(element_partition[e]) == my_rank) n_my_elements++;

          sum_element_partition += element_partition[e];
        }
        if (sum_element_partition == 0)
        {
          oomph_info << "INFO: using METIS to partition elements" << std::endl;
          partition_global_mesh(global_mesh_pt, doc_info, element_domain);
          used_preset_partitioning = false;
        }
        else
        {
          oomph_info << "INFO: using pre-set partition of elements"
                     << std::endl;
          used_preset_partitioning = true;
          element_domain = element_partition;
        }

        // Set the GLOBAL Mesh as being distributed
        global_mesh_pt->set_communicator_pt(this->communicator_pt());

        double t_end = 0.0;
        if (Doc_time_in_distribute)
        {
          t_end = TimingHelpers::timer();
          oomph_info << "Time for partitioning of global mesh: "
                     << t_end - t_start << std::endl;
          t_start = TimingHelpers::timer();
        }

        // Store how many elements we had in the various sub-meshes
        // before actions_before_distribute() (which may empty some of
        // them).
        Vector<unsigned> n_element_in_old_submesh(n_mesh);
        if (n_mesh != 0)
        {
          for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
          {
            unsigned nsub_elem = mesh_pt(i_mesh)->nelement();
            n_element_in_old_submesh[i_mesh] = nsub_elem;
          }
        }

        // Partitioning complete; call actions before distribute
        actions_before_distribute();

        if (Doc_time_in_distribute)
        {
          t_end = TimingHelpers::timer();
          oomph_info << "Time for actions before distribute: "
                     << t_end - t_start << std::endl;
        }

        // This next bit is cheap -- omit timing
        // t_start = TimingHelpers::timer();

        // Number of submeshes (NB: some may have been deleted in
        //                          actions_after_distribute())
        n_mesh = nsub_mesh();


        // Prepare vector of vectors for submesh element domains
        Vector<Vector<unsigned>> submesh_element_domain(n_mesh);

        // The submeshes need to know their own element domains.
        // Also if any meshes have been emptied we ignore their
        // partitioning in the vector that we return from here
        return_element_domain.reserve(element_domain.size());
        if (n_mesh != 0)
        {
          unsigned count = 0;
          for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
          {
            unsigned nsub_elem = mesh_pt(i_mesh)->nelement();
            submesh_element_domain[i_mesh].resize(nsub_elem);
            unsigned nsub_elem_old = n_element_in_old_submesh[i_mesh];
            for (unsigned e = 0; e < nsub_elem_old; e++)
            {
              if (nsub_elem_old == nsub_elem)
              {
                submesh_element_domain[i_mesh][e] = element_domain[count];
                return_element_domain.push_back(element_domain[count]);
              }
              // return_element_domain.push_back(element_domain[count]);
              count++;
            }
          }
        }
        else
        {
          return_element_domain = element_domain;
        }

        if (Doc_time_in_distribute)
        {
          t_start = TimingHelpers::timer();
        }

        // Setup the map between "root" element and number in global mesh
        // (currently used in the load_balance() routines)

        // This map is only established for structured meshes, then we
        // need to check here the type of mesh
        if (n_mesh == 0)
        {
          // Check if the only one mesh is an structured mesh
          bool structured_mesh = true;
          TriangleMeshBase* tri_mesh_pt =
            dynamic_cast<TriangleMeshBase*>(mesh_pt(0));
          if (tri_mesh_pt != 0)
          {
            structured_mesh = false;
          } // if (tri_mesh_pt != 0)
          if (structured_mesh)
          {
            const unsigned n_ele = global_mesh_pt->nelement();
            Base_mesh_element_pt.resize(n_ele);
            Base_mesh_element_number_plus_one.clear();
            for (unsigned e = 0; e < n_ele; e++)
            {
              GeneralisedElement* el_pt = global_mesh_pt->element_pt(e);
              Base_mesh_element_number_plus_one[el_pt] = e + 1;
              Base_mesh_element_pt[e] = el_pt;
            } // for (e<n_ele)
          } // A TreeBaseMesh mesh
        } // if (n_mesh==0)
        else
        {
          // If we have submeshes then we only add those elements that
          // belong to structured meshes, but first compute the number
          // of total elements in the structured meshes
          unsigned nglobal_element = 0;
          // Store which submeshes are structured
          std::vector<bool> is_structured_mesh(n_mesh);
          for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
          {
            TriangleMeshBase* tri_mesh_pt =
              dynamic_cast<TriangleMeshBase*>(mesh_pt(i_mesh));
            if (tri_mesh_pt != 0)
            {
              // Set the flags to indicate this is not an structured
              // mesh
              is_structured_mesh[i_mesh] = false;
            } // if (tri_mesh_pt != 0)
            else
            {
              // Set the flags to indicate this is an structured
              // mesh
              is_structured_mesh[i_mesh] = true;
            } // else if (tri_mesh_pt != 0)
            // Check if mesh is an structured mesh
            if (is_structured_mesh[i_mesh])
            {
              nglobal_element += mesh_pt(i_mesh)->nelement();
            } // A TreeBaseMesh mesh
          } // for (i_mesh<n_mesh)

          // Once computed the number of elements, then resize the
          // structure
          Base_mesh_element_pt.resize(nglobal_element);
          Base_mesh_element_number_plus_one.clear();
          unsigned counter = 0;
          for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
          {
            // Check if mesh is an structured mesh
            if (is_structured_mesh[i_mesh])
            {
              const unsigned n_ele = mesh_pt(i_mesh)->nelement();
              for (unsigned e = 0; e < n_ele; e++)
              {
                GeneralisedElement* el_pt = mesh_pt(i_mesh)->element_pt(e);
                Base_mesh_element_number_plus_one[el_pt] = counter + 1;
                Base_mesh_element_pt[counter] = el_pt;
                // Inrease the global element number
                counter++;
              } // for (e<n_ele)
            } // An structured mesh
          } // for (i_mesh<n_mesh)

#ifdef PARANOID
          if (counter != nglobal_element)
          {
            std::ostringstream error_stream;
            error_stream
              << "The number of global elements (" << nglobal_element
              << ") is not the sameas the number of\nadded elements ("
              << counter << ") to the Base_mesh_element_pt data "
              << "structure!!!\n\n";
            throw OomphLibError(error_stream.str(),
                                "Problem::distribute()",
                                OOMPH_EXCEPTION_LOCATION);
          } // if (counter != nglobal_element)
#endif // #ifdef PARANOID

        } // else if (n_mesh==0)

        // Wipe everything if a pre-determined partitioning
        // didn't specify ANY elements for this processor
        // (typically happens during restarts with larger number
        // of processors -- in this case we really want an empty
        // processor rather than one with any "kept" halo elements)
        bool overrule_keep_as_halo_element_status = false;
        if ((n_my_elements == 0) && (used_preset_partitioning))
        {
          oomph_info << "INFO: We're over-ruling the \"keep as halo element\"\n"
                     << "      status because the preset partitioning\n"
                     << "      didn't place ANY elements on this processor,\n"
                     << "      probably because of a restart on a larger \n"
                     << "      number of processors\n";
          overrule_keep_as_halo_element_status = true;
        }


        // Distribute the (sub)meshes (i.e. sort out their halo lookup schemes)
        Vector<GeneralisedElement*> deleted_element_pt;
        if (n_mesh == 0)
        {
          global_mesh_pt->distribute(this->communicator_pt(),
                                     element_domain,
                                     deleted_element_pt,
                                     doc_info,
                                     report_stats,
                                     overrule_keep_as_halo_element_status);
        }
        else // There are submeshes, "distribute" each one separately
        {
          for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
          {
            if (report_stats)
            {
              oomph_info << "Distributing submesh " << i_mesh << std::endl
                         << "--------------------" << std::endl;
            }
            // Set the doc_info number to reflect the submesh
            doc_info.number() = i_mesh;
            mesh_pt(i_mesh)->distribute(this->communicator_pt(),
                                        submesh_element_domain[i_mesh],
                                        deleted_element_pt,
                                        doc_info,
                                        report_stats,
                                        overrule_keep_as_halo_element_status);
          }
          // Rebuild the global mesh
          rebuild_global_mesh();
        }

        // Null out information associated with deleted elements
        unsigned n_del = deleted_element_pt.size();
        for (unsigned e = 0; e < n_del; e++)
        {
          GeneralisedElement* el_pt = deleted_element_pt[e];
          unsigned old_el_number = Base_mesh_element_number_plus_one[el_pt] - 1;
          Base_mesh_element_number_plus_one[el_pt] = 0;
          Base_mesh_element_pt[old_el_number] = 0;
        }

        if (Doc_time_in_distribute)
        {
          t_end = TimingHelpers::timer();
          oomph_info << "Time for mesh-level distribution: " << t_end - t_start
                     << std::endl;
          t_start = TimingHelpers::timer();
        }

        // Now the problem has been distributed
        Problem_has_been_distributed = true;

        // Call actions after distribute
        actions_after_distribute();

        if (Doc_time_in_distribute)
        {
          t_end = TimingHelpers::timer();
          oomph_info << "Time for actions after distribute: " << t_end - t_start
                     << std::endl;
          t_start = TimingHelpers::timer();
        }

        // Re-assign the equation numbers (incl synchronisation if reqd)
        unsigned n_dof = assign_eqn_numbers();
        oomph_info << "Number of equations: " << n_dof << std::endl;

        if (Doc_time_in_distribute)
        {
          t_end = TimingHelpers::timer();
          oomph_info << "Time for re-assigning eqn numbers (in distribute): "
                     << t_end - t_start << std::endl;
        }


#ifdef PARANOID
        if (n_dof != old_ndof)
        {
          std::ostringstream error_stream;
          error_stream
            << "Number of dofs in distribute() has changed "
            << "from " << old_ndof << " to " << n_dof << "\n"
            << "Check that you've implemented any necessary "
               "actions_before/after\n"
            << "distribute functions, e.g. to pin redundant pressure dofs"
            << " etc.\n";
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      } // end if to check for uniformly refined mesh(es)

    } // end if to check number of processors vs. number of elements etc.


    // Force re-analysis of time spent on assembly each
    // elemental Jacobian
    Must_recompute_load_balance_for_assembly = true;
    Elemental_assembly_time.clear();

    // Return the partition vector used in the distribution
    return return_element_domain;
  }

  //==================================================================
  /// Partition the global mesh, return vector specifying the processor
  /// number for each element. Virtual so that it can be overloaded by
  /// any user; the default is to use METIS to perform the partitioning
  /// (with a bit of cleaning up afterwards to sort out "special cases").
  //==================================================================
  void Problem::partition_global_mesh(Mesh*& global_mesh_pt,
                                      DocInfo& doc_info,
                                      Vector<unsigned>& element_domain,
                                      const bool& report_stats)
  {
    // Storage for number of processors and current processor
    int n_proc = this->communicator_pt()->nproc();
    int rank = this->communicator_pt()->my_rank();

    std::ostringstream filename;
    std::ofstream some_file;

    // Doc the original mesh on proc 0
    //--------------------------------
    if (doc_info.is_doc_enabled())
    {
      if (rank == 0)
      {
        filename << doc_info.directory() << "/complete_mesh"
                 << doc_info.number() << ".dat";
        global_mesh_pt->output(filename.str().c_str(), 5);
      }
    }

    // Partition the mesh
    //-------------------
    // METIS Objective (0: minimise edge cut; 1: minimise total comm volume)
    unsigned objective = 0;

    // Do the partitioning
    unsigned nelem = 0;
    if (this->communicator_pt()->my_rank() == 0)
    {
      METIS::partition_mesh(this, n_proc, objective, element_domain);
      nelem = element_domain.size();
    }
    MPI_Bcast(&nelem, 1, MPI_UNSIGNED, 0, this->communicator_pt()->mpi_comm());
    element_domain.resize(nelem);
    MPI_Bcast(&element_domain[0],
              nelem,
              MPI_UNSIGNED,
              0,
              this->communicator_pt()->mpi_comm());

    // On very coarse meshes with larger numbers of processors, METIS
    // occasionally returns an element_domain Vector for which a particular
    // processor has no elements affiliated to it; the following fixes this

    // Convert element_domain to integer storage
    Vector<int> int_element_domain(nelem);
    for (unsigned e = 0; e < nelem; e++)
    {
      int_element_domain[e] = element_domain[e];
    }

    // Global storage for number of elements on each process
    int my_number_of_elements = 0;
    Vector<int> number_of_elements(n_proc, 0);

    for (unsigned e = 0; e < nelem; e++)
    {
      if (int_element_domain[e] == rank)
      {
        my_number_of_elements++;
      }
    }

    // Communicate the correct value for each single process into
    // the global storage vector
    MPI_Allgather(&my_number_of_elements,
                  1,
                  MPI_INT,
                  &number_of_elements[0],
                  1,
                  MPI_INT,
                  this->communicator_pt()->mpi_comm());

    // If a process has no elements then switch an element with the
    // process with the largest number of elements, assuming
    // that it still has enough elements left to share
    int max_number_of_elements = 0;
    int process_with_max_elements = 0;
    for (int d = 0; d < n_proc; d++)
    {
      if (number_of_elements[d] == 0)
      {
        // Find the process with maximum number of elements
        if (max_number_of_elements <= 1)
        {
          for (int dd = 0; dd < n_proc; dd++)
          {
            if (number_of_elements[dd] > max_number_of_elements)
            {
              max_number_of_elements = number_of_elements[dd];
              process_with_max_elements = dd;
            }
          }
        }

        // Check that this number of elements is okay for sharing
        if (max_number_of_elements <= 1)
        {
          // Throw an error if elements can't be shared
          std::ostringstream error_stream;
          error_stream << "No process has more than 1 element, and\n"
                       << "at least one process has no elements!\n"
                       << "Suggest rerunning with more refinement.\n"
                       << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        // Loop over the element domain vector and switch
        // one value for process "process_with_max_elements" with d
        for (unsigned e = 0; e < nelem; e++)
        {
          if (int_element_domain[e] == process_with_max_elements)
          {
            int_element_domain[e] = d;
            // Change the numbers associated with these processes
            number_of_elements[d]++;
            number_of_elements[process_with_max_elements]--;
            // Reduce the number of elements available on "max" process
            max_number_of_elements--;
            // Inform the user that a switch has taken place
            if (report_stats)
            {
              oomph_info << "INFO: Switched element domain at position " << e
                         << std::endl
                         << "from process " << process_with_max_elements
                         << " to process " << d << std::endl
                         << "which was given no elements by METIS partition"
                         << std::endl;
            }
            // Only need to do this once for this element loop, otherwise
            // this will take all the elements from "max" process and put them
            // in process d, thus leaving essentially the same problem!
            break;
          }
        }
      }
    }

    // Reassign new values to the element_domain vector
    for (unsigned e = 0; e < nelem; e++)
    {
      element_domain[e] = int_element_domain[e];
    }

    unsigned count_elements = 0;
    for (unsigned e = 0; e < nelem; e++)
    {
      if (int(element_domain[e]) == rank)
      {
        count_elements++;
      }
    }

    if (report_stats)
    {
      oomph_info << "I have " << count_elements
                 << " elements from this partition" << std::endl
                 << std::endl;
    }
  }

  //==================================================================
  /// (Irreversibly) prune halo(ed) elements and nodes, usually
  /// after another round of refinement, to get rid of
  /// excessively wide halo layers. Note that the current
  /// mesh will be now regarded as the base mesh and no unrefinement
  /// relative to it will be possible once this function
  /// has been called.
  //==================================================================
  void Problem::prune_halo_elements_and_nodes(DocInfo& doc_info,
                                              const bool& report_stats)
  {
    // Storage for number of processors and current processor
    int n_proc = this->communicator_pt()->nproc();

    // Has the problem been distributed yet?
    if (!Problem_has_been_distributed)
    {
      oomph_info
        << "WARNING: Problem::prune_halo_elements_and_nodes() was called on a "
        << "non-distributed Problem!" << std::endl;
      oomph_info << "Ignoring your request..." << std::endl;
    }
    else
    {
      // There are no halo layers to prune if it's a single-process job
      if (n_proc == 1)
      {
        oomph_info
          << "WARNING: You've tried to prune halo layers on a problem\n"
          << "with only one processor: this is unnecessary.\n"
          << "Ignoring your request." << std::endl
          << std::endl;
      }
      else
      {
#ifdef PARANOID
        unsigned old_ndof = ndof();
#endif

        double t_start = 0.0;
        if (Global_timings::Doc_comprehensive_timings)
        {
          t_start = TimingHelpers::timer();
        }

        // Call actions before distribute
        actions_before_distribute();

        double t_end = 0.0;
        if (Global_timings::Doc_comprehensive_timings)
        {
          t_end = TimingHelpers::timer();
          oomph_info << "Time for actions_before_distribute() in "
                     << "Problem::prune_halo_elements_and_nodes(): "
                     << t_end - t_start << std::endl;
          t_start = TimingHelpers::timer();
        }

        // Associate all elements with root in current Base mesh
        unsigned nel = Base_mesh_element_pt.size();
        std::map<GeneralisedElement*, unsigned>
          old_base_element_number_plus_one;
        std::vector<bool> old_root_is_halo_or_non_existent(nel, true);
        for (unsigned e = 0; e < nel; e++)
        {
          // Get the base element
          GeneralisedElement* base_el_pt = Base_mesh_element_pt[e];

          // Does it exist locally?
          if (base_el_pt != 0)
          {
            // Check if it's a halo element
            if (!base_el_pt->is_halo())
            {
              old_root_is_halo_or_non_existent[e] = false;
            }

            // Not refineable: It's only the element iself
            RefineableElement* ref_el_pt = 0;
            ref_el_pt = dynamic_cast<RefineableElement*>(base_el_pt);
            if (ref_el_pt == 0)
            {
              old_base_element_number_plus_one[base_el_pt] = e + 1;
            }
            // Refineable: Get entire tree of elements
            else
            {
              Vector<Tree*> tree_pt;
              ref_el_pt->tree_pt()->stick_all_tree_nodes_into_vector(tree_pt);
              unsigned ntree = tree_pt.size();
              for (unsigned t = 0; t < ntree; t++)
              {
                old_base_element_number_plus_one[tree_pt[t]->object_pt()] =
                  e + 1;
              }
            }
          }
        }


        if (Global_timings::Doc_comprehensive_timings)
        {
          t_end = TimingHelpers::timer();
          oomph_info << "Time for setup old root elements in "
                     << "Problem::prune_halo_elements_and_nodes(): "
                     << t_end - t_start << std::endl;
          t_start = TimingHelpers::timer();
        }


        // Now remember the old number of base elements
        unsigned nel_base_old = nel;


        // Prune the halo elements and nodes of the mesh(es)
        Vector<GeneralisedElement*> deleted_element_pt;
        unsigned n_mesh = nsub_mesh();
        if (n_mesh == 0)
        {
          // Prune halo elements and nodes for the (single) global mesh
          mesh_pt()->prune_halo_elements_and_nodes(
            deleted_element_pt, doc_info, report_stats);
        }
        else
        {
          // Loop over individual submeshes and prune separately
          for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
          {
            mesh_pt(i_mesh)->prune_halo_elements_and_nodes(
              deleted_element_pt, doc_info, report_stats);
          }

          // Rebuild the global mesh
          rebuild_global_mesh();
        }

        if (Global_timings::Doc_comprehensive_timings)
        {
          t_end = TimingHelpers::timer();
          oomph_info << "Total time for all mesh-level prunes in "
                     << "Problem::prune_halo_elements_and_nodes(): "
                     << t_end - t_start << std::endl;
          t_start = TimingHelpers::timer();
        }

        // Loop over all elements in newly rebuilt mesh (which contains
        // all element in "tree order"), find the roots
        // (which are either non-refineable elements or refineable elements
        // whose tree representations are TreeRoots)
        std::map<FiniteElement*, bool> root_el_done;

        // Vector storing vectors of pointers to new base elements associated
        // with the same old base element
        Vector<Vector<GeneralisedElement*>>
          new_base_element_associated_with_old_base_element(nel_base_old);

        unsigned n_meshes = n_mesh;
        // Change the value for the number of submeshes if there is only
        // one mesh so that the loop below works if we have only one
        // mesh
        if (n_meshes == 0)
        {
          n_meshes = 1;
        }

        // Store which submeshes, if there are some are structured
        // meshes
        std::vector<bool> is_structured_mesh(n_meshes);

        // Loop over all elements in the rebuilt mesh, but make sure
        // that we are only looping over the structured meshes
        nel = 0;
        for (unsigned i_mesh = 0; i_mesh < n_meshes; i_mesh++)
        {
          TriangleMeshBase* tri_mesh_pt =
            dynamic_cast<TriangleMeshBase*>(mesh_pt(i_mesh));
          if (!(tri_mesh_pt != 0))
          {
            // Mark the mesh as structured mesh
            is_structured_mesh[i_mesh] = true;
            // Add the number of elements
            nel += mesh_pt(i_mesh)->nelement();
          } // if (!(tri_mesh_pt!=0))
          else
          {
            // Mark the mesh as nonstructured mesh
            is_structured_mesh[i_mesh] = false;
          } // else if (!(tri_mesh_pt!=0))
        } // for (i_mesh < n_mesh)

        // Go for all the meshes (if there are submeshes)
        for (unsigned i_mesh = 0; i_mesh < n_meshes; i_mesh++)
        {
          // Only work with the elements in the mesh if it is an
          // structured mesh
          if (is_structured_mesh[i_mesh])
          {
            // Get the number of elements in the submesh
            const unsigned nele_submesh = mesh_pt(i_mesh)->nelement();
            for (unsigned e = 0; e < nele_submesh; e++)
            {
              // Get the element
              GeneralisedElement* el_pt = mesh_pt(i_mesh)->element_pt(e);

              // Not refineable: It's definitely a new base element
              RefineableElement* ref_el_pt = 0;
              ref_el_pt = dynamic_cast<RefineableElement*>(el_pt);
              if (ref_el_pt == 0)
              {
                unsigned old_base_el_no =
                  old_base_element_number_plus_one[el_pt] - 1;
                new_base_element_associated_with_old_base_element
                  [old_base_el_no]
                    .push_back(el_pt);
              }
              // Refineable
              else
              {
                // Is it a tree root (after pruning)? In that case it's
                // a new base element
                if (dynamic_cast<TreeRoot*>(ref_el_pt->tree_pt()))
                {
                  unsigned old_base_el_no =
                    old_base_element_number_plus_one[el_pt] - 1;
                  new_base_element_associated_with_old_base_element
                    [old_base_el_no]
                      .push_back(el_pt);
                }
                else
                {
                  // Get associated root element
                  FiniteElement* root_el_pt =
                    ref_el_pt->tree_pt()->root_pt()->object_pt();

                  if (!root_el_done[root_el_pt])
                  {
                    root_el_done[root_el_pt] = true;
                    unsigned old_base_el_no =
                      old_base_element_number_plus_one[el_pt] - 1;
                    new_base_element_associated_with_old_base_element
                      [old_base_el_no]
                        .push_back(root_el_pt);
                  }
                }
              }
            } // for (e < nele_submesh)
          } // if (is_structured_mesh[i_mesh])
        } // for (i_mesh < n_mesh)

        // Create a vector that stores how many new root/base elements
        // got spawned from each old root/base element in the global mesh
        Vector<unsigned> local_n_new_root(nel_base_old);
#ifdef PARANOID
        Vector<unsigned> n_new_root_back(nel_base_old);
#endif
        for (unsigned e = 0; e < nel_base_old; e++)
        {
          local_n_new_root[e] =
            new_base_element_associated_with_old_base_element[e].size();

#ifdef PARANOID
          // Backup so we can check that halo data was consistent
          n_new_root_back[e] = local_n_new_root[e];
#endif
        }

        if (Global_timings::Doc_comprehensive_timings)
        {
          t_end = TimingHelpers::timer();
          oomph_info << "Time for setup of new base elements in "
                     << "Problem::prune_halo_elements_and_nodes(): "
                     << t_end - t_start << std::endl;
          t_start = TimingHelpers::timer();
        }

        // Now do reduce operation to get information for all
        // old root/base elements -- the pruned (halo!) base elements contain
        // fewer associated new roots.
        Vector<unsigned> n_new_root(nel_base_old);
        MPI_Allreduce(&local_n_new_root[0],
                      &n_new_root[0],
                      nel_base_old,
                      MPI_UNSIGNED,
                      MPI_MAX,
                      this->communicator_pt()->mpi_comm());


        if (Global_timings::Doc_comprehensive_timings)
        {
          t_end = TimingHelpers::timer();
          oomph_info << "Time for allreduce in "
                     << "Problem::prune_halo_elements_and_nodes(): "
                     << t_end - t_start << std::endl;
          t_start = TimingHelpers::timer();
        }

        // Find out total number of base elements
        unsigned nel_base_new = 0;
        for (unsigned e = 0; e < nel_base_old; e++)
        {
          // Increment
          nel_base_new += n_new_root[e];

#ifdef PARANOID
          // If we already had data for this root previously then
          // the data ought to be consistent afterwards (since taking
          // the max of consistent numbers shouldn't change things -- this
          // deals with halo/haloed elements)
          if (!old_root_is_halo_or_non_existent[e])
          {
            if (n_new_root_back[e] != 0)
            {
              if (n_new_root_back[e] != n_new_root[e])
              {
                std::ostringstream error_stream;
                error_stream
                  << "Number of new root elements spawned from old root " << e
                  << ": " << n_new_root[e] << "\nis not consistent"
                  << " with previous value: " << n_new_root_back[e]
                  << std::endl;
                throw OomphLibError(error_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
            }
          }

#endif
        }

        // Reset base_mesh information
        Base_mesh_element_pt.clear();
        Base_mesh_element_pt.resize(nel_base_new, 0);
        Base_mesh_element_number_plus_one.clear();

        // Now enumerate the new base/root elements consistently
        unsigned count = 0;
        for (unsigned e = 0; e < nel_base_old; e++)
        {
          // Old root is non-halo: Just add the new roots into the
          // new lookup scheme consecutively
          if (!old_root_is_halo_or_non_existent[e])
          {
            // Loop over new root/base element
            unsigned n_new_root =
              new_base_element_associated_with_old_base_element[e].size();
            for (unsigned j = 0; j < n_new_root; j++)
            {
              // Store new root/base element
              GeneralisedElement* el_pt =
                new_base_element_associated_with_old_base_element[e][j];
              Base_mesh_element_pt[count] = el_pt;
              Base_mesh_element_number_plus_one[el_pt] = count + 1;

              // Bump counter
              count++;
            }
          }
          // Old root element is halo so skip insertion (i.e. leave
          // entries in lookup schemes nulled) but increase counter to
          // ensure consistency between processors
          else
          {
            unsigned nskip = n_new_root[e];
            count += nskip;
          }
        }

        // Re-setup the map between "root" element and number in global mesh
        // (used in the load_balance() routines)
        setup_base_mesh_info_after_pruning();


        if (Global_timings::Doc_comprehensive_timings)
        {
          t_end = TimingHelpers::timer();
          oomph_info << "Time for finishing off base mesh info "
                     << "Problem::prune_halo_elements_and_nodes(): "
                     << t_end - t_start << std::endl;
          t_start = TimingHelpers::timer();
        }


        // Call actions after distribute
        actions_after_distribute();


        if (Global_timings::Doc_comprehensive_timings)
        {
          t_end = TimingHelpers::timer();
          oomph_info << "Time for actions_after_distribute() "
                     << "Problem::prune_halo_elements_and_nodes(): "
                     << t_end - t_start << std::endl;
          t_start = TimingHelpers::timer();
        }


        // Re-assign the equation numbers (incl synchronisation if reqd)
#ifdef PARANOID
        unsigned n_dof = assign_eqn_numbers();
#else
        assign_eqn_numbers();
#endif


        if (Global_timings::Doc_comprehensive_timings)
        {
          t_end = TimingHelpers::timer();
          oomph_info << "Time for assign_eqn_numbers() "
                     << "Problem::prune_halo_elements_and_nodes(): "
                     << t_end - t_start << std::endl;
          t_start = TimingHelpers::timer();
        }


#ifdef PARANOID
        if (!Bypass_increase_in_dof_check_during_pruning)
        {
          if (n_dof != old_ndof)
          {
            std::ostringstream error_stream;
            error_stream
              << "Number of dofs in prune_halo_elements_and_nodes() has "
                 "changed "
              << "from " << old_ndof << " to " << n_dof << "\n"
              << "Check that you've implemented any necessary "
                 "actions_before/after"
              << "\nadapt/distribute functions, e.g. to pin redundant pressure"
              << " dofs etc.\n";
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
#endif
      }
    }
  }


#endif


  //===================================================================
  /// Build a single (global) mesh from a number
  /// of submeshes which are passed as a vector of pointers to the
  /// submeshes. The ordering is not necessarily optimal.
  //==============================================================
  void Problem::build_global_mesh()
  {
#ifdef PARANOID
    // Has a global mesh already been built
    if (Mesh_pt != 0)
    {
      std::string error_message = "Problem::build_global_mesh() called,\n";
      error_message += " but a global mesh has already been built:\n";
      error_message += "Problem::Mesh_pt is not zero!\n";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // Check that there are submeshes
    if (Sub_mesh_pt.size() == 0)
    {
      std::string error_message = "Problem::build_global_mesh() called,\n";
      error_message += " but there are no submeshes:\n";
      error_message += "Problem::Sub_mesh_pt has no entries\n";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Create an empty mesh
    Mesh_pt = new Mesh();

    // Call the rebuild function to construct the mesh
    rebuild_global_mesh();
  }

  //====================================================================
  /// If one of the submeshes has changed (e.g. by
  /// mesh adaptation) we need to update the global mesh.
  /// \b Note: The nodes boundary information refers to the
  /// boundary numbers within the submesh!
  /// N.B. This is essentially the same function as the Mesh constructor
  /// that assembles a single global mesh from submeshes
  //=====================================================================
  void Problem::rebuild_global_mesh()
  {
    // Use the function in mesh to merge the submeshes into this one
    Mesh_pt->merge_meshes(Sub_mesh_pt);
  }


  //================================================================
  ///  Add a timestepper to the problem. The function will automatically
  /// create or resize the Time object so that it contains the appropriate
  /// number of levels of storage.
  //================================================================
  void Problem::add_time_stepper_pt(TimeStepper* const& time_stepper_pt)
  {
    // Add the timestepper to the vector
    Time_stepper_pt.push_back(time_stepper_pt);

    // Find the number of timesteps required by the timestepper
    unsigned ndt = time_stepper_pt->ndt();

    // If time has not been allocated, create time object with the
    // required number of time steps
    if (Time_pt == 0)
    {
      Time_pt = new Time(ndt);
      oomph_info << "Created Time with " << ndt << " timesteps" << std::endl;
    }
    else
    {
      // If the required number of time steps is greater than currently stored
      // resize the time storage
      if (ndt > Time_pt->ndt())
      {
        Time_pt->resize(ndt);
        oomph_info << "Resized Time to include " << ndt << " timesteps"
                   << std::endl;
      }
      // Otherwise report that we are OK
      else
      {
        oomph_info << "Time object already has storage for " << ndt
                   << " timesteps" << std::endl;
      }
    }

    // Pass the pointer to time to the timestepper
    time_stepper_pt->time_pt() = Time_pt;
  }

  //================================================================
  /// Set the explicit time stepper for the problem and also
  /// ensure that a time object has been created.
  //================================================================
  void Problem::set_explicit_time_stepper_pt(
    ExplicitTimeStepper* const& explicit_time_stepper_pt)
  {
    // Set the explicit time stepper
    Explicit_time_stepper_pt = explicit_time_stepper_pt;

    // If time has not been allocated, create time object with the
    // required number of time steps
    if (Time_pt == 0)
    {
      Time_pt = new Time(0);
      oomph_info << "Created Time with storage for no previous timestep"
                 << std::endl;
    }
    else
    {
      oomph_info << "Time object already exists " << std::endl;
    }
  }


#ifdef OOMPH_HAS_MPI

  //================================================================
  /// Set default first and last elements for parallel assembly
  /// of non-distributed problem.
  //================================================================
  void Problem::set_default_first_and_last_element_for_assembly()
  {
    if (!Problem_has_been_distributed)
    {
      // Minimum number of elements per processor if there are fewer elements
      // than processors
      unsigned min_el = 10;

      // Resize and make default assignments
      int n_proc = this->communicator_pt()->nproc();
      unsigned n_elements = Mesh_pt->nelement();
      First_el_for_assembly.resize(n_proc, 0);
      Last_el_plus_one_for_assembly.resize(n_proc, 0);

      // In the absence of any better knowledge distribute work evenly
      // over elements
      unsigned range = 0;
      unsigned lo_proc = 0;
      unsigned hi_proc = n_proc - 1;
      if (int(n_elements) >= n_proc)
      {
        range = unsigned(double(n_elements) / double(n_proc));
      }
      else
      {
        range = min_el;
        lo_proc = 0;
        hi_proc = unsigned(double(n_elements) / double(min_el));
      }

      for (int p = lo_proc; p <= int(hi_proc); p++)
      {
        First_el_for_assembly[p] = p * range;

        unsigned last_el_plus_one = (p + 1) * range;
        if (last_el_plus_one > n_elements) last_el_plus_one = n_elements;
        Last_el_plus_one_for_assembly[p] = last_el_plus_one;
      }

      // Last one needs to incorporate any dangling elements
      if (int(n_elements) >= n_proc)
      {
        Last_el_plus_one_for_assembly[n_proc - 1] = n_elements;
      }

      // Doc
      if (n_proc > 1)
      {
        if (!Shut_up_in_newton_solve)
        {
          oomph_info << "Problem is not distributed. Parallel assembly of "
                     << "Jacobian uses default partitioning: " << std::endl;
          for (int p = 0; p < n_proc; p++)
          {
            if (Last_el_plus_one_for_assembly[p] != 0)
            {
              oomph_info << "Proc " << p << " assembles from element "
                         << First_el_for_assembly[p] << " to "
                         << Last_el_plus_one_for_assembly[p] - 1 << " \n";
            }
            else
            {
              oomph_info << "Proc " << p << " assembles no elements\n";
            }
          }
        }
      }
    }
  }


  //=======================================================================
  /// Helper function to re-assign the first and last elements to be
  /// assembled by each processor during parallel assembly for
  /// non-distributed problem.
  //=======================================================================
  void Problem::recompute_load_balanced_assembly()
  {
    // Wait until all processes have completed/timed their assembly
    MPI_Barrier(this->communicator_pt()->mpi_comm());

    // Storage for number of processors and current processor
    int n_proc = this->communicator_pt()->nproc();
    int rank = this->communicator_pt()->my_rank();

    // Don't bother to update if we've got fewer elements than
    // processors
    unsigned nel = Elemental_assembly_time.size();
    if (int(nel) < n_proc)
    {
      oomph_info << "Not re-computing distribution of elemental assembly\n"
                 << "because there are fewer elements than processors\n";
      return;
    }

    // Setup vectors storing the number of element timings to be sent
    // and the offset in the final vector
    Vector<int> receive_count(n_proc);
    Vector<int> displacement(n_proc);
    int offset = 0;
    for (int p = 0; p < n_proc; p++)
    {
      // Default distribution of labour
      unsigned el_lo = First_el_for_assembly[p];
      unsigned el_hi = Last_el_plus_one_for_assembly[p] - 1;

      // Number of timings to be sent and offset from start in
      // final vector
      receive_count[p] = el_hi - el_lo + 1;
      displacement[p] = offset;
      offset += el_hi - el_lo + 1;
    }

    // Make temporary c-style array to avoid over-writing in Gatherv below
    double* el_ass_time = new double[nel];
    for (unsigned e = 0; e < nel; e++)
    {
      el_ass_time[e] = Elemental_assembly_time[e];
    }

    // Gather timings on root processor
    unsigned nel_local =
      Last_el_plus_one_for_assembly[rank] - 1 - First_el_for_assembly[rank] + 1;
    MPI_Gatherv(&el_ass_time[First_el_for_assembly[rank]],
                nel_local,
                MPI_DOUBLE,
                &Elemental_assembly_time[0],
                &receive_count[0],
                &displacement[0],
                MPI_DOUBLE,
                0,
                this->communicator_pt()->mpi_comm());
    delete[] el_ass_time;

    // Vector of first and last elements for each processor
    Vector<Vector<int>> first_and_last_element(n_proc);
    for (int p = 0; p < n_proc; p++)
    {
      first_and_last_element[p].resize(2);
    }

    // Re-distribute work
    if (rank == 0)
    {
      if (!Shut_up_in_newton_solve)
      {
        oomph_info
          << std::endl
          << "Re-assigning distribution of element assembly over processors:"
          << std::endl;
      }

      // Get total assembly time
      double total = 0.0;
      unsigned n_elements = Mesh_pt->nelement();
      for (unsigned e = 0; e < n_elements; e++)
      {
        total += Elemental_assembly_time[e];
      }

      // Target load per processor
      double target_load = total / double(n_proc);

      // We're on the root processor: Always start with the first element
      int proc = 0;
      first_and_last_element[0][0] = 0;

      // Highest element we can help ourselves to if we want to leave
      // at least one element for all subsequent processors
      unsigned max_el_avail = n_elements - n_proc;

      // Initialise total work allocated
      total = 0.0;
      for (unsigned e = 0; e < n_elements; e++)
      {
        total += Elemental_assembly_time[e];

        // Once we have reached the target load or we've used up practically
        // all the elements...
        if ((total > target_load) || (e == max_el_avail))
        {
          // Last element for current processor
          first_and_last_element[proc][1] = e;

          // Provided that we are not on the last processor
          if (proc < (n_proc - 1))
          {
            // Set first element for next one
            first_and_last_element[proc + 1][0] = e + 1;

            // Move on to the next processor
            proc++;
          }

          // Can have one more...
          max_el_avail++;

          // Re-initialise the time
          total = 0.0;
        } // end of test for "total exceeds target"
      }


      // Last element for last processor
      first_and_last_element[n_proc - 1][1] = n_elements - 1;


      // The following block should probably be paranoidified away
      // but we've screwed the logic up so many times that I feel
      // it's safer to keep it...
      bool wrong = false;
      std::ostringstream error_stream;
      for (int p = 0; p < n_proc - 1; p++)
      {
        unsigned first_of_current = first_and_last_element[p][0];
        unsigned last_of_current = first_and_last_element[p][1];
        if (first_of_current > last_of_current)
        {
          wrong = true;
          error_stream << "Error: First/last element of proc " << p << ": "
                       << first_of_current << " " << last_of_current
                       << std::endl;
        }
        unsigned first_of_next = first_and_last_element[p + 1][0];
        if (first_of_next != (last_of_current + 1))
        {
          wrong = true;
          error_stream << "Error: First element of proc " << p + 1 << ": "
                       << first_of_next << " and last element of proc " << p
                       << ": " << last_of_current << std::endl;
        }
      }
      if (wrong)
      {
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }


      // THIS TIDY UP SHOULD NO LONGER BE REQUIRED AND CAN GO AT SOME POINT

      // //If we haven't got to the end of the processor list then
      // //need to shift things about slightly because the processors
      // //at the end will be empty.
      // //This can occur when you have very fast assembly times and the
      // //rounding errors mean that the targets are achieved before all
      // processors
      // //have been visited.
      // //Happens a lot when you massively oversubscribe the CPUs (which was
      // //only ever for testing!)
      // if (proc!=n_proc-1)
      //  {
      //   oomph_info
      //    << "First pass did not allocate elements on every processor\n";
      //   oomph_info <<
      //    "Moving elements so that each processor has at least one\n";

      //   //Work out number of empty processos
      //   unsigned n_empty_processors = n_proc - proc + 1;

      //   //Loop over the processors that do have elements
      //   //and work out how many we need to steal elements from
      //   unsigned n_element_on_processors=0;
      //   do
      //    {
      //     //Step down the processors
      //     --proc;
      //     //Add the current processor to the number of empty processors
      //     //because the elements have to be shared between processors
      //     //including the one(s) on which they are currently stored.
      //     ++n_empty_processors;
      //     n_element_on_processors +=
      //      (first_and_last_element[proc][1] -
      //       first_and_last_element[proc][0] + 1);
      //    }
      //   while(n_element_on_processors < n_empty_processors);

      //   //Should now be able to put one element on each processor
      //   //Start from the end and do so
      //   unsigned current_element  = n_elements-1;
      //   for(int p=n_proc-1;p>proc;p--)
      //    {
      //     first_and_last_element[p][1] = current_element;
      //     first_and_last_element[p][0] = --current_element;
      //    }

      //   //Now for the last processor we touched, just adjust the final value
      //   first_and_last_element[proc][1] = current_element;
      //  }
      // //Otherwise just put the rest of the elements on the final
      // //processor
      // else
      //  {
      //   // Last one
      //   first_and_last_element[n_proc-1][1]=n_elements-1;
      //  }


      // END PRESUMED-TO-BE-UNNECESSARY BLOCK...


      // Now communicate the information

      // Set local informationt for this (root) processor
      First_el_for_assembly[0] = first_and_last_element[0][0];
      Last_el_plus_one_for_assembly[0] = first_and_last_element[0][1] + 1;

      if (!Shut_up_in_newton_solve)
      {
        oomph_info << "Processor " << 0 << " assembles Jacobians"
                   << " from elements " << first_and_last_element[0][0]
                   << " to " << first_and_last_element[0][1] << " "
                   << std::endl;
      }

      // Only now can we send the information to the other processors
      for (int p = 1; p < n_proc; ++p)
      {
        MPI_Send(&first_and_last_element[p][0],
                 2,
                 MPI_INT,
                 p,
                 0,
                 this->communicator_pt()->mpi_comm());


        if (!Shut_up_in_newton_solve)
        {
          oomph_info << "Processor " << p << " assembles Jacobians"
                     << " from elements " << first_and_last_element[p][0]
                     << " to " << first_and_last_element[p][1] << " "
                     << std::endl;
        }
      }
    }
    // Receive first and last element from root on non-master processors
    else
    {
      Vector<int> aux(2);
      MPI_Status status;
      MPI_Recv(&aux[0],
               2,
               MPI_INT,
               0,
               0,
               this->communicator_pt()->mpi_comm(),
               &status);
      First_el_for_assembly[rank] = aux[0];
      Last_el_plus_one_for_assembly[rank] = aux[1] + 1;
    }

    // Wipe all others
    for (int p = 0; p < n_proc; p++)
    {
      if (p != rank)
      {
        First_el_for_assembly[p] = 0;
        Last_el_plus_one_for_assembly[p] = 1;
      }
    }

    // The equations assembled by this processor may have changed so
    // we must resize the sparse assemble with arrays previous allocation
    Sparse_assemble_with_arrays_previous_allocation.resize(0);
  }

#endif

  //================================================================
  /// Assign all equation numbers for problem: Deals with global
  /// data (= data that isn't attached to any elements) and then
  /// does the equation numbering for the elements.  Bool argument
  /// can be set to false to ignore assigning local equation numbers
  /// (necessary in the parallel implementation of locate_zeta
  /// between multiple meshes).
  //================================================================
  unsigned long Problem::assign_eqn_numbers(
    const bool& assign_local_eqn_numbers)
  {
    // Check that the global mesh has been build
#ifdef PARANOID
    if (Mesh_pt == 0)
    {
      std::ostringstream error_stream;
      error_stream << "Global mesh does not exist, so equation numbers cannot "
                      "be assigned.\n";
      // Check for sub meshes
      if (nsub_mesh() == 0)
      {
        error_stream << "There aren't even any sub-meshes in the Problem.\n"
                     << "You can set the global mesh directly by using\n"
                     << "Problem::mesh_pt() = my_mesh_pt;\n"
                     << "OR you can use Problem::add_sub_mesh(mesh_pt); "
                     << "to add a sub mesh.\n";
      }
      else
      {
        error_stream << "There are " << nsub_mesh() << " sub-meshes.\n";
      }
      error_stream << "You need to call Problem::build_global_mesh() to create "
                      "a global mesh\n"
                   << "from the sub-meshes.\n\n";

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Number of submeshes
    unsigned n_sub_mesh = Sub_mesh_pt.size();

#ifdef OOMPH_HAS_MPI

    // Storage for number of processors
    int n_proc = this->communicator_pt()->nproc();


    if (n_proc > 1)
    {
      // Force re-analysis of time spent on assembly each
      // elemental Jacobian
      Must_recompute_load_balance_for_assembly = true;
      Elemental_assembly_time.clear();
    }
    else
    {
      Must_recompute_load_balance_for_assembly = false;
    }

    // Re-distribution of elements over processors during assembly
    // must be recomputed
    if (!Problem_has_been_distributed)
    {
      // Set default first and last elements for parallel assembly
      // of non-distributed problem.
      set_default_first_and_last_element_for_assembly();
    }

#endif


    double t_start = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start = TimingHelpers::timer();
    }

    // Loop over all elements in the mesh and set up any additional
    // dependencies that they may have (e.g. storing the geometric
    // Data, i.e. Data that affects an element's shape in elements
    // with algebraic node-update functions
    unsigned nel = Mesh_pt->nelement();
    for (unsigned e = 0; e < nel; e++)
    {
      Mesh_pt->element_pt(e)->complete_setup_of_dependencies();
    }

#ifdef OOMPH_HAS_MPI
    // Complete setup of dependencies for external halo elements too
    unsigned n_mesh = this->nsub_mesh();
    for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
    {
      for (int iproc = 0; iproc < n_proc; iproc++)
      {
        unsigned n_ext_halo_el = mesh_pt(i_mesh)->nexternal_halo_element(iproc);
        for (unsigned e = 0; e < n_ext_halo_el; e++)
        {
          mesh_pt(i_mesh)
            ->external_halo_element_pt(iproc, e)
            ->complete_setup_of_dependencies();
        }
      }
    }
#endif


    double t_end = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info
        << "Time for complete setup of dependencies in assign_eqn_numbers: "
        << t_end - t_start << std::endl;
    }


    // Initialise number of dofs for reserve below
    unsigned n_dof = 0;

    // Potentially loop over remainder of routine, possible re-visiting all
    // those parts that must be redone, following the removal of duplicate
    // external halo data.
    for (unsigned loop_count = 0; loop_count < 2; loop_count++)
    {
      //(Re)-set the dof pointer to zero length because entries are
      // pushed back onto it -- if it's not reset here then we get into
      // trouble during mesh refinement when we reassign all dofs
      Dof_pt.resize(0);

      // Reserve from previous allocation if we're going around again
      Dof_pt.reserve(n_dof);

      // Reset the equation number
      unsigned long equation_number = 0;

      // Now set equation numbers for the global Data
      unsigned Nglobal_data = nglobal_data();
      for (unsigned i = 0; i < Nglobal_data; i++)
      {
        Global_data_pt[i]->assign_eqn_numbers(equation_number, Dof_pt);
      }

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_start = TimingHelpers::timer();
      }

      // Call assign equation numbers on the global mesh
      n_dof = Mesh_pt->assign_global_eqn_numbers(Dof_pt);

      // Deal with the spine meshes additional numbering
      // If there is only one mesh
      if (n_sub_mesh == 0)
      {
        if (SpineMesh* const spine_mesh_pt = dynamic_cast<SpineMesh*>(Mesh_pt))
        {
          n_dof = spine_mesh_pt->assign_global_spine_eqn_numbers(Dof_pt);
        }
      }
      // Otherwise loop over the sub meshes
      else
      {
        // Assign global equation numbers first
        for (unsigned i = 0; i < n_sub_mesh; i++)
        {
          if (SpineMesh* const spine_mesh_pt =
                dynamic_cast<SpineMesh*>(Sub_mesh_pt[i]))
          {
            n_dof = spine_mesh_pt->assign_global_spine_eqn_numbers(Dof_pt);
          }
        }
      }

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info
          << "Time for assign_global_eqn_numbers in assign_eqn_numbers: "
          << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }


#ifdef OOMPH_HAS_MPI

      // reset previous allocation
      Parallel_sparse_assemble_previous_allocation = 0;

      // Only synchronise if the problem has actually been
      // distributed.
      if (Problem_has_been_distributed)
      {
        // Synchronise the equation numbers and return the total
        // number of degrees of freedom in the overall problem
        // Do not assign local equation numbers -- we're doing this
        // below.
        n_dof = synchronise_eqn_numbers(false);
      }
      // ..else just setup the Dof_distribution_pt
      // NOTE: this is setup by synchronise_eqn_numbers(...)
      // if Problem_has_been_distributed
      else
#endif
      {
        Dof_distribution_pt->build(Communicator_pt, n_dof, false);
      }

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for Problem::synchronise_eqn_numbers in "
                   << "Problem::assign_eqn_numbers: " << t_end - t_start
                   << std::endl;
      }


#ifdef OOMPH_HAS_MPI


      // Now remove duplicate data in external halo elements
      if (Problem_has_been_distributed)
      {
        if (Global_timings::Doc_comprehensive_timings)
        {
          t_start = TimingHelpers::timer();
        }

        // Monitor if we've actually changed anything
        bool actually_removed_some_data = false;

        // Only do it once!
        if (loop_count == 0)
        {
          if (n_sub_mesh == 0)
          {
            remove_duplicate_data(Mesh_pt, actually_removed_some_data);
          }
          else
          {
            for (unsigned i = 0; i < n_sub_mesh; i++)
            {
              bool tmp_actually_removed_some_data = false;
              remove_duplicate_data(Sub_mesh_pt[i],
                                    tmp_actually_removed_some_data);
              if (tmp_actually_removed_some_data)
                actually_removed_some_data = true;
            }
          }
        }


        if (Global_timings::Doc_comprehensive_timings)
        {
          t_end = TimingHelpers::timer();
          std::stringstream tmp;
          tmp << "Time for calls to Problem::remove_duplicate_data in "
              << "Problem::assign_eqn_numbers: " << t_end - t_start
              << " ; have ";
          if (!actually_removed_some_data)
          {
            tmp << " not ";
          }
          tmp << " removed some/any data.\n";
          oomph_info << tmp.str();
          t_start = TimingHelpers::timer();
        }

        // Break out of the loop if we haven't done anything here.
        unsigned status = 0;
        if (actually_removed_some_data) status = 1;

        // Allreduce to check if anyone has removed any data
        unsigned overall_status = 0;
        MPI_Allreduce(&status,
                      &overall_status,
                      1,
                      MPI_UNSIGNED,
                      MPI_MAX,
                      this->communicator_pt()->mpi_comm());


        if (Global_timings::Doc_comprehensive_timings)
        {
          t_end = TimingHelpers::timer();
          std::stringstream tmp;
          tmp
            << "Time for MPI_Allreduce after Problem::remove_duplicate_data in "
            << "Problem::assign_eqn_numbers: " << t_end - t_start << std::endl;
          oomph_info << tmp.str();
          t_start = TimingHelpers::timer();
        }

        // Bail out if we haven't done anything here
        if (overall_status != 1)
        {
          break;
        }

        // Big tidy up: Remove null pointers from halo/haloed node storage
        // for all meshes (this involves comms and therefore must be
        // performed outside loop over meshes so the all-to-all is only
        // done once)
        remove_null_pointers_from_external_halo_node_storage();

        // Time it...
        if (Global_timings::Doc_comprehensive_timings)
        {
          double t_end = TimingHelpers::timer();
          oomph_info << "Total time for "
                     << "Problem::remove_null_pointers_from_external_halo_node_"
                        "storage(): "
                     << t_end - t_start << std::endl;
        }
      }
      else
      {
        // Problem not distributed; no need for another loop
        break;
      }

#else

      // Serial run: Again no need for a second loop
      break;

#endif

    } // end of loop over fcts that need to be re-executed if
    // we've removed duplicate external data


    // Resize the sparse assemble with arrays previous allocation
    Sparse_assemble_with_arrays_previous_allocation.resize(0);


    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start = TimingHelpers::timer();
    }

    // Finally assign local equations
    if (assign_local_eqn_numbers)
    {
      if (n_sub_mesh == 0)
      {
        Mesh_pt->assign_local_eqn_numbers(Store_local_dof_pt_in_elements);
      }
      else
      {
        for (unsigned i = 0; i < n_sub_mesh; i++)
        {
          Sub_mesh_pt[i]->assign_local_eqn_numbers(
            Store_local_dof_pt_in_elements);
        }
      }
    }

    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Total time for all Mesh::assign_local_eqn_numbers in "
                 << "Problem::assign_eqn_numbers: " << t_end - t_start
                 << std::endl;
    }


    // and return the total number of DOFs
    return n_dof;
  }
  //================================================================
  /// \short Function to describe the dofs in terms of the global
  /// equation number, i.e. what type of value (nodal value of
  /// a Node; value in a Data object; value of internal Data in an
  /// element; etc) is the unknown with a certain global equation number.
  /// Output stream defaults to oomph_info.
  //================================================================
  void Problem::describe_dofs(std::ostream& out) const
  {
    // Check that the global mesh has been build
#ifdef PARANOID
    if (Mesh_pt == 0)
    {
      std::ostringstream error_stream;
      error_stream
        << "Global mesh does not exist, so equation numbers cannot be found.\n";
      // Check for sub meshes
      if (nsub_mesh() == 0)
      {
        error_stream << "There aren't even any sub-meshes in the Problem.\n"
                     << "You can set the global mesh directly by using\n"
                     << "Problem::mesh_pt() = my_mesh_pt;\n"
                     << "OR you can use Problem::add_sub_mesh(mesh_pt); "
                     << "to add a sub mesh.\n";
      }
      else
      {
        error_stream << "There are " << nsub_mesh() << " sub-meshes.\n";
      }
      error_stream << "You need to call Problem::build_global_mesh() to create "
                      "a global mesh\n"
                   << "from the sub-meshes.\n\n";

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    out
      << "Although this program will describe the degrees of freedom in the \n"
      << "problem, it will do so using the typedef for the elements. This is \n"
      << "not neccesarily human readable, but there is a solution.\n"
      << "Pipe your program's output through c++filt, with the argument -t.\n"
      << "e.g. \"./two_d_multi_poisson | c++filt -t > ReadableOutput.txt\".\n "
      << "(Disregarding the quotes)\n\n\n";

    out << "Classifying Global Equation Numbers" << std::endl;
    out << std::string(80, '-') << std::endl;

    // Number of submeshes
    unsigned n_sub_mesh = Sub_mesh_pt.size();

    // Classify Global dofs
    unsigned Nglobal_data = nglobal_data();
    for (unsigned i = 0; i < Nglobal_data; i++)
    {
      std::stringstream conversion;
      conversion << " in Global Data " << i << ".";
      std::string in(conversion.str());
      Global_data_pt[i]->describe_dofs(out, in);
    }

    // Put string in limiting scope.
    {
      // Descend into assignment for mesh.
      std::string in(" in Problem's Only Mesh.");
      Mesh_pt->describe_dofs(out, in);
    }

    // Deal with the spine meshes additional numbering:
    // If there is only one mesh:
    if (n_sub_mesh == 0)
    {
      if (SpineMesh* const spine_mesh_pt = dynamic_cast<SpineMesh*>(Mesh_pt))
      {
        std::string in(" in Problem's Only SpineMesh.");
        spine_mesh_pt->describe_spine_dofs(out, in);
      }
    }
    // Otherwise loop over the sub meshes
    else
    {
      // Assign global equation numbers first
      for (unsigned i = 0; i < n_sub_mesh; i++)
      {
        if (SpineMesh* const spine_mesh_pt =
              dynamic_cast<SpineMesh*>(Sub_mesh_pt[i]))
        {
          std::stringstream conversion;
          conversion << " in Sub-SpineMesh " << i << ".";
          std::string in(conversion.str());
          spine_mesh_pt->describe_spine_dofs(out, in);
        } // end if.
      } // end for.
    } // end else.


    out << std::string(80, '\\') << std::endl;
    out << std::string(80, '\\') << std::endl;
    out << std::string(80, '\\') << std::endl;
    out << "Classifying global eqn numbers in terms of elements." << std::endl;
    out << std::string(80, '-') << std::endl;
    out << "Eqns   | Source" << std::endl;
    out << std::string(80, '-') << std::endl;

    if (n_sub_mesh == 0)
    {
      std::string in(" in Problem's Only Mesh.");
      Mesh_pt->describe_local_dofs(out, in);
    }
    else
    {
      for (unsigned i = 0; i < n_sub_mesh; i++)
      {
        std::stringstream conversion;
        conversion << " in Sub-Mesh " << i << ".";
        std::string in(conversion.str());
        Sub_mesh_pt[i]->describe_local_dofs(out, in);
      } // End for
    } // End else
  } // End problem::describe_dofs(...)


  //================================================================
  /// Get the vector of dofs, i.e. a vector containing the current
  /// values of all unknowns.
  //================================================================
  void Problem::get_dofs(DoubleVector& dofs) const
  {
    // Find number of dofs
    const unsigned long n_dof = ndof();

    // Resize the vector
    dofs.build(Dof_distribution_pt, 0.0);

    // Copy dofs into vector
    for (unsigned long l = 0; l < n_dof; l++)
    {
      dofs[l] = *Dof_pt[l];
    }
  }

  /// Get history values of dofs in a double vector.
  void Problem::get_dofs(const unsigned& t, DoubleVector& dofs) const
  {
#ifdef PARANOID
    if (distributed())
    {
      throw OomphLibError("Not designed for distributed problems",
                          OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
      // might work, not sure
    }
#endif

    // Resize the vector
    dofs.build(Dof_distribution_pt, 0.0);

    // First deal with global data
    unsigned Nglobal_data = nglobal_data();
    for (unsigned i = 0; i < Nglobal_data; i++)
    {
      for (unsigned j = 0, nj = Global_data_pt[i]->nvalue(); j < nj; j++)
      {
        // For each data get the equation number and copy out the value.
        int eqn_number = Global_data_pt[i]->eqn_number(j);
        if (eqn_number >= 0)
        {
          dofs[eqn_number] = Global_data_pt[i]->value(t, j);
        }
      }
    }

    // Next element internal data
    for (unsigned i = 0, ni = mesh_pt()->nelement(); i < ni; i++)
    {
      GeneralisedElement* ele_pt = mesh_pt()->element_pt(i);
      for (unsigned j = 0, nj = ele_pt->ninternal_data(); j < nj; j++)
      {
        Data* d_pt = ele_pt->internal_data_pt(j);
        for (unsigned k = 0, nk = d_pt->nvalue(); k < nk; k++)
        {
          int eqn_number = d_pt->eqn_number(k);
          if (eqn_number >= 0)
          {
            dofs[eqn_number] = d_pt->value(t, k);
          }
        }
      }
    }

    // Now the nodes
    for (unsigned i = 0, ni = mesh_pt()->nnode(); i < ni; i++)
    {
      Node* node_pt = mesh_pt()->node_pt(i);
      for (unsigned j = 0, nj = node_pt->nvalue(); j < nj; j++)
      {
        // For each node get the equation number and copy out the value.
        int eqn_number = node_pt->eqn_number(j);
        if (eqn_number >= 0)
        {
          dofs[eqn_number] = node_pt->value(t, j);
        }
      }
    }
  }


#ifdef OOMPH_HAS_MPI

  //=======================================================================
  /// Private helper function to remove repeated data
  /// in external haloed elements associated with specified mesh.
  /// Bool is true if some data was removed -- this usually requires
  /// re-running through certain parts of the equation numbering procedure.
  //======================================================================
  void Problem::remove_duplicate_data(Mesh* const& mesh_pt,
                                      bool& actually_removed_some_data)
  {
    //    //   // Taken out again by MH -- clutters up output
    //    // Doc timings if required
    //    double t_start=0.0;
    //    if (Global_timings::Doc_comprehensive_timings)
    //     {
    //      t_start=TimingHelpers::timer();
    //     }

    int n_proc = this->communicator_pt()->nproc();
    int my_rank = this->communicator_pt()->my_rank();

    // Initialise
    actually_removed_some_data = false;

    // Each individual container of external halo nodes has unique
    // nodes/equation numbers, but there may be some duplication between
    // two or more different containers; the following code checks for this
    // and removes the duplication by overwriting any data point with an already
    // existing eqn number with the original data point which had the eqn no.

    // // Storage for existing nodes, enumerated by first non-negative
    // // global equation number
    // unsigned n_dof=ndof();

    // Note: This used to be
    // Vector<Node*> global_node_pt(n_dof,0);
    // but this is a total killer! Memory allocation is extremely
    // costly and only relatively few entries are used so use
    // map:
    std::map<unsigned, Node*> global_node_pt;

    // Only do each retained node once
    std::map<Node*, bool> node_done;

    // Loop over existing "normal" elements in mesh
    unsigned n_element = mesh_pt->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      FiniteElement* el_pt =
        dynamic_cast<FiniteElement*>(mesh_pt->element_pt(e));
      if (el_pt != 0)
      {
        // Loop over nodes
        unsigned n_node = el_pt->nnode();
        for (unsigned j = 0; j < n_node; j++)
        {
          Node* nod_pt = el_pt->node_pt(j);

          // Have we already done the node?
          if (!node_done[nod_pt])
          {
            node_done[nod_pt] = true;

            // Loop over values stored at node (if any) to find
            // the first non-negative eqn number
            unsigned first_non_negative_eqn_number_plus_one = 0;
            unsigned n_val = nod_pt->nvalue();
            for (unsigned i_val = 0; i_val < n_val; i_val++)
            {
              int eqn_no = nod_pt->eqn_number(i_val);
              if (eqn_no >= 0)
              {
                first_non_negative_eqn_number_plus_one = eqn_no + 1;
                break;
              }
            }

            // If we haven't found a non-negative eqn number check
            // eqn numbers associated with solid data (if any)
            if (first_non_negative_eqn_number_plus_one == 0)
            {
              // Is it a solid node?
              SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
              if (solid_nod_pt != 0)
              {
                // Loop over values stored at node (if any) to find
                // the first non-negative eqn number
                unsigned n_val = solid_nod_pt->variable_position_pt()->nvalue();
                for (unsigned i_val = 0; i_val < n_val; i_val++)
                {
                  int eqn_no =
                    solid_nod_pt->variable_position_pt()->eqn_number(i_val);
                  if (eqn_no >= 0)
                  {
                    first_non_negative_eqn_number_plus_one = eqn_no + 1;
                    break;
                  }
                }
              }
            }

            // Associate node with first non negative global eqn number
            if (first_non_negative_eqn_number_plus_one > 0)
            {
              global_node_pt[first_non_negative_eqn_number_plus_one - 1] =
                nod_pt;
            }


            // Take into account master nodes too
            if (dynamic_cast<RefineableElement*>(el_pt) != 0)
            {
              int n_cont_int_values = dynamic_cast<RefineableElement*>(el_pt)
                                        ->ncont_interpolated_values();
              for (int i_cont = -1; i_cont < n_cont_int_values; i_cont++)
              {
                if (nod_pt->is_hanging(i_cont))
                {
                  HangInfo* hang_pt = nod_pt->hanging_pt(i_cont);
                  unsigned n_master = hang_pt->nmaster();
                  for (unsigned m = 0; m < n_master; m++)
                  {
                    Node* master_nod_pt = hang_pt->master_node_pt(m);
                    if (!node_done[master_nod_pt])
                    {
                      node_done[master_nod_pt] = true;

                      // Loop over values stored at node (if any) to find
                      // the first non-negative eqn number
                      unsigned first_non_negative_eqn_number_plus_one = 0;
                      unsigned n_val = master_nod_pt->nvalue();
                      for (unsigned i_val = 0; i_val < n_val; i_val++)
                      {
                        int eqn_no = master_nod_pt->eqn_number(i_val);
                        if (eqn_no >= 0)
                        {
                          first_non_negative_eqn_number_plus_one = eqn_no + 1;
                          break;
                        }
                      }

                      // If we haven't found a non-negative eqn number check
                      // eqn numbers associated with solid data (if any)
                      if (first_non_negative_eqn_number_plus_one == 0)
                      {
                        // If this master is a SolidNode then add its extra
                        // eqn numbers
                        SolidNode* master_solid_nod_pt =
                          dynamic_cast<SolidNode*>(master_nod_pt);
                        if (master_solid_nod_pt != 0)
                        {
                          // Loop over values stored at node (if any) to find
                          // the first non-negative eqn number
                          unsigned n_val =
                            master_solid_nod_pt->variable_position_pt()
                              ->nvalue();
                          for (unsigned i_val = 0; i_val < n_val; i_val++)
                          {
                            int eqn_no =
                              master_solid_nod_pt->variable_position_pt()
                                ->eqn_number(i_val);
                            if (eqn_no >= 0)
                            {
                              first_non_negative_eqn_number_plus_one =
                                eqn_no + 1;
                              break;
                            }
                          }
                        }
                      }
                      // Associate node with first non negative global
                      // eqn number
                      if (first_non_negative_eqn_number_plus_one > 0)
                      {
                        global_node_pt[first_non_negative_eqn_number_plus_one -
                                       1] = master_nod_pt;
                      }

                    } // End of not-yet-done hang node
                  }
                }
              }
            }
          } // endif for node already done
        } // End of loop over nodes
      } // End of FiniteElement

      // Internal data equation numbers do not need to be added since
      // internal data cannot be shared between distinct elements, so
      // internal data on locally-stored elements can never be halo.
    }

    // Set to record duplicate nodes scheduled to be killed
    std::set<Node*> killed_nodes;

    // Now loop over the other processors from highest to lowest
    // (i.e. if there is a duplicate between these containers
    //  then this will use the node on the highest numbered processor)
    for (int iproc = n_proc - 1; iproc >= 0; iproc--)
    {
      // Don't have external halo elements with yourself!
      if (iproc != my_rank)
      {
        // Loop over external halo elements with iproc
        // to remove the duplicates
        unsigned n_element = mesh_pt->nexternal_halo_element(iproc);
        for (unsigned e_ext = 0; e_ext < n_element; e_ext++)
        {
          FiniteElement* finite_ext_el_pt = dynamic_cast<FiniteElement*>(
            mesh_pt->external_halo_element_pt(iproc, e_ext));
          if (finite_ext_el_pt != 0)
          {
            // Loop over nodes
            unsigned n_node = finite_ext_el_pt->nnode();
            for (unsigned j = 0; j < n_node; j++)
            {
              Node* nod_pt = finite_ext_el_pt->node_pt(j);

              // Loop over values stored at node (if any) to find
              // the first non-negative eqn number
              unsigned first_non_negative_eqn_number_plus_one = 0;
              unsigned n_val = nod_pt->nvalue();
              for (unsigned i_val = 0; i_val < n_val; i_val++)
              {
                int eqn_no = nod_pt->eqn_number(i_val);
                if (eqn_no >= 0)
                {
                  first_non_negative_eqn_number_plus_one = eqn_no + 1;
                  break;
                }
              }

              // If we haven't found a non-negative eqn number check
              // eqn numbers associated with solid data (if any)
              if (first_non_negative_eqn_number_plus_one == 0)
              {
                // Is it a solid node?
                SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
                if (solid_nod_pt != 0)
                {
                  // Loop over values stored at node (if any) to find
                  // the first non-negative eqn number
                  unsigned n_val =
                    solid_nod_pt->variable_position_pt()->nvalue();
                  for (unsigned i_val = 0; i_val < n_val; i_val++)
                  {
                    int eqn_no =
                      solid_nod_pt->variable_position_pt()->eqn_number(i_val);
                    if (eqn_no >= 0)
                    {
                      first_non_negative_eqn_number_plus_one = eqn_no + 1;
                      break;
                    }
                  }
                }
              }

              // Identified which node we're dealing with via first non-negative
              // global eqn number (if there is none, everything is pinned
              // and we don't give a damn...)
              if (first_non_negative_eqn_number_plus_one > 0)
              {
                Node* existing_node_pt =
                  global_node_pt[first_non_negative_eqn_number_plus_one - 1];

                // Does this node already exist?
                if (existing_node_pt != 0)
                {
                  // Record that we're about to cull one
                  actually_removed_some_data = true;

                  // It's a duplicate, so store the duplicated one for
                  // later killing...
                  Node* duplicated_node_pt = nod_pt;
                  if (!node_done[duplicated_node_pt])
                  {
                    // Remove node from all boundaries
                    std::set<unsigned>* boundaries_pt;
                    duplicated_node_pt->get_boundaries_pt(boundaries_pt);
                    if (boundaries_pt != 0)
                    {
                      Vector<unsigned> bound;
                      unsigned nb = (*boundaries_pt).size();
                      bound.reserve(nb);
                      for (std::set<unsigned>::iterator it =
                             (*boundaries_pt).begin();
                           it != (*boundaries_pt).end();
                           it++)
                      {
                        bound.push_back((*it));
                      }
                      for (unsigned i = 0; i < nb; i++)
                      {
                        mesh_pt->remove_boundary_node(bound[i],
                                                      duplicated_node_pt);
                      }
                    }

                    // Get ready to kill it
                    killed_nodes.insert(duplicated_node_pt);
                    unsigned i_proc = unsigned(iproc);
                    mesh_pt->null_external_halo_node(i_proc,
                                                     duplicated_node_pt);
                  }


                  // Note: For now we're leaving the "dangling" (no longer
                  // accessed masters where they are; they get cleaned
                  // up next time we delete all the external storage
                  // for the meshes so it's a temporary "leak" only...
                  // At some point we should probably delete them properly too

#ifdef PARANOID

                  // Check that hang status of exiting and replacement node
                  // matches
                  if (dynamic_cast<RefineableElement*>(finite_ext_el_pt) != 0)
                  {
                    int n_cont_inter_values =
                      dynamic_cast<RefineableElement*>(finite_ext_el_pt)
                        ->ncont_interpolated_values();
                    for (int i_cont = -1; i_cont < n_cont_inter_values;
                         i_cont++)
                    {
                      unsigned n_master_orig = 0;
                      if (finite_ext_el_pt->node_pt(j)->is_hanging(i_cont))
                      {
                        n_master_orig = finite_ext_el_pt->node_pt(j)
                                          ->hanging_pt(i_cont)
                                          ->nmaster();

                        // Temporary leak: Resolve like this:
                        // loop over all external halo nodes and identify the
                        // the ones that are still reached by any of the
                        // external elements. Kill the dangling ones.
                      }
                      unsigned n_master_replace = 0;
                      if (existing_node_pt->is_hanging(i_cont))
                      {
                        n_master_replace =
                          existing_node_pt->hanging_pt(i_cont)->nmaster();
                      }

                      if (n_master_orig != n_master_replace)
                      {
                        std::ostringstream error_stream;
                        error_stream
                          << "Number of master nodes for node to be replaced, "
                          << n_master_orig << ", doesn't match"
                          << "those of replacement node, " << n_master_replace
                          << " for i_cont=" << i_cont << std::endl;
                        {
                          error_stream
                            << "Nodal coordinates of replacement node:";
                          unsigned ndim = existing_node_pt->ndim();
                          for (unsigned i = 0; i < ndim; i++)
                          {
                            error_stream << existing_node_pt->x(i) << " ";
                          }
                          error_stream << "\n";
                          error_stream << "The coordinates of its "
                                       << n_master_replace
                                       << " master nodes are: \n";
                          for (unsigned k = 0; k < n_master_replace; k++)
                          {
                            Node* master_nod_pt =
                              existing_node_pt->hanging_pt(i_cont)
                                ->master_node_pt(k);
                            unsigned ndim = master_nod_pt->ndim();
                            for (unsigned i = 0; i < ndim; i++)
                            {
                              error_stream << master_nod_pt->x(i) << " ";
                            }
                            error_stream << "\n";
                          }
                        }

                        {
                          error_stream
                            << "Nodal coordinates of node to be replaced:";
                          unsigned ndim = finite_ext_el_pt->node_pt(j)->ndim();
                          for (unsigned i = 0; i < ndim; i++)
                          {
                            error_stream << finite_ext_el_pt->node_pt(j)->x(i)
                                         << " ";
                          }
                          error_stream << "\n";
                          error_stream << "The coordinates of its "
                                       << n_master_orig
                                       << " master nodes are: \n";
                          for (unsigned k = 0; k < n_master_orig; k++)
                          {
                            Node* master_nod_pt = finite_ext_el_pt->node_pt(j)
                                                    ->hanging_pt(i_cont)
                                                    ->master_node_pt(k);
                            unsigned ndim = master_nod_pt->ndim();
                            for (unsigned i = 0; i < ndim; i++)
                            {
                              error_stream << master_nod_pt->x(i) << " ";
                            }
                            error_stream << "\n";
                          }
                        }


                        throw OomphLibError(error_stream.str(),
                                            OOMPH_CURRENT_FUNCTION,
                                            OOMPH_EXCEPTION_LOCATION);
                      }
                    }
                  }
#endif
                  // ...and point to the existing one
                  finite_ext_el_pt->node_pt(j) = existing_node_pt;
                }
                // If it doesn't add it to the list of existing ones
                else
                {
                  global_node_pt[first_non_negative_eqn_number_plus_one - 1] =
                    nod_pt;
                  node_done[nod_pt] = true;
                }
              }


              // Do the same for any master nodes of that (possibly replaced)
              // node
              if (dynamic_cast<RefineableElement*>(finite_ext_el_pt) != 0)
              {
                int n_cont_inter_values =
                  dynamic_cast<RefineableElement*>(finite_ext_el_pt)
                    ->ncont_interpolated_values();
                for (int i_cont = -1; i_cont < n_cont_inter_values; i_cont++)
                {
                  if (finite_ext_el_pt->node_pt(j)->is_hanging(i_cont))
                  {
                    HangInfo* hang_pt =
                      finite_ext_el_pt->node_pt(j)->hanging_pt(i_cont);
                    unsigned n_master = hang_pt->nmaster();
                    for (unsigned m = 0; m < n_master; m++)
                    {
                      Node* master_nod_pt = hang_pt->master_node_pt(m);
                      unsigned n_val = master_nod_pt->nvalue();
                      unsigned first_non_negative_eqn_number_plus_one = 0;
                      for (unsigned i_val = 0; i_val < n_val; i_val++)
                      {
                        int eqn_no = master_nod_pt->eqn_number(i_val);
                        if (eqn_no >= 0)
                        {
                          first_non_negative_eqn_number_plus_one = eqn_no + 1;
                          break;
                        }
                      }

                      // If we haven't found a non-negative eqn number check
                      // eqn numbers associated with solid data (if any)
                      if (first_non_negative_eqn_number_plus_one == 0)
                      {
                        SolidNode* solid_master_nod_pt =
                          dynamic_cast<SolidNode*>(master_nod_pt);
                        if (solid_master_nod_pt != 0)
                        {
                          // Loop over values stored at node (if any) to find
                          // the first non-negative eqn number
                          unsigned n_val =
                            solid_master_nod_pt->variable_position_pt()
                              ->nvalue();
                          for (unsigned i_val = 0; i_val < n_val; i_val++)
                          {
                            int eqn_no =
                              solid_master_nod_pt->variable_position_pt()
                                ->eqn_number(i_val);
                            if (eqn_no >= 0)
                            {
                              first_non_negative_eqn_number_plus_one =
                                eqn_no + 1;
                              break;
                            }
                          }
                        }
                      }

                      // Identified which node we're dealing with via
                      // first non-negative global eqn number (if there
                      // is none, everything is pinned and we don't give a
                      // damn...)
                      if (first_non_negative_eqn_number_plus_one > 0)
                      {
                        Node* existing_node_pt = global_node_pt
                          [first_non_negative_eqn_number_plus_one - 1];

                        // Does this node already exist?
                        if (existing_node_pt != 0)
                        {
                          // Record that we're about to cull one
                          actually_removed_some_data = true;

                          // It's a duplicate, so store the duplicated one for
                          // later killing...
                          Node* duplicated_node_pt = master_nod_pt;

                          if (!node_done[duplicated_node_pt])
                          {
                            // Remove node from all boundaries
                            std::set<unsigned>* boundaries_pt;
                            duplicated_node_pt->get_boundaries_pt(
                              boundaries_pt);
                            if (boundaries_pt != 0)
                            {
                              for (std::set<unsigned>::iterator it =
                                     (*boundaries_pt).begin();
                                   it != (*boundaries_pt).end();
                                   it++)
                              {
                                mesh_pt->remove_boundary_node(
                                  (*it), duplicated_node_pt);
                              }
                            }

                            killed_nodes.insert(duplicated_node_pt);
                            unsigned i_proc = unsigned(iproc);
                            mesh_pt->null_external_halo_node(
                              i_proc, duplicated_node_pt);
                          }

                          // Weight of the original node
                          double m_weight = hang_pt->master_weight(m);


#ifdef PARANOID
                          // Sanity check: setting replacement master
                          // node for non-hanging node? Sign of really
                          // f***ed up code.
                          Node* tmp_nod_pt = finite_ext_el_pt->node_pt(j);
                          if (!tmp_nod_pt->is_hanging(i_cont))
                          {
                            std::ostringstream error_stream;
                            error_stream
                              << "About to re-set master for i_cont= " << i_cont
                              << " for external node (with proc " << iproc
                              << " )" << tmp_nod_pt << " at ";
                            unsigned n = tmp_nod_pt->ndim();
                            for (unsigned jj = 0; jj < n; jj++)
                            {
                              error_stream << tmp_nod_pt->x(jj) << " ";
                            }
                            error_stream
                              << " which is not hanging --> About to die!"
                              << "Outputting offending element into oomph-info "
                              << "stream. \n\n";
                            oomph_info << "\n\n";
                            finite_ext_el_pt->output(*(oomph_info.stream_pt()));
                            oomph_info << "\n\n";
                            oomph_info.stream_pt()->flush();
                            throw OomphLibError(error_stream.str(),
                                                OOMPH_CURRENT_FUNCTION,
                                                OOMPH_EXCEPTION_LOCATION);
                          }
#endif


                          // And re-set master
                          finite_ext_el_pt->node_pt(j)
                            ->hanging_pt(i_cont)
                            ->set_master_node_pt(m, existing_node_pt, m_weight);
                        }
                        // If it doesn't, add it to the list of existing ones
                        else
                        {
                          global_node_pt
                            [first_non_negative_eqn_number_plus_one - 1] =
                              master_nod_pt;
                          node_done[master_nod_pt] = true;
                        }
                      }
                    } // End of loop over master nodes
                  } // end of hanging
                } // end of loop over continously interpolated variables
              } // end refineable element (with potentially hanging node

            } // end loop over nodes on external halo elements

          } // End of check for finite element

        } // end loop over external halo elements
      }
    } // end loop over processors


    // Now kill all the deleted nodes
    for (std::set<Node*>::iterator it = killed_nodes.begin();
         it != killed_nodes.end();
         it++)
    {
      delete (*it);
    }


    //   oomph_info << "Number of nonzero entries in global_node_pt: "
    //              << global_node_pt.size() << std::endl;

    //    // Time it...
    //    // Taken out again by MH -- clutters up output
    //    if (Global_timings::Doc_comprehensive_timings)
    //     {
    //      double t_end = TimingHelpers::timer();
    //      oomph_info
    //       << "Total time for Problem::remove_duplicate_data: "
    //       << t_end-t_start << std::endl;
    //     }
  }


  //========================================================================
  /// Consolidate external halo node storage by removing nulled out
  /// pointers in external halo and haloed schemes for all meshes.
  //========================================================================
  void Problem::remove_null_pointers_from_external_halo_node_storage()
  {
    // Do we have submeshes?
    unsigned n_mesh_loop = 1;
    unsigned nmesh = nsub_mesh();
    if (nmesh > 0)
    {
      n_mesh_loop = nmesh;
    }

    // Storage for number of processors and current processor
    int n_proc = this->communicator_pt()->nproc();
    int my_rank = this->communicator_pt()->my_rank();

    // If only one processor then return
    if (n_proc == 1)
    {
      return;
    }

    // Loop over all (other) processors and store index of any nulled-out
    // external halo nodes in storage scheme.

    // Data to be sent to each processor
    Vector<int> send_n(n_proc, 0);

    // Storage for all values to be sent to all processors
    Vector<int> send_data;

    // Start location within send_data for data to be sent to each processor
    Vector<int> send_displacement(n_proc, 0);

    // Check missing ones
    for (int domain = 0; domain < n_proc; domain++)
    {
      // Set the offset for the current processor
      send_displacement[domain] = send_data.size();

      // Don't bother to do anything if the processor in the loop is the
      // current processor
      if (domain != my_rank)
      {
        // Deal with sub-meshes one-by-one if required
        Mesh* my_mesh_pt = 0;

        // Loop over submeshes
        for (unsigned imesh = 0; imesh < n_mesh_loop; imesh++)
        {
          if (nmesh == 0)
          {
            my_mesh_pt = mesh_pt();
          }
          else
          {
            my_mesh_pt = mesh_pt(imesh);
          }

          // Make backup of external halo node pointers with this domain
          Vector<Node*> backup_pt(my_mesh_pt->external_halo_node_pt(domain));

          // How many do we have currently?
          unsigned nnod = backup_pt.size();

          // Prepare storage for updated halo nodes
          Vector<Node*> new_external_halo_node_pt;
          new_external_halo_node_pt.reserve(nnod);

          // Loop over external halo nodes with this domain
          for (unsigned j = 0; j < nnod; j++)
          {
            // Get pointer to node
            Node* nod_pt = backup_pt[j];

            // Has it been nulled out?
            if (nod_pt == 0)
            {
              // Save index of nulled out one
              send_data.push_back(j);
            }
            else
            {
              // Still alive: Copy across
              new_external_halo_node_pt.push_back(nod_pt);
            }
          }

          // Set new external halo node vector
          my_mesh_pt->set_external_halo_node_pt(domain,
                                                new_external_halo_node_pt);

          // End of data for this mesh
          send_data.push_back(-1);

        } // end of loop over meshes

      } // end skip own domain

      // Find the number of data added to the vector
      send_n[domain] = send_data.size() - send_displacement[domain];
    }

    // Storage for the number of data to be received from each processor
    Vector<int> receive_n(n_proc, 0);

    // Now send numbers of data to be sent between all processors
    MPI_Alltoall(&send_n[0],
                 1,
                 MPI_INT,
                 &receive_n[0],
                 1,
                 MPI_INT,
                 this->communicator_pt()->mpi_comm());


    // We now prepare the data to be received
    // by working out the displacements from the received data
    Vector<int> receive_displacement(n_proc, 0);
    int receive_data_count = 0;
    for (int rank = 0; rank < n_proc; ++rank)
    {
      // Displacement is number of data received so far
      receive_displacement[rank] = receive_data_count;
      receive_data_count += receive_n[rank];
    }

    // Now resize the receive buffer for all data from all processors
    // Make sure that it has a size of at least one
    if (receive_data_count == 0)
    {
      ++receive_data_count;
    }
    Vector<int> receive_data(receive_data_count);

    // Make sure that the send buffer has size at least one
    // so that we don't get a segmentation fault
    if (send_data.size() == 0)
    {
      send_data.resize(1);
    }

    // Now send the data between all the processors
    MPI_Alltoallv(&send_data[0],
                  &send_n[0],
                  &send_displacement[0],
                  MPI_INT,
                  &receive_data[0],
                  &receive_n[0],
                  &receive_displacement[0],
                  MPI_INT,
                  this->communicator_pt()->mpi_comm());

    // Now use the received data
    for (int send_rank = 0; send_rank < n_proc; send_rank++)
    {
      // Don't bother to do anything for the processor corresponding to the
      // current processor or if no data were received from this processor
      if ((send_rank != my_rank) && (receive_n[send_rank] != 0))
      {
        // Counter for the data within the large array
        unsigned count = receive_displacement[send_rank];

        // Deal with sub-meshes one-by-one if required
        Mesh* my_mesh_pt = 0;

        // Loop over submeshes
        for (unsigned imesh = 0; imesh < n_mesh_loop; imesh++)
        {
          if (nmesh == 0)
          {
            my_mesh_pt = mesh_pt();
          }
          else
          {
            my_mesh_pt = mesh_pt(imesh);
          }

          // Make backup of external haloed node pointers with this domain
          Vector<Node*> backup_pt =
            my_mesh_pt->external_haloed_node_pt(send_rank);

          // Unpack until we reach "end of data" indicator (-1) for this mesh
          while (true)
          {
            // Read next entry
            int next_one = receive_data[count++];

            if (next_one == -1)
            {
              break;
            }
            else
            {
              // Null out the entry
              backup_pt[next_one] = 0;
            }
          }

          // How many do we have currently?
          unsigned nnod = backup_pt.size();

          // Prepare storage for updated haloed nodes
          Vector<Node*> new_external_haloed_node_pt;
          new_external_haloed_node_pt.reserve(nnod);

          // Loop over external haloed nodes with this domain
          for (unsigned j = 0; j < nnod; j++)
          {
            // Get pointer to node
            Node* nod_pt = backup_pt[j];

            // Has it been nulled out?
            if (nod_pt != 0)
            {
              // Still alive: Copy across
              new_external_haloed_node_pt.push_back(nod_pt);
            }
          }

          // Set new external haloed node vector
          my_mesh_pt->set_external_haloed_node_pt(send_rank,
                                                  new_external_haloed_node_pt);
        }
      }

    } // End of data is received
  }

#endif


  //=======================================================================
  /// Function that sets the values of the dofs in the object
  //======================================================================
  void Problem::set_dofs(const DoubleVector& dofs)
  {
    const unsigned long n_dof = this->ndof();
#ifdef PARANOID
    if (n_dof != dofs.nrow())
    {
      std::ostringstream error_stream;
      error_stream << "Number of degrees of freedom in vector argument "
                   << dofs.nrow() << "\n"
                   << "does not equal number of degrees of freedom in problem "
                   << n_dof;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif
    for (unsigned long l = 0; l < n_dof; l++)
    {
      *Dof_pt[l] = dofs[l];
    }
  }

  /// Set history values of dofs
  void Problem::set_dofs(const unsigned& t, DoubleVector& dofs)
  {
#ifdef PARANOID
    if (distributed())
    {
      throw OomphLibError("Not designed for distributed problems",
                          OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
      // might work if the dofs vector is distributed in the right way...
    }
#endif

    // First deal with global data
    unsigned Nglobal_data = nglobal_data();
    for (unsigned i = 0; i < Nglobal_data; i++)
    {
      for (unsigned j = 0, nj = Global_data_pt[i]->nvalue(); j < nj; j++)
      {
        // For each data get the equation number and copy out the value.
        int eqn_number = Global_data_pt[i]->eqn_number(j);
        if (eqn_number >= 0)
        {
          Global_data_pt[i]->set_value(t, j, dofs[eqn_number]);
        }
      }
    }

    // Next element internal data
    for (unsigned i = 0, ni = mesh_pt()->nelement(); i < ni; i++)
    {
      GeneralisedElement* ele_pt = mesh_pt()->element_pt(i);
      for (unsigned j = 0, nj = ele_pt->ninternal_data(); j < nj; j++)
      {
        Data* d_pt = ele_pt->internal_data_pt(j);
        for (unsigned k = 0, nk = d_pt->nvalue(); k < nk; k++)
        {
          int eqn_number = d_pt->eqn_number(k);
          if (eqn_number >= 0)
          {
            d_pt->set_value(t, k, dofs[eqn_number]);
          }
        }
      }
    }

    // Now the nodes
    for (unsigned i = 0, ni = mesh_pt()->nnode(); i < ni; i++)
    {
      Node* node_pt = mesh_pt()->node_pt(i);
      for (unsigned j = 0, nj = node_pt->nvalue(); j < nj; j++)
      {
        // For each node get the equation number and copy out the value.
        int eqn_number = node_pt->eqn_number(j);
        if (eqn_number >= 0)
        {
          node_pt->set_value(t, j, dofs[eqn_number]);
        }
      }
    }
  }


  /// Set history values of dofs from the type of vector stored in
  /// problem::Dof_pt.
  void Problem::set_dofs(const unsigned& t, Vector<double*>& dof_pt)
  {
#ifdef PARANOID
    if (distributed())
    {
      throw OomphLibError("Not implemented for distributed problems!",
                          OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
    }
#endif

    // If we have any spine meshes I think there might be more degrees
    // of freedom there. I don't use them though so I'll let someone who
    // knows what they are doing handle it. --David Shepherd

    // First deal with global data
    unsigned Nglobal_data = nglobal_data();
    for (unsigned i = 0; i < Nglobal_data; i++)
    {
      for (unsigned j = 0, nj = Global_data_pt[i]->nvalue(); j < nj; j++)
      {
        // For each data get the equation number and copy in the value.
        int eqn_number = Global_data_pt[i]->eqn_number(j);
        if (eqn_number >= 0)
        {
          Global_data_pt[i]->set_value(t, j, *(dof_pt[eqn_number]));
        }
      }
    }

    // Now the mesh data
    // nodes
    for (unsigned i = 0, ni = mesh_pt()->nnode(); i < ni; i++)
    {
      Node* node_pt = mesh_pt()->node_pt(i);
      for (unsigned j = 0, nj = node_pt->nvalue(); j < nj; j++)
      {
        // For each node get the equation number and copy in the value.
        int eqn_number = node_pt->eqn_number(j);
        if (eqn_number >= 0)
        {
          node_pt->set_value(t, j, *(dof_pt[eqn_number]));
        }
      }
    }

    // and non-nodal data inside elements
    for (unsigned i = 0, ni = mesh_pt()->nelement(); i < ni; i++)
    {
      GeneralisedElement* ele_pt = mesh_pt()->element_pt(i);
      for (unsigned j = 0, nj = ele_pt->ninternal_data(); j < nj; j++)
      {
        Data* data_pt = ele_pt->internal_data_pt(j);
        // For each node get the equation number and copy in the value.
        int eqn_number = data_pt->eqn_number(j);
        if (eqn_number >= 0)
        {
          data_pt->set_value(t, j, *(dof_pt[eqn_number]));
        }
      }
    }
  }


  //===================================================================
  /// Function that adds the values to the dofs
  //==================================================================
  void Problem::add_to_dofs(const double& lambda,
                            const DoubleVector& increment_dofs)
  {
    const unsigned long n_dof = this->ndof();
    for (unsigned long l = 0; l < n_dof; l++)
    {
      *Dof_pt[l] += lambda * increment_dofs[l];
    }
  }


  //=========================================================================
  /// Return the residual vector multiplied by the inverse mass matrix
  /// Virtual so that it can be overloaded for mpi problems
  //=========================================================================
  void Problem::get_inverse_mass_matrix_times_residuals(DoubleVector& Mres)
  {
    // This function does not make sense for assembly handlers other than the
    // default, so complain if we try to call it with another handler

#ifdef PARANOID
    // If we are not the default, then complain
    if (this->assembly_handler_pt() != Default_assembly_handler_pt)
    {
      std::ostringstream error_stream;
      error_stream << "The function get_inverse_mass_matrix_times_residuals() "
                      "can only be\n"
                   << "used with the default assembly handler\n\n";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find the number of degrees of freedom in the problem
    const unsigned n_dof = this->ndof();

    // Resize the vector
    LinearAlgebraDistribution dist(this->communicator_pt(), n_dof, false);
    Mres.build(&dist, 0.0);

    // If we have discontinuous formulation
    // We can invert the mass matrix element by element
    if (Discontinuous_element_formulation)
    {
      // Loop over the elements and get their residuals
      const unsigned n_element = Problem::mesh_pt()->nelement();
      Vector<double> element_Mres;
      for (unsigned e = 0; e < n_element; e++)
      {
        // Cache the element
        DGElement* const elem_pt =
          dynamic_cast<DGElement*>(Problem::mesh_pt()->element_pt(e));

        // Find the elemental inverse mass matrix times residuals
        const unsigned n_el_dofs = elem_pt->ndof();
        elem_pt->get_inverse_mass_matrix_times_residuals(element_Mres);

        // Add contribution to global matrix
        for (unsigned i = 0; i < n_el_dofs; i++)
        {
          Mres[elem_pt->eqn_number(i)] = element_Mres[i];
        }
      }
    }
    // Otherwise it's continous and we must invert the full
    // mass matrix via a global linear solve.
    else
    {
      // Now do the linear solve -- recycling Mass matrix if requested
      // If we already have the factorised mass matrix, then resolve
      if (Mass_matrix_reuse_is_enabled && Mass_matrix_has_been_computed)
      {
        if (!Shut_up_in_newton_solve)
        {
          oomph_info << "Not recomputing Mass Matrix " << std::endl;
        }

        // Get the residuals
        DoubleVector residuals(&dist, 0.0);
        this->get_residuals(residuals);

        // Resolve the linear system
        this->mass_matrix_solver_for_explicit_timestepper_pt()->resolve(
          residuals, Mres);
      }
      // Otherwise solve for the first time
      else
      {
        // If we wish to reuse the mass matrix, then enable resolve
        if (Mass_matrix_reuse_is_enabled)
        {
          if (!Shut_up_in_newton_solve)
          {
            oomph_info << "Enabling resolve in explicit timestep" << std::endl;
          }
          this->mass_matrix_solver_for_explicit_timestepper_pt()
            ->enable_resolve();
        }

        // Use a custom assembly handler to assemble and invert the mass matrix

        // Store the old assembly handler
        AssemblyHandler* old_assembly_handler_pt = this->assembly_handler_pt();
        // Set the assembly handler to the explicit timestep handler
        this->assembly_handler_pt() = new ExplicitTimeStepHandler;

        // Solve the linear system
        this->mass_matrix_solver_for_explicit_timestepper_pt()->solve(this,
                                                                      Mres);
        // The mass matrix has now been computed
        Mass_matrix_has_been_computed = true;

        // Delete the Explicit Timestep handler
        delete this->assembly_handler_pt();
        // Reset the assembly handler to the original handler
        this->assembly_handler_pt() = old_assembly_handler_pt;
      }
    }
  }

  void Problem::get_dvaluesdt(DoubleVector& f)
  {
    // Loop over timesteppers: make them (temporarily) steady and store their
    // is_steady status.
    unsigned n_time_steppers = this->ntime_stepper();
    std::vector<bool> was_steady(n_time_steppers);
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      was_steady[i] = time_stepper_pt(i)->is_steady();
      time_stepper_pt(i)->make_steady();
    }

    // Calculate f using the residual/jacobian machinary.
    get_inverse_mass_matrix_times_residuals(f);

    // Reset the is_steady status of all timesteppers that weren't already
    // steady when we came in here and reset their weights
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      if (!was_steady[i])
      {
        time_stepper_pt(i)->undo_make_steady();
      }
    }
  }


  //================================================================
  /// Get the total residuals Vector for the problem
  //================================================================
  void Problem::get_residuals(DoubleVector& residuals)
  {
    // Three different cases; if MPI_Helpers::MPI_has_been_initialised=true
    // this means MPI_Helpers::init() has been called.  This could happen on a
    // code compiled with MPI but run serially; in this instance the
    // get_residuals function still works on one processor.
    //
    // Secondly, if a code has been compiled with MPI, but MPI_Helpers::init()
    // has not been called, then MPI_Helpers::MPI_has_been_initialised=false
    // and the code calls...
    //
    // Thirdly, the serial version (compiled by all, but only run when compiled
    // with MPI if MPI_Helpers::MPI_has_been_initialised=false

    // Check that the residuals has the correct number of rows if it has been
    // setup
#ifdef PARANOID
    if (residuals.built())
    {
      if (residuals.distribution_pt()->nrow() != this->ndof())
      {
        std::ostringstream error_stream;
        error_stream << "The distribution of the residuals vector does not "
                        "have the correct\n"
                     << "number of global rows\n";

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Find the number of rows
    const unsigned nrow = this->ndof();

    // Determine the distribution for the residuals vector
    // IF the vector has distribution setup then use that
    // ELSE determine the distribution based on the
    // distributed_matrix_distribution enum
    LinearAlgebraDistribution* dist_pt = 0;
    if (residuals.built())
    {
      dist_pt = new LinearAlgebraDistribution(residuals.distribution_pt());
    }
    else
    {
#ifdef OOMPH_HAS_MPI
      // if problem is only one one processor assemble non-distributed
      // distribution
      if (Communicator_pt->nproc() == 1)
      {
        dist_pt = new LinearAlgebraDistribution(Communicator_pt, nrow, false);
      }
      // if the problem is not distributed then assemble the jacobian with
      // a uniform distributed distribution
      else if (!Problem_has_been_distributed)
      {
        dist_pt = new LinearAlgebraDistribution(Communicator_pt, nrow, true);
      }
      // otherwise the problem is a distributed problem
      else
      {
        switch (Dist_problem_matrix_distribution)
        {
          case Uniform_matrix_distribution:
            dist_pt =
              new LinearAlgebraDistribution(Communicator_pt, nrow, true);
            break;
          case Problem_matrix_distribution:
            dist_pt = new LinearAlgebraDistribution(Dof_distribution_pt);
            break;
          case Default_matrix_distribution:
            LinearAlgebraDistribution* uniform_dist_pt =
              new LinearAlgebraDistribution(Communicator_pt, nrow, true);
            bool use_problem_dist = true;
            unsigned nproc = Communicator_pt->nproc();
            for (unsigned p = 0; p < nproc; p++)
            {
              if ((double)Dof_distribution_pt->nrow_local(p) >
                  ((double)uniform_dist_pt->nrow_local(p)) * 1.1)
              {
                use_problem_dist = false;
              }
            }
            if (use_problem_dist)
            {
              dist_pt = new LinearAlgebraDistribution(Dof_distribution_pt);
            }
            else
            {
              dist_pt = new LinearAlgebraDistribution(uniform_dist_pt);
            }
            delete uniform_dist_pt;
            break;
        }
      }
#else
      dist_pt = new LinearAlgebraDistribution(Communicator_pt, nrow, false);
#endif
    }

    // Locally cache pointer to assembly handler
    AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

    // Build and zero the residuals
    residuals.build(dist_pt, 0.0);

    // Serial (or one processor case)
#ifdef OOMPH_HAS_MPI
    if (this->communicator_pt()->nproc() == 1)
    {
#endif // OOMPH_HAS_MPI
      // Loop over all the elements
      unsigned long Element_pt_range = Mesh_pt->nelement();
      for (unsigned long e = 0; e < Element_pt_range; e++)
      {
        // Get the pointer to the element
        GeneralisedElement* elem_pt = Mesh_pt->element_pt(e);
        // Find number of dofs in the element
        unsigned n_element_dofs = assembly_handler_pt->ndof(elem_pt);
        // Set up an array
        Vector<double> element_residuals(n_element_dofs);
        // Fill the array
        assembly_handler_pt->get_residuals(elem_pt, element_residuals);
        // Now loop over the dofs and assign values to global Vector
        for (unsigned l = 0; l < n_element_dofs; l++)
        {
          residuals[assembly_handler_pt->eqn_number(elem_pt, l)] +=
            element_residuals[l];
        }
      }
      // Otherwise parallel case
#ifdef OOMPH_HAS_MPI
    }
    else
    {
      // Store the current assembly handler
      AssemblyHandler* const old_assembly_handler_pt = Assembly_handler_pt;
      // Create a new assembly handler that only assembles the residuals
      Assembly_handler_pt =
        new ParallelResidualsHandler(old_assembly_handler_pt);

      // Setup memory for parallel sparse assemble
      // No matrix so all size zero
      Vector<int*> column_index;
      Vector<int*> row_start;
      Vector<double*> value;
      Vector<unsigned> nnz;
      // One set of residuals of sizer one
      Vector<double*> res(1);

      // Call the parallel sparse assemble, that should only assemble residuals
      parallel_sparse_assemble(
        dist_pt, column_index, row_start, value, nnz, res);
      // Fill in the residuals data
      residuals.set_external_values(res[0], true);

      // Delete new assembly handler
      delete Assembly_handler_pt;
      // Reset the assembly handler to the original
      Assembly_handler_pt = old_assembly_handler_pt;
    }
#endif

    // Delete the distribution
    delete dist_pt;
  }

  //=============================================================================
  /// Get the fully assembled residual vector and Jacobian matrix
  /// in dense storage. The DoubleVector residuals returned will be
  /// non-distributed. If on calling this method the DoubleVector residuals is
  /// setup then it must be non-distributed and of the correct length.
  /// The matrix type DenseDoubleMatrix is not distributable and therefore
  /// the residual vector is also assumed to be non distributable.
  //=============================================================================
  void Problem::get_jacobian(DoubleVector& residuals,
                             DenseDoubleMatrix& jacobian)
  {
    // get the number of degrees of freedom
    unsigned n_dof = ndof();

#ifdef PARANOID
    // PARANOID checks : if the distribution of residuals is setup then it must
    // must not be distributed, have the right number of rows, and the same
    // communicator as the problem
    if (residuals.built())
    {
      if (residuals.distribution_pt()->distributed())
      {
        std::ostringstream error_stream;
        error_stream
          << "If the DoubleVector residuals is setup then it must not "
          << "be distributed.";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      if (residuals.distribution_pt()->nrow() != n_dof)
      {
        std::ostringstream error_stream;
        error_stream
          << "If the DoubleVector residuals is setup then it must have"
          << " the correct number of rows";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      if (!(*Communicator_pt ==
            *residuals.distribution_pt()->communicator_pt()))
      {
        std::ostringstream error_stream;
        error_stream
          << "If the DoubleVector residuals is setup then it must have"
          << " the same communicator as the problem.";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // set the residuals distribution if it is not setup
    if (!residuals.built())
    {
      LinearAlgebraDistribution dist(Communicator_pt, n_dof, false);
      residuals.build(&dist, 0.0);
    }
    // else just zero the residuals
    else
    {
      residuals.initialise(0.0);
    }

    // Resize the matrices -- this cannot always be done externally
    // because get_jacobian exists in many different versions for
    // different storage formats -- resizing a CC or CR matrix doesn't
    // make sense.

    // resize the jacobian
    jacobian.resize(n_dof, n_dof);
    jacobian.initialise(0.0);

    // Locally cache pointer to assembly handler
    AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

    // Loop over all the elements
    unsigned long n_element = Mesh_pt->nelement();
    for (unsigned long e = 0; e < n_element; e++)
    {
      // Get the pointer to the element
      GeneralisedElement* elem_pt = Mesh_pt->element_pt(e);
      // Find number of dofs in the element
      unsigned n_element_dofs = assembly_handler_pt->ndof(elem_pt);
      // Set up an array
      Vector<double> element_residuals(n_element_dofs);
      // Set up a matrix
      DenseMatrix<double> element_jacobian(n_element_dofs);
      // Fill the array
      assembly_handler_pt->get_jacobian(
        elem_pt, element_residuals, element_jacobian);
      // Now loop over the dofs and assign values to global Vector
      for (unsigned l = 0; l < n_element_dofs; l++)
      {
        unsigned long eqn_number = assembly_handler_pt->eqn_number(elem_pt, l);
        residuals[eqn_number] += element_residuals[l];
        for (unsigned l2 = 0; l2 < n_element_dofs; l2++)
        {
          jacobian(eqn_number, assembly_handler_pt->eqn_number(elem_pt, l2)) +=
            element_jacobian(l, l2);
        }
      }
    }
  }

  //=============================================================================
  /// Return the fully-assembled Jacobian and residuals for the problem,
  /// in the case where the Jacobian matrix is in a distributable
  /// row compressed storage format.
  /// 1. If the distribution of the jacobian and residuals is setup then, they
  /// will be returned with that distribution.
  /// Note. the jacobian and residuals must have the same distribution.
  /// 2. If the distribution of the jacobian and residuals are not setup then
  /// their distribution will computed based on:
  /// Distributed_problem_matrix_distribution.
  //=============================================================================
  void Problem::get_jacobian(DoubleVector& residuals, CRDoubleMatrix& jacobian)
  {
    // Three different cases; if MPI_Helpers::MPI_has_been_initialised=true
    // this means MPI_Helpers::setup() has been called.  This could happen on a
    // code compiled with MPI but run serially; in this instance the
    // get_residuals function still works on one processor.
    //
    // Secondly, if a code has been compiled with MPI, but MPI_Helpers::setup()
    // has not been called, then MPI_Helpers::MPI_has_been_initialised=false
    // and the code calls...
    //
    // Thirdly, the serial version (compiled by all, but only run when compiled
    // with MPI if MPI_Helpers::MPI_has_been_initialised=false
    //
    // The only case where an MPI code cannot run serially at present
    // is one where the distribute function is used (i.e. METIS is called)

    // Allocate storage for the matrix entries
    // The generalised Vector<Vector<>> structure is required
    // for the most general interface to sparse_assemble() which allows
    // the assembly of multiple matrices at once.
    Vector<int*> column_index(1);
    Vector<int*> row_start(1);
    Vector<double*> value(1);
    Vector<unsigned> nnz(1);

#ifdef PARANOID
    // PARANOID checks that the distribution of the jacobian matches that of the
    // residuals (if they are setup) and that they have the right number of rows
    if (residuals.built() && jacobian.distribution_built())
    {
      if (!(*residuals.distribution_pt() == *jacobian.distribution_pt()))
      {
        std::ostringstream error_stream;
        error_stream << "The distribution of the residuals must "
                     << "be the same as the distribution of the jacobian.";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      if (jacobian.distribution_pt()->nrow() != this->ndof())
      {
        std::ostringstream error_stream;
        error_stream
          << "The distribution of the jacobian and residuals does not"
          << "have the correct number of global rows.";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
    else if (residuals.built() != jacobian.distribution_built())
    {
      std::ostringstream error_stream;
      error_stream << "The distribution of the jacobian and residuals must "
                   << "both be setup or both not setup";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // Allocate generalised storage format for passing to sparse_assemble()
    Vector<double*> res(1);

    // number of rows
    unsigned nrow = this->ndof();

    // determine the distribution for the jacobian.
    // IF the jacobian has distribution setup then use that
    // ELSE determine the distribution based on the
    // distributed_matrix_distribution enum
    LinearAlgebraDistribution* dist_pt = 0;
    if (jacobian.distribution_built())
    {
      dist_pt = new LinearAlgebraDistribution(jacobian.distribution_pt());
    }
    else
    {
#ifdef OOMPH_HAS_MPI
      // if problem is only one one processor
      if (Communicator_pt->nproc() == 1)
      {
        dist_pt = new LinearAlgebraDistribution(Communicator_pt, nrow, false);
      }
      // if the problem is not distributed then assemble the jacobian with
      // a uniform distributed distribution
      else if (!Problem_has_been_distributed)
      {
        dist_pt = new LinearAlgebraDistribution(Communicator_pt, nrow, true);
      }
      // otherwise the problem is a distributed problem
      else
      {
        switch (Dist_problem_matrix_distribution)
        {
          case Uniform_matrix_distribution:
            dist_pt =
              new LinearAlgebraDistribution(Communicator_pt, nrow, true);
            break;
          case Problem_matrix_distribution:
            dist_pt = new LinearAlgebraDistribution(Dof_distribution_pt);
            break;
          case Default_matrix_distribution:
            LinearAlgebraDistribution* uniform_dist_pt =
              new LinearAlgebraDistribution(Communicator_pt, nrow, true);
            bool use_problem_dist = true;
            unsigned nproc = Communicator_pt->nproc();
            for (unsigned p = 0; p < nproc; p++)
            {
              if ((double)Dof_distribution_pt->nrow_local(p) >
                  ((double)uniform_dist_pt->nrow_local(p)) * 1.1)
              {
                use_problem_dist = false;
              }
            }
            if (use_problem_dist)
            {
              dist_pt = new LinearAlgebraDistribution(Dof_distribution_pt);
            }
            else
            {
              dist_pt = new LinearAlgebraDistribution(uniform_dist_pt);
            }
            delete uniform_dist_pt;
            break;
        }
      }
#else
      dist_pt = new LinearAlgebraDistribution(Communicator_pt, nrow, false);
#endif
    }


    // The matrix is in compressed row format
    bool compressed_row_flag = true;

#ifdef OOMPH_HAS_MPI
    //
    if (Communicator_pt->nproc() == 1)
    {
#endif
      sparse_assemble_row_or_column_compressed(
        column_index, row_start, value, nnz, res, compressed_row_flag);
      jacobian.build(dist_pt);
      jacobian.build_without_copy(
        dist_pt->nrow(), nnz[0], value[0], column_index[0], row_start[0]);
      residuals.build(dist_pt, 0.0);
      residuals.set_external_values(res[0], true);
#ifdef OOMPH_HAS_MPI
    }
    else
    {
      if (dist_pt->distributed())
      {
        parallel_sparse_assemble(
          dist_pt, column_index, row_start, value, nnz, res);
        jacobian.build(dist_pt);
        jacobian.build_without_copy(
          dist_pt->nrow(), nnz[0], value[0], column_index[0], row_start[0]);
        residuals.build(dist_pt, 0.0);
        residuals.set_external_values(res[0], true);
      }
      else
      {
        LinearAlgebraDistribution* temp_dist_pt =
          new LinearAlgebraDistribution(Communicator_pt, dist_pt->nrow(), true);
        parallel_sparse_assemble(
          temp_dist_pt, column_index, row_start, value, nnz, res);
        jacobian.build(temp_dist_pt);
        jacobian.build_without_copy(
          dist_pt->nrow(), nnz[0], value[0], column_index[0], row_start[0]);
        jacobian.redistribute(dist_pt);
        residuals.build(temp_dist_pt, 0.0);
        residuals.set_external_values(res[0], true);
        residuals.redistribute(dist_pt);
        delete temp_dist_pt;
      }
    }
#endif

    // clean up dist_pt and residuals_vector pt
    delete dist_pt;
  }

  //=============================================================================
  /// Return the fully-assembled Jacobian and residuals for the problem,
  /// in the case when the jacobian matrix is in column-compressed storage
  /// format.
  //=============================================================================
  void Problem::get_jacobian(DoubleVector& residuals, CCDoubleMatrix& jacobian)
  {
    // Three different cases; if MPI_Helpers::MPI_has_been_initialised=true
    // this means MPI_Helpers::setup() has been called.  This could happen on a
    // code compiled with MPI but run serially; in this instance the
    // get_residuals function still works on one processor.
    //
    // Secondly, if a code has been compiled with MPI, but MPI_Helpers::setup()
    // has not been called, then MPI_Helpers::MPI_has_been_initialised=false
    // and the code calls...
    //
    // Thirdly, the serial version (compiled by all, but only run when compiled
    // with MPI if MPI_Helpers::MPI_has_been_5Binitialised=false
    //
    // The only case where an MPI code cannot run serially at present
    // is one where the distribute function is used (i.e. METIS is called)

    // get the number of degrees of freedom
    unsigned n_dof = ndof();

#ifdef PARANOID
    // PARANOID checks : if the distribution of residuals is setup then it must
    // must not be distributed, have the right number of rows, and the same
    // communicator as the problem
    if (residuals.built())
    {
      if (residuals.distribution_pt()->distributed())
      {
        std::ostringstream error_stream;
        error_stream
          << "If the DoubleVector residuals is setup then it must not "
          << "be distributed.";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      if (residuals.distribution_pt()->nrow() != n_dof)
      {
        std::ostringstream error_stream;
        error_stream
          << "If the DoubleVector residuals is setup then it must have"
          << " the correct number of rows";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      if (!(*Communicator_pt ==
            *residuals.distribution_pt()->communicator_pt()))
      {
        std::ostringstream error_stream;
        error_stream
          << "If the DoubleVector residuals is setup then it must have"
          << " the same communicator as the problem.";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Allocate storage for the matrix entries
    // The generalised Vector<Vector<>> structure is required
    // for the most general interface to sparse_assemble() which allows
    // the assembly of multiple matrices at once.
    Vector<int*> row_index(1);
    Vector<int*> column_start(1);
    Vector<double*> value(1);

    // Allocate generalised storage format for passing to sparse_assemble()
    Vector<double*> res(1);

    // allocate storage for the number of non-zeros in each matrix
    Vector<unsigned> nnz(1);

    // The matrix is in compressed column format
    bool compressed_row_flag = false;

    // get the distribution for the residuals
    LinearAlgebraDistribution* dist_pt;
    if (!residuals.built())
    {
      dist_pt =
        new LinearAlgebraDistribution(Communicator_pt, this->ndof(), false);
    }
    else
    {
      dist_pt = new LinearAlgebraDistribution(residuals.distribution_pt());
    }

#ifdef OOMPH_HAS_MPI
    if (communicator_pt()->nproc() == 1)
    {
#endif
      sparse_assemble_row_or_column_compressed(
        row_index, column_start, value, nnz, res, compressed_row_flag);
      jacobian.build_without_copy(
        value[0], row_index[0], column_start[0], nnz[0], n_dof, n_dof);
      residuals.build(dist_pt, 0.0);
      residuals.set_external_values(res[0], true);
#ifdef OOMPH_HAS_MPI
    }
    else
    {
      std::ostringstream error_stream;
      error_stream << "Cannot assemble a CCDoubleMatrix Jacobian on more "
                   << "than one processor.";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // clean up
    delete dist_pt;
  }


  //===================================================================
  /// \short Set all pinned values to zero.
  /// Used to set boundary conditions to be homogeneous in the copy
  /// of the problem  used in adaptive bifurcation tracking
  /// (ALH: TEMPORARY HACK, WILL BE FIXED)
  //==================================================================
  void Problem::set_pinned_values_to_zero()
  {
    // NOTE THIS DOES NOT ZERO ANY SPINE DATA, but otherwise everything else
    // should be zeroed

    // Zero any pinned global Data
    const unsigned n_global_data = nglobal_data();
    for (unsigned i = 0; i < n_global_data; i++)
    {
      Data* const local_data_pt = Global_data_pt[i];
      const unsigned n_value = local_data_pt->nvalue();
      for (unsigned j = 0; j < n_value; j++)
      {
        // If the data value is pinned set the value to zero
        if (local_data_pt->is_pinned(j))
        {
          local_data_pt->set_value(j, 0.0);
        }
      }
    }

    // Loop over the submeshes:
    const unsigned n_sub_mesh = Sub_mesh_pt.size();
    if (n_sub_mesh == 0)
    {
      // Loop over the nodes in the element
      const unsigned n_node = Mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        Node* const local_node_pt = Mesh_pt->node_pt(n);
        const unsigned n_value = local_node_pt->nvalue();
        for (unsigned j = 0; j < n_value; j++)
        {
          // If the data value is pinned set the value to zero
          if (local_node_pt->is_pinned(j))
          {
            local_node_pt->set_value(j, 0.0);
          }
        }

        // Try to cast to a solid node
        SolidNode* const local_solid_node_pt =
          dynamic_cast<SolidNode*>(local_node_pt);
        // If we are successful
        if (local_solid_node_pt)
        {
          // Find the dimension of the node
          const unsigned n_dim = local_solid_node_pt->ndim();
          // Find number of positions
          const unsigned n_position_type =
            local_solid_node_pt->nposition_type();

          for (unsigned k = 0; k < n_position_type; k++)
          {
            for (unsigned i = 0; i < n_dim; i++)
            {
              // If the generalised position is pinned,
              // set the value to zero
              if (local_solid_node_pt->position_is_pinned(k, i))
              {
                local_solid_node_pt->x_gen(k, i) = 0.0;
              }
            }
          }
        }
      }

      // Now loop over the element's and zero the internal data
      const unsigned n_element = Mesh_pt->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        GeneralisedElement* const local_element_pt = Mesh_pt->element_pt(e);
        const unsigned n_internal = local_element_pt->ninternal_data();
        for (unsigned i = 0; i < n_internal; i++)
        {
          Data* const local_data_pt = local_element_pt->internal_data_pt(i);
          const unsigned n_value = local_data_pt->nvalue();
          for (unsigned j = 0; j < n_value; j++)
          {
            // If the data value is pinned set the value to zero
            if (local_data_pt->is_pinned(j))
            {
              local_data_pt->set_value(j, 0.0);
            }
          }
        }
      } // End of loop over elements
    }
    else
    {
      // Alternatively loop over all sub meshes
      for (unsigned m = 0; m < n_sub_mesh; m++)
      {
        // Loop over the nodes in the element
        const unsigned n_node = Sub_mesh_pt[m]->nnode();
        for (unsigned n = 0; n < n_node; n++)
        {
          Node* const local_node_pt = Sub_mesh_pt[m]->node_pt(n);
          const unsigned n_value = local_node_pt->nvalue();
          for (unsigned j = 0; j < n_value; j++)
          {
            // If the data value is pinned set the value to zero
            if (local_node_pt->is_pinned(j))
            {
              local_node_pt->set_value(j, 0.0);
            }
          }

          // Try to cast to a solid node
          SolidNode* const local_solid_node_pt =
            dynamic_cast<SolidNode*>(local_node_pt);
          // If we are successful
          if (local_solid_node_pt)
          {
            // Find the dimension of the node
            const unsigned n_dim = local_solid_node_pt->ndim();
            // Find number of positions
            const unsigned n_position_type =
              local_solid_node_pt->nposition_type();

            for (unsigned k = 0; k < n_position_type; k++)
            {
              for (unsigned i = 0; i < n_dim; i++)
              {
                // If the generalised position is pinned,
                // set the value to zero
                if (local_solid_node_pt->position_is_pinned(k, i))
                {
                  local_solid_node_pt->x_gen(k, i) = 0.0;
                }
              }
            }
          }
        }

        // Now loop over the element's and zero the internal data
        const unsigned n_element = Sub_mesh_pt[m]->nelement();
        for (unsigned e = 0; e < n_element; e++)
        {
          GeneralisedElement* const local_element_pt =
            Sub_mesh_pt[m]->element_pt(e);
          const unsigned n_internal = local_element_pt->ninternal_data();
          for (unsigned i = 0; i < n_internal; i++)
          {
            Data* const local_data_pt = local_element_pt->internal_data_pt(i);
            const unsigned n_value = local_data_pt->nvalue();
            for (unsigned j = 0; j < n_value; j++)
            {
              // If the data value is pinned set the value to zero
              if (local_data_pt->is_pinned(j))
              {
                local_data_pt->set_value(j, 0.0);
              }
            }
          }
        } // End of loop over elements
      }
    }
  }


  //=====================================================================
  /// This is a (private) helper function that is used to assemble system
  /// matrices in compressed row or column format
  /// and compute residual vectors.
  /// The default action is to assemble the jacobian matrix and
  /// residuals for the Newton method. The action can be
  /// overloaded at an elemental level by changing the default
  /// behaviour of the function Element::get_all_vectors_and_matrices().
  /// column_or_row_index: Column [or row] index of given entry
  /// row_or_column_start: Index of first entry for given row [or column]
  /// value              : Vector of nonzero entries
  /// residuals          : Residual vector
  /// compressed_row_flag: Bool flag to indicate if storage format is
  ///                      compressed row [if false interpretation of
  ///                      arguments is as stated in square brackets].
  /// We provide four different assembly methods, each with different
  /// memory requirements/execution speeds. The method is set by
  /// the public flag Problem::Sparse_assembly_method.
  //=====================================================================
  void Problem::sparse_assemble_row_or_column_compressed(
    Vector<int*>& column_or_row_index,
    Vector<int*>& row_or_column_start,
    Vector<double*>& value,
    Vector<unsigned>& nnz,
    Vector<double*>& residuals,
    bool compressed_row_flag)
  {
    // Choose the actual method
    switch (Sparse_assembly_method)
    {
      case Perform_assembly_using_vectors_of_pairs:

        sparse_assemble_row_or_column_compressed_with_vectors_of_pairs(
          column_or_row_index,
          row_or_column_start,
          value,
          nnz,
          residuals,
          compressed_row_flag);

        break;

      case Perform_assembly_using_two_vectors:

        sparse_assemble_row_or_column_compressed_with_two_vectors(
          column_or_row_index,
          row_or_column_start,
          value,
          nnz,
          residuals,
          compressed_row_flag);

        break;

      case Perform_assembly_using_maps:

        sparse_assemble_row_or_column_compressed_with_maps(column_or_row_index,
                                                           row_or_column_start,
                                                           value,
                                                           nnz,
                                                           residuals,
                                                           compressed_row_flag);

        break;

      case Perform_assembly_using_lists:

        sparse_assemble_row_or_column_compressed_with_lists(
          column_or_row_index,
          row_or_column_start,
          value,
          nnz,
          residuals,
          compressed_row_flag);

        break;

      case Perform_assembly_using_two_arrays:

        sparse_assemble_row_or_column_compressed_with_two_arrays(
          column_or_row_index,
          row_or_column_start,
          value,
          nnz,
          residuals,
          compressed_row_flag);

        break;

      default:

        std::ostringstream error_stream;
        error_stream
          << "Error: Incorrect value for Problem::Sparse_assembly_method"
          << Sparse_assembly_method << std::endl
          << "It should be one of the enumeration Problem::Assembly_method"
          << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //=====================================================================
  /// This is a (private) helper function that is used to assemble system
  /// matrices in compressed row or column format
  /// and compute residual vectors, using maps
  /// The default action is to assemble the jacobian matrix and
  /// residuals for the Newton method. The action can be
  /// overloaded at an elemental level by chaging the default
  /// behaviour of the function Element::get_all_vectors_and_matrices().
  /// column_or_row_index: Column [or row] index of given entry
  /// row_or_column_start: Index of first entry for given row [or column]
  /// value              : Vector of nonzero entries
  /// residuals          : Residual vector
  /// compressed_row_flag: Bool flag to indicate if storage format is
  ///                      compressed row [if false interpretation of
  ///                      arguments is as stated in square brackets].
  //=====================================================================
  void Problem::sparse_assemble_row_or_column_compressed_with_maps(
    Vector<int*>& column_or_row_index,
    Vector<int*>& row_or_column_start,
    Vector<double*>& value,
    Vector<unsigned>& nnz,
    Vector<double*>& residuals,
    bool compressed_row_flag)
  {
    // Total number of elements
    const unsigned long n_elements = mesh_pt()->nelement();

    // Default range of elements for distributed problems
    unsigned long el_lo = 0;
    unsigned long el_hi = n_elements - 1;

#ifdef OOMPH_HAS_MPI
    // Otherwise just loop over a fraction of the elements
    // (This will either have been initialised in
    // Problem::set_default_first_and_last_element_for_assembly() or
    // will have been re-assigned during a previous assembly loop
    // Note that following the re-assignment only the entries
    // for the current processor are relevant.
    if (!Problem_has_been_distributed)
    {
      el_lo = First_el_for_assembly[Communicator_pt->my_rank()];
      el_hi = Last_el_plus_one_for_assembly[Communicator_pt->my_rank()] - 1;
    }
#endif

    // number of dofs
    unsigned ndof = this->ndof();

    // Find the number of vectors to be assembled
    const unsigned n_vector = residuals.size();

    // Find the number of matrices to be assembled
    const unsigned n_matrix = column_or_row_index.size();

    // Locally cache pointer to assembly handler
    AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

#ifdef OOMPH_HAS_MPI
    bool doing_residuals = false;
    if (dynamic_cast<ParallelResidualsHandler*>(Assembly_handler_pt) != 0)
    {
      doing_residuals = true;
    }
#endif

// Error check dimensions
#ifdef PARANOID
    if (row_or_column_start.size() != n_matrix)
    {
      std::ostringstream error_stream;
      error_stream << "Error: " << std::endl
                   << "row_or_column_start.size() "
                   << row_or_column_start.size() << " does not equal "
                   << "column_or_row_index.size() "
                   << column_or_row_index.size() << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    if (value.size() != n_matrix)
    {
      std::ostringstream error_stream;
      error_stream
        << "Error in Problem::sparse_assemble_row_or_column_compressed "
        << std::endl
        << "value.size() " << value.size() << " does not equal "
        << "column_or_row_index.size() " << column_or_row_index.size()
        << std::endl
        << std::endl
        << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // The idea behind this sparse assembly routine is to use a vector of
    // maps for the entries in each row or column of the complete matrix.
    // The key for each map is the global row or column number and
    // the default comparison operator for integers means that each map
    // is ordered by the global row or column number. Thus, we need not
    // sort the maps, that happens at each insertion of a new entry. The
    // price we pay  is that for large maps, inseration is not a
    // cheap operation. Hash maps can be used to increase the speed, but then
    // the ordering is lost and we would have to sort anyway. The solution if
    // speed is required is to use lists, see below.


    // Set up a vector of vectors of maps of entries of each  matrix,
    // indexed by either the column or row. The entries of the vector for
    // each matrix correspond to all the rows or columns of that matrix.
    // The use of the map storage
    // scheme, with its implicit ordering on the first index, gives
    // a sparse ordered list of the entries in the given row or column.
    Vector<Vector<std::map<unsigned, double>>> matrix_data_map(n_matrix);
    // Loop over the number of matrices being assembled and resize
    // each vector of maps to the number of rows or columns of the matrix
    for (unsigned m = 0; m < n_matrix; m++)
    {
      matrix_data_map[m].resize(ndof);
    }

    // Resize the residuals vectors
    for (unsigned v = 0; v < n_vector; v++)
    {
      residuals[v] = new double[ndof];
      for (unsigned i = 0; i < ndof; i++)
      {
        residuals[v][i] = 0;
      }
    }


#ifdef OOMPH_HAS_MPI


    // Storage for assembly time for elements
    double t_assemble_start = 0.0;

    // Storage for assembly times
    if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
    {
      Elemental_assembly_time.resize(n_elements);
    }

#endif

    //----------------Assemble and populate the maps-------------------------
    {
      // Allocate local storage for the element's contribution to the
      // residuals vectors and system matrices of the size of the maximum
      // number of dofs in any element.
      // This means that the storage is only allocated (and deleted) once
      Vector<Vector<double>> el_residuals(n_vector);
      Vector<DenseMatrix<double>> el_jacobian(n_matrix);

      // Loop over the elements for this processor
      for (unsigned long e = el_lo; e <= el_hi; e++)
      {
#ifdef OOMPH_HAS_MPI
        // Time it?
        if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
        {
          t_assemble_start = TimingHelpers::timer();
        }
#endif

        // Get the pointer to the element
        GeneralisedElement* elem_pt = mesh_pt()->element_pt(e);

#ifdef OOMPH_HAS_MPI
        // Ignore halo elements
        if (!elem_pt->is_halo())
        {
#endif

          // Find number of degrees of freedom in the element
          const unsigned nvar = assembly_handler_pt->ndof(elem_pt);

          // Resize the storage for elemental jacobian and residuals
          for (unsigned v = 0; v < n_vector; v++)
          {
            el_residuals[v].resize(nvar);
          }
          for (unsigned m = 0; m < n_matrix; m++)
          {
            el_jacobian[m].resize(nvar);
          }

          // Now get the residuals and jacobian for the element
          assembly_handler_pt->get_all_vectors_and_matrices(
            elem_pt, el_residuals, el_jacobian);

          //---------------Insert the values into the maps--------------

          // Loop over the first index of local variables
          for (unsigned i = 0; i < nvar; i++)
          {
            // Get the local equation number
            unsigned eqn_number = assembly_handler_pt->eqn_number(elem_pt, i);

            // Add the contribution to the residuals
            for (unsigned v = 0; v < n_vector; v++)
            {
              // Fill in each residuals vector
              residuals[v][eqn_number] += el_residuals[v][i];
            }

            // Now loop over the other index
            for (unsigned j = 0; j < nvar; j++)
            {
              // Get the number of the unknown
              unsigned unknown = assembly_handler_pt->eqn_number(elem_pt, j);

              // Loop over the matrices
              for (unsigned m = 0; m < n_matrix; m++)
              {
                // Get the value of the matrix at this point
                double value = el_jacobian[m](i, j);
                // Only bother to add to the map if it's non-zero
                if (std::fabs(value) > Numerical_zero_for_sparse_assembly)
                {
                  // If it's compressed row storage, then our vector of maps
                  // is indexed by row (equation number)
                  if (compressed_row_flag)
                  {
                    // Add the data into the map using the unknown as the map
                    // key
                    matrix_data_map[m][eqn_number][unknown] += value;
                  }
                  // Otherwise it's compressed column storage and our vector is
                  // indexed by column (the unknown)
                  else
                  {
                    // Add the data into the map using the eqn_numbe as the map
                    // key
                    matrix_data_map[m][unknown][eqn_number] += value;
                  }
                }
              } // End of loop over matrices
            }
          }

#ifdef OOMPH_HAS_MPI
        } // endif halo element
#endif


#ifdef OOMPH_HAS_MPI
        // Time it?
        if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
        {
          Elemental_assembly_time[e] =
            TimingHelpers::timer() - t_assemble_start;
        }
#endif

      } // End of loop over the elements

    } // End of map assembly


#ifdef OOMPH_HAS_MPI

    // Postprocess timing information and re-allocate distribution of
    // elements during subsequent assemblies.
    if ((!doing_residuals) && (!Problem_has_been_distributed) &&
        Must_recompute_load_balance_for_assembly)
    {
      recompute_load_balanced_assembly();
    }

    // We have determined load balancing for current setup.
    // This can remain the same until assign_eqn_numbers() is called
    // again -- the flag is re-set to true there.
    if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
    {
      Must_recompute_load_balance_for_assembly = false;
    }

#endif


    //-----------Finally we need to convert the beautiful map storage scheme
    //------------------------to the containers required by SuperLU

    // Loop over the number of matrices
    for (unsigned m = 0; m < n_matrix; m++)
    {
      // Set the number of rows or columns
      row_or_column_start[m] = new int[ndof + 1];
      // Counter for the total number of entries in the storage scheme
      unsigned long entry_count = 0;
      row_or_column_start[m][0] = entry_count;

      // first we compute the number of non-zeros
      nnz[m] = 0;
      for (unsigned long i_global = 0; i_global < ndof; i_global++)
      {
        nnz[m] += matrix_data_map[m][i_global].size();
      }

      // and then resize the storage
      column_or_row_index[m] = new int[nnz[m]];
      value[m] = new double[nnz[m]];

      // Now we merely loop over the number of rows or columns
      for (unsigned long i_global = 0; i_global < ndof; i_global++)
      {
        // Start index for the present row
        row_or_column_start[m][i_global] = entry_count;
        // If there are no entries in the map then skip the rest of the loop
        if (matrix_data_map[m][i_global].empty())
        {
          continue;
        }

        // Loop over all the entries in the map corresponding to the given
        // row or column. It will be ordered

        for (std::map<unsigned, double>::iterator it =
               matrix_data_map[m][i_global].begin();
             it != matrix_data_map[m][i_global].end();
             ++it)
        {
          // The first value is the column or row index
          column_or_row_index[m][entry_count] = it->first;
          // The second value is the actual data value
          value[m][entry_count] = it->second;
          // Increase the value of the counter
          entry_count++;
        }
      }

      // Final entry in the row/column start vector
      row_or_column_start[m][ndof] = entry_count;
    } // End of the loop over the matrices

    if (Pause_at_end_of_sparse_assembly)
    {
      oomph_info << "Pausing at end of sparse assembly." << std::endl;
      pause("Check memory usage now.");
    }
  }


  //=====================================================================
  /// This is a (private) helper function that is used to assemble system
  /// matrices in compressed row or column format
  /// and compute residual vectors using lists
  /// The default action is to assemble the jacobian matrix and
  /// residuals for the Newton method. The action can be
  /// overloaded at an elemental level by chaging the default
  /// behaviour of the function Element::get_all_vectors_and_matrices().
  /// column_or_row_index: Column [or row] index of given entry
  /// row_or_column_start: Index of first entry for given row [or column]
  /// value              : Vector of nonzero entries
  /// residuals          : Residual vector
  /// compressed_row_flag: Bool flag to indicate if storage format is
  ///                      compressed row [if false interpretation of
  ///                      arguments is as stated in square brackets].
  //=====================================================================
  void Problem::sparse_assemble_row_or_column_compressed_with_lists(
    Vector<int*>& column_or_row_index,
    Vector<int*>& row_or_column_start,
    Vector<double*>& value,
    Vector<unsigned>& nnz,
    Vector<double*>& residuals,
    bool compressed_row_flag)
  {
    // Total number of elements
    const unsigned long n_elements = mesh_pt()->nelement();

    // Default range of elements for distributed problems
    unsigned long el_lo = 0;
    unsigned long el_hi = n_elements - 1;

#ifdef OOMPH_HAS_MPI
    // Otherwise just loop over a fraction of the elements
    // (This will either have been initialised in
    // Problem::set_default_first_and_last_element_for_assembly() or
    // will have been re-assigned during a previous assembly loop
    // Note that following the re-assignment only the entries
    // for the current processor are relevant.
    if (!Problem_has_been_distributed)
    {
      el_lo = First_el_for_assembly[Communicator_pt->my_rank()];
      el_hi = Last_el_plus_one_for_assembly[Communicator_pt->my_rank()] - 1;
    }
#endif

    // number of dofs
    unsigned ndof = this->ndof();

    // Find the number of vectors to be assembled
    const unsigned n_vector = residuals.size();

    // Find the number of matrices to be assembled
    const unsigned n_matrix = column_or_row_index.size();

    // Locally cache pointer to assembly handler
    AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

#ifdef OOMPH_HAS_MPI
    bool doing_residuals = false;
    if (dynamic_cast<ParallelResidualsHandler*>(Assembly_handler_pt) != 0)
    {
      doing_residuals = true;
    }
#endif

// Error check dimensions
#ifdef PARANOID
    if (row_or_column_start.size() != n_matrix)
    {
      std::ostringstream error_stream;
      error_stream << "Error: " << std::endl
                   << "row_or_column_start.size() "
                   << row_or_column_start.size() << " does not equal "
                   << "column_or_row_index.size() "
                   << column_or_row_index.size() << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    if (value.size() != n_matrix)
    {
      std::ostringstream error_stream;
      error_stream
        << "Error in Problem::sparse_assemble_row_or_column_compressed "
        << std::endl
        << "value.size() " << value.size() << " does not equal "
        << "column_or_row_index.size() " << column_or_row_index.size()
        << std::endl
        << std::endl
        << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // The idea behind this sparse assembly routine is to use a vector of
    // lists for the entries in each row or column of the complete matrix.
    // The lists contain pairs of entries (global row/column number, value).
    // All non-zero contributions from each element are added to the lists.
    // We then sort each list by global row/column number and then combine
    // the entries corresponding to each row/column before adding to the
    // vectors column_or_row_index and value.

    // Note the trade off for "fast assembly" is that we will require
    // more memory during the assembly phase. Then again, if we can
    // only just assemble the sparse matrix, we're in real trouble.

    // Set up a vector of lists of paired entries of
    //(row/column index, jacobian matrix entry).
    // The entries of the vector correspond to all the rows or columns.
    // The use of the list storage scheme, should give fast insertion
    // and fast sorts later.
    Vector<Vector<std::list<std::pair<unsigned, double>>>> matrix_data_list(
      n_matrix);
    // Loop over the number of matrices and resize
    for (unsigned m = 0; m < n_matrix; m++)
    {
      matrix_data_list[m].resize(ndof);
    }

    // Resize the residuals vectors
    for (unsigned v = 0; v < n_vector; v++)
    {
      residuals[v] = new double[ndof];
      for (unsigned i = 0; i < ndof; i++)
      {
        residuals[v][i] = 0;
      }
    }

#ifdef OOMPH_HAS_MPI


    // Storage for assembly time for elements
    double t_assemble_start = 0.0;

    // Storage for assembly times
    if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
    {
      Elemental_assembly_time.resize(n_elements);
    }

#endif

    //------------Assemble and populate the lists-----------------------
    {
      // Allocate local storage for the element's contribution to the
      // residuals vectors and system matrices of the size of the maximum
      // number of dofs in any element.
      // This means that the stored is only allocated (and deleted) once
      Vector<Vector<double>> el_residuals(n_vector);
      Vector<DenseMatrix<double>> el_jacobian(n_matrix);


      // Pointer to a single list to be used during the assembly
      std::list<std::pair<unsigned, double>>* list_pt;

      // Loop over the all elements
      for (unsigned long e = el_lo; e <= el_hi; e++)
      {
#ifdef OOMPH_HAS_MPI
        // Time it?
        if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
        {
          t_assemble_start = TimingHelpers::timer();
        }
#endif

        // Get the pointer to the element
        GeneralisedElement* elem_pt = mesh_pt()->element_pt(e);

#ifdef OOMPH_HAS_MPI
        // Ignore halo elements
        if (!elem_pt->is_halo())
        {
#endif

          // Find number of degrees of freedom in the element
          const unsigned nvar = assembly_handler_pt->ndof(elem_pt);

          // Resize the storage for the elemental jacobian and residuals
          for (unsigned v = 0; v < n_vector; v++)
          {
            el_residuals[v].resize(nvar);
          }
          for (unsigned m = 0; m < n_matrix; m++)
          {
            el_jacobian[m].resize(nvar);
          }

          // Now get the residuals and jacobian for the element
          assembly_handler_pt->get_all_vectors_and_matrices(
            elem_pt, el_residuals, el_jacobian);

          //---------------- Insert the values into the lists -----------

          // Loop over the first index of local variables
          for (unsigned i = 0; i < nvar; i++)
          {
            // Get the local equation number
            unsigned eqn_number = assembly_handler_pt->eqn_number(elem_pt, i);

            // Add the contribution to the residuals
            for (unsigned v = 0; v < n_vector; v++)
            {
              // Fill in the residuals vector
              residuals[v][eqn_number] += el_residuals[v][i];
            }

            // Now loop over the other index
            for (unsigned j = 0; j < nvar; j++)
            {
              // Get the number of the unknown
              unsigned unknown = assembly_handler_pt->eqn_number(elem_pt, j);

              // Loop over the matrices
              for (unsigned m = 0; m < n_matrix; m++)
              {
                // Get the value of the matrix at this point
                double value = el_jacobian[m](i, j);
                // Only add to theif it's non-zero
                if (std::fabs(value) > Numerical_zero_for_sparse_assembly)
                {
                  // If it's compressed row storage, then our vector is indexed
                  // by row (the equation number)
                  if (compressed_row_flag)
                  {
                    // Find the list that corresponds to the desired row
                    list_pt = &matrix_data_list[m][eqn_number];
                    // Insert the data into the list, the first entry
                    // in the pair is the unknown (column index),
                    // the second is the value itself.
                    list_pt->insert(list_pt->end(),
                                    std::make_pair(unknown, value));
                  }
                  // Otherwise it's compressed column storage, and our
                  // vector is indexed by column (the unknown)
                  else
                  {
                    // Find the list that corresponds to the desired column
                    list_pt = &matrix_data_list[m][unknown];
                    // Insert the data into the list, the first entry
                    // in the pair is the equation number (row index),
                    // the second is the value itself.
                    list_pt->insert(list_pt->end(),
                                    std::make_pair(eqn_number, value));
                  }
                }
              }
            }
          }

#ifdef OOMPH_HAS_MPI
        } // endif halo element
#endif


#ifdef OOMPH_HAS_MPI
        // Time it?
        if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
        {
          Elemental_assembly_time[e] =
            TimingHelpers::timer() - t_assemble_start;
        }
#endif

      } // End of loop over the elements

    } // list_pt goes out of scope


#ifdef OOMPH_HAS_MPI

    // Postprocess timing information and re-allocate distribution of
    // elements during subsequent assemblies.
    if ((!doing_residuals) && (!Problem_has_been_distributed) &&
        Must_recompute_load_balance_for_assembly)
    {
      recompute_load_balanced_assembly();
    }

    // We have determined load balancing for current setup.
    // This can remain the same until assign_eqn_numbers() is called
    // again -- the flag is re-set to true there.
    if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
    {
      Must_recompute_load_balance_for_assembly = false;
    }

#endif


    //----Finally we need to convert the beautiful list storage scheme---
    //----------to the containers required by SuperLU--------------------

    // Loop over the number of matrices
    for (unsigned m = 0; m < n_matrix; m++)
    {
      // Set the number of rows or columns
      row_or_column_start[m] = new int[ndof + 1];
      // Counter for the total number of entries in the storage scheme
      unsigned long entry_count = 0;
      // The first entry is 0
      row_or_column_start[m][0] = entry_count;

      // first we compute the number of non-zeros
      nnz[m] = 0;
      for (unsigned long i_global = 0; i_global < ndof; i_global++)
      {
        nnz[m] += matrix_data_list[m][i_global].size();
      }

      // and then resize the storage
      column_or_row_index[m] = new int[nnz[m]];
      value[m] = new double[nnz[m]];

      // Now we merely loop over the number of rows or columns
      for (unsigned long i_global = 0; i_global < ndof; i_global++)
      {
        // Start index for the present row is the number of entries so far
        row_or_column_start[m][i_global] = entry_count;
        // If there are no entries in the list then skip the loop
        if (matrix_data_list[m][i_global].empty())
        {
          continue;
        }

        // Sort the list corresponding to this row or column by the
        // column or row index (first entry in the pair).
        // This might be inefficient, but we only have to do the sort ONCE
        // for each list. This is faster than using a map storage scheme, where
        // we are sorting for every insertion (although the map structure
        // is cleaner and more memory efficient)
        matrix_data_list[m][i_global].sort();

        // Set up an iterator for start of the list
        std::list<std::pair<unsigned, double>>::iterator it =
          matrix_data_list[m][i_global].begin();

        // Get the first row or column index in the list...
        unsigned current_index = it->first;
        //...and the corresponding value
        double current_value = it->second;

        // Loop over all the entries in the sorted list
        // Increase the iterator so that we start at the second entry
        for (++it; it != matrix_data_list[m][i_global].end(); ++it)
        {
          // If the index has not changed, then we must add the contribution
          // of the present entry to the value.
          // Additionally check that the entry is non-zero
          if ((it->first == current_index) &&
              (std::fabs(it->second) > Numerical_zero_for_sparse_assembly))
          {
            current_value += it->second;
          }
          // Otherwise, we have added all the contributions to the index
          // to current_value, so add it to the SuperLU data structure
          else
          {
            // Add the row or column index to the vector
            column_or_row_index[m][entry_count] = current_index;
            // Add the actual value to the vector
            value[m][entry_count] = current_value;
            // Increase the counter for the number of entries in each vector
            entry_count++;

            // Set the index and value to be those of the current entry in the
            // list
            current_index = it->first;
            current_value = it->second;
          }
        } // End of loop over all list entries for this global row or column

        // There are TWO special cases to consider.
        // If there is only one equation number in the list, then it
        // will NOT have been added. We test this case by comparing the
        // number of entries with those stored in row_or_column_start[i_global]
        // Otherwise
        // If the final entry in the list has the same index as the penultimate
        // entry, then it will NOT have been added to the SuperLU storage scheme
        // Check this by comparing the current_index with the final index
        // stored in the SuperLU scheme. If they are not the same, then
        // add the current_index and value.

        // If single equation number in list
        if ((static_cast<int>(entry_count) == row_or_column_start[m][i_global])
            // If we have a single equation number, this will not be evaluated.
            // If we don't then we do the test to check that the final
            // entry is added
            || (static_cast<int>(current_index) !=
                column_or_row_index[m][entry_count - 1]))
        {
          // Add the row or column index to the vector
          column_or_row_index[m][entry_count] = current_index;
          // Add the actual value to the vector
          value[m][entry_count] = current_value;
          // Increase the counter for the number of entries in each vector
          entry_count++;
        }

      } // End of loop over the rows or columns of the entire matrix

      // Final entry in the row/column start vector
      row_or_column_start[m][ndof] = entry_count;
    } // End of loop over matrices

    if (Pause_at_end_of_sparse_assembly)
    {
      oomph_info << "Pausing at end of sparse assembly." << std::endl;
      pause("Check memory usage now.");
    }
  }


  //=====================================================================
  /// This is a (private) helper function that is used to assemble system
  /// matrices in compressed row or column format
  /// and compute residual vectors using vectors of pairs
  /// The default action is to assemble the jacobian matrix and
  /// residuals for the Newton method. The action can be
  /// overloaded at an elemental level by chaging the default
  /// behaviour of the function Element::get_all_vectors_and_matrices().
  /// column_or_row_index: Column [or row] index of given entry
  /// row_or_column_start: Index of first entry for given row [or column]
  /// value              : Vector of nonzero entries
  /// residuals          : Residual vector
  /// compressed_row_flag: Bool flag to indicate if storage format is
  ///                      compressed row [if false interpretation of
  ///                      arguments is as stated in square brackets].
  //=====================================================================
  void Problem::sparse_assemble_row_or_column_compressed_with_vectors_of_pairs(
    Vector<int*>& column_or_row_index,
    Vector<int*>& row_or_column_start,
    Vector<double*>& value,
    Vector<unsigned>& nnz,
    Vector<double*>& residuals,
    bool compressed_row_flag)
  {
    // Total number of elements
    const unsigned long n_elements = mesh_pt()->nelement();

    // Default range of elements for distributed problems
    unsigned long el_lo = 0;
    unsigned long el_hi = n_elements - 1;

#ifdef OOMPH_HAS_MPI
    // Otherwise just loop over a fraction of the elements
    // (This will either have been initialised in
    // Problem::set_default_first_and_last_element_for_assembly() or
    // will have been re-assigned during a previous assembly loop
    // Note that following the re-assignment only the entries
    // for the current processor are relevant.
    if (!Problem_has_been_distributed)
    {
      el_lo = First_el_for_assembly[Communicator_pt->my_rank()];
      el_hi = Last_el_plus_one_for_assembly[Communicator_pt->my_rank()] - 1;
    }
#endif

    // number of local eqns
    unsigned ndof = this->ndof();

    // Find the number of vectors to be assembled
    const unsigned n_vector = residuals.size();

    // Find the number of matrices to be assembled
    const unsigned n_matrix = column_or_row_index.size();

    // Locally cache pointer to assembly handler
    AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

#ifdef OOMPH_HAS_MPI
    bool doing_residuals = false;
    if (dynamic_cast<ParallelResidualsHandler*>(Assembly_handler_pt) != 0)
    {
      doing_residuals = true;
    }
#endif

// Error check dimensions
#ifdef PARANOID
    if (row_or_column_start.size() != n_matrix)
    {
      std::ostringstream error_stream;
      error_stream << "Error: " << std::endl
                   << "row_or_column_start.size() "
                   << row_or_column_start.size() << " does not equal "
                   << "column_or_row_index.size() "
                   << column_or_row_index.size() << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    if (value.size() != n_matrix)
    {
      std::ostringstream error_stream;
      error_stream << "Error: " << std::endl
                   << "value.size() " << value.size() << " does not equal "
                   << "column_or_row_index.size() "
                   << column_or_row_index.size() << std::endl
                   << std::endl
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // The idea behind this sparse assembly routine is to use a Vector of
    // Vectors of pairs for each complete matrix.
    // Each inner Vector stores pairs and holds the row (or column) index
    // and the value of the matrix entry.

    // Set up Vector of Vectors to store the entries of each matrix,
    // indexed by either the column or row.
    Vector<Vector<Vector<std::pair<unsigned, double>>>> matrix_data(n_matrix);

    // Loop over the number of matrices being assembled and resize
    // each Vector of Vectors to the number of rows or columns of the matrix
    for (unsigned m = 0; m < n_matrix; m++)
    {
      matrix_data[m].resize(ndof);
    }

    // Resize the residuals vectors
    for (unsigned v = 0; v < n_vector; v++)
    {
      residuals[v] = new double[ndof];
      for (unsigned i = 0; i < ndof; i++)
      {
        residuals[v][i] = 0;
      }
    }

#ifdef OOMPH_HAS_MPI

    // Storage for assembly time for elements
    double t_assemble_start = 0.0;

    // Storage for assembly times
    if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
    {
      Elemental_assembly_time.resize(n_elements);
    }

#endif

    //----------------Assemble and populate the vector storage scheme--------
    {
      // Allocate local storage for the element's contribution to the
      // residuals vectors and system matrices of the size of the maximum
      // number of dofs in any element
      // This means that the storage is only allocated (and deleted) once
      Vector<Vector<double>> el_residuals(n_vector);
      Vector<DenseMatrix<double>> el_jacobian(n_matrix);

      // Loop over the elements
      for (unsigned long e = el_lo; e <= el_hi; e++)
      {
#ifdef OOMPH_HAS_MPI
        // Time it?
        if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
        {
          t_assemble_start = TimingHelpers::timer();
        }
#endif

        // Get the pointer to the element
        GeneralisedElement* elem_pt = mesh_pt()->element_pt(e);

#ifdef OOMPH_HAS_MPI
        // Ignore halo elements
        if (!elem_pt->is_halo())
        {
#endif

          // Find number of degrees of freedom in the element
          const unsigned nvar = assembly_handler_pt->ndof(elem_pt);

          // Resize the storage for elemental jacobian and residuals
          for (unsigned v = 0; v < n_vector; v++)
          {
            el_residuals[v].resize(nvar);
          }
          for (unsigned m = 0; m < n_matrix; m++)
          {
            el_jacobian[m].resize(nvar);
          }

          // Now get the residuals and jacobian for the element
          assembly_handler_pt->get_all_vectors_and_matrices(
            elem_pt, el_residuals, el_jacobian);

          //---------------Insert the values into the vectors--------------

          // Loop over the first index of local variables
          for (unsigned i = 0; i < nvar; i++)
          {
            // Get the local equation number
            unsigned eqn_number = assembly_handler_pt->eqn_number(elem_pt, i);

            // Add the contribution to the residuals
            for (unsigned v = 0; v < n_vector; v++)
            {
              // Fill in each residuals vector
              residuals[v][eqn_number] += el_residuals[v][i];
            }

            // Now loop over the other index
            for (unsigned j = 0; j < nvar; j++)
            {
              // Get the number of the unknown
              unsigned unknown = assembly_handler_pt->eqn_number(elem_pt, j);

              // Loop over the matrices
              // If it's compressed row storage, then our vector of maps
              // is indexed by row (equation number)
              for (unsigned m = 0; m < n_matrix; m++)
              {
                // Get the value of the matrix at this point
                double value = el_jacobian[m](i, j);
                // Only bother to add to the vector if it's non-zero
                if (std::fabs(value) > Numerical_zero_for_sparse_assembly)
                {
                  // If it's compressed row storage, then our vector of maps
                  // is indexed by row (equation number)
                  if (compressed_row_flag)
                  {
                    // Find the correct position and add the data into the
                    // vectors
                    const unsigned size = matrix_data[m][eqn_number].size();
                    for (unsigned k = 0; k <= size; k++)
                    {
                      if (k == size)
                      {
                        matrix_data[m][eqn_number].push_back(
                          std::make_pair(unknown, value));
                        break;
                      }
                      else if (matrix_data[m][eqn_number][k].first == unknown)
                      {
                        matrix_data[m][eqn_number][k].second += value;
                        break;
                      }
                    }
                  }
                  // Otherwise it's compressed column storage and our vector is
                  // indexed by column (the unknown)
                  else
                  {
                    // Add the data into the vectors in the correct position
                    const unsigned size = matrix_data[m][unknown].size();
                    for (unsigned k = 0; k <= size; k++)
                    {
                      if (k == size)
                      {
                        matrix_data[m][unknown].push_back(
                          std::make_pair(eqn_number, value));
                        break;
                      }
                      else if (matrix_data[m][unknown][k].first == eqn_number)
                      {
                        matrix_data[m][unknown][k].second += value;
                        break;
                      }
                    }
                  }
                }
              } // End of loop over matrices
            }
          }

#ifdef OOMPH_HAS_MPI
        } // endif halo element
#endif


#ifdef OOMPH_HAS_MPI
        // Time it?
        if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
        {
          Elemental_assembly_time[e] =
            TimingHelpers::timer() - t_assemble_start;
        }
#endif

      } // End of loop over the elements


    } // End of vector assembly


#ifdef OOMPH_HAS_MPI

    // Postprocess timing information and re-allocate distribution of
    // elements during subsequent assemblies.
    if ((!doing_residuals) && (!Problem_has_been_distributed) &&
        Must_recompute_load_balance_for_assembly)
    {
      recompute_load_balanced_assembly();
    }

    // We have determined load balancing for current setup.
    // This can remain the same until assign_eqn_numbers() is called
    // again -- the flag is re-set to true there.
    if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
    {
      Must_recompute_load_balance_for_assembly = false;
    }

#endif


    //-----------Finally we need to convert this vector storage scheme
    //------------------------to the containers required by SuperLU

    // Loop over the number of matrices
    for (unsigned m = 0; m < n_matrix; m++)
    {
      // Set the number of rows or columns
      row_or_column_start[m] = new int[ndof + 1];

      // fill row_or_column_start and find the number of entries
      row_or_column_start[m][0] = 0;
      for (unsigned long i = 0; i < ndof; i++)
      {
        row_or_column_start[m][i + 1] =
          row_or_column_start[m][i] + matrix_data[m][i].size();
      }
      const unsigned entries = row_or_column_start[m][ndof];

      // resize vectors
      column_or_row_index[m] = new int[entries];
      value[m] = new double[entries];
      nnz[m] = entries;

      // Now we merely loop over the number of rows or columns
      for (unsigned long i_global = 0; i_global < ndof; i_global++)
      {
        // If there are no entries in the vector then skip the rest of the loop
        if (matrix_data[m][i_global].empty())
        {
          continue;
        }

        // Loop over all the entries in the vectors corresponding to the given
        // row or column. It will NOT be ordered
        unsigned p = 0;
        for (int j = row_or_column_start[m][i_global];
             j < row_or_column_start[m][i_global + 1];
             j++)
        {
          column_or_row_index[m][j] = matrix_data[m][i_global][p].first;
          value[m][j] = matrix_data[m][i_global][p].second;
          ++p;
        }
      }
    } // End of the loop over the matrices

    if (Pause_at_end_of_sparse_assembly)
    {
      oomph_info << "Pausing at end of sparse assembly." << std::endl;
      pause("Check memory usage now.");
    }
  }


  //=====================================================================
  /// This is a (private) helper function that is used to assemble system
  /// matrices in compressed row or column format
  /// and compute residual vectors using two vectors.
  /// The default action is to assemble the jacobian matrix and
  /// residuals for the Newton method. The action can be
  /// overloaded at an elemental level by chaging the default
  /// behaviour of the function Element::get_all_vectors_and_matrices().
  /// column_or_row_index: Column [or row] index of given entry
  /// row_or_column_start: Index of first entry for given row [or column]
  /// value              : Vector of nonzero entries
  /// residuals          : Residual vector
  /// compressed_row_flag: Bool flag to indicate if storage format is
  ///                      compressed row [if false interpretation of
  ///                      arguments is as stated in square brackets].
  //=====================================================================
  void Problem::sparse_assemble_row_or_column_compressed_with_two_vectors(
    Vector<int*>& column_or_row_index,
    Vector<int*>& row_or_column_start,
    Vector<double*>& value,
    Vector<unsigned>& nnz,
    Vector<double*>& residuals,
    bool compressed_row_flag)
  {
    // Total number of elements
    const unsigned long n_elements = mesh_pt()->nelement();

    // Default range of elements for distributed problems
    unsigned long el_lo = 0;
    unsigned long el_hi = n_elements - 1;


#ifdef OOMPH_HAS_MPI
    // Otherwise just loop over a fraction of the elements
    // (This will either have been initialised in
    // Problem::set_default_first_and_last_element_for_assembly() or
    // will have been re-assigned during a previous assembly loop
    // Note that following the re-assignment only the entries
    // for the current processor are relevant.
    if (!Problem_has_been_distributed)
    {
      el_lo = First_el_for_assembly[Communicator_pt->my_rank()];
      el_hi = Last_el_plus_one_for_assembly[Communicator_pt->my_rank()] - 1;
    }
#endif

    // number of local eqns
    unsigned ndof = this->ndof();

    // Find the number of vectors to be assembled
    const unsigned n_vector = residuals.size();

    // Find the number of matrices to be assembled
    const unsigned n_matrix = column_or_row_index.size();

    // Locally cache pointer to assembly handler
    AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

#ifdef OOMPH_HAS_MPI
    bool doing_residuals = false;
    if (dynamic_cast<ParallelResidualsHandler*>(Assembly_handler_pt) != 0)
    {
      doing_residuals = true;
    }
#endif

// Error check dimensions
#ifdef PARANOID
    if (row_or_column_start.size() != n_matrix)
    {
      std::ostringstream error_stream;
      error_stream << "Error: " << std::endl
                   << "row_or_column_start.size() "
                   << row_or_column_start.size() << " does not equal "
                   << "column_or_row_index.size() "
                   << column_or_row_index.size() << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    if (value.size() != n_matrix)
    {
      std::ostringstream error_stream;
      error_stream << "Error: " << std::endl
                   << "value.size() " << value.size() << " does not equal "
                   << "column_or_row_index.size() "
                   << column_or_row_index.size() << std::endl
                   << std::endl
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // The idea behind this sparse assembly routine is to use Vectors of
    // Vectors for the entries in each complete matrix. And a second
    // Vector of Vectors stores the global row (or column) indeces. This
    // will not have the memory overheads associated with the methods using
    // lists or maps, but insertion will be more costly.

    // Set up two vector of vectors to store the entries of each  matrix,
    // indexed by either the column or row. The entries of the vector for
    // each matrix correspond to all the rows or columns of that matrix.
    Vector<Vector<Vector<unsigned>>> matrix_row_or_col_indices(n_matrix);
    Vector<Vector<Vector<double>>> matrix_values(n_matrix);

    // Loop over the number of matrices being assembled and resize
    // each vector of vectors to the number of rows or columns of the matrix
    for (unsigned m = 0; m < n_matrix; m++)
    {
      matrix_row_or_col_indices[m].resize(ndof);
      matrix_values[m].resize(ndof);
    }

    // Resize the residuals vectors
    for (unsigned v = 0; v < n_vector; v++)
    {
      residuals[v] = new double[ndof];
      for (unsigned i = 0; i < ndof; i++)
      {
        residuals[v][i] = 0;
      }
    }

#ifdef OOMPH_HAS_MPI

    // Storage for assembly time for elements
    double t_assemble_start = 0.0;

    // Storage for assembly times
    if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
    {
      Elemental_assembly_time.resize(n_elements);
    }

#endif


    //----------------Assemble and populate the vector storage scheme-------
    {
      // Allocate local storage for the element's contribution to the
      // residuals vectors and system matrices of the size of the maximum
      // number of dofs in any element
      // This means that the storage will only be allocated (and deleted) once
      Vector<Vector<double>> el_residuals(n_vector);
      Vector<DenseMatrix<double>> el_jacobian(n_matrix);

      // Loop over the elements
      for (unsigned long e = el_lo; e <= el_hi; e++)
      {
#ifdef OOMPH_HAS_MPI
        // Time it?
        if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
        {
          t_assemble_start = TimingHelpers::timer();
        }
#endif

        // Get the pointer to the element
        GeneralisedElement* elem_pt = mesh_pt()->element_pt(e);

#ifdef OOMPH_HAS_MPI
        // Ignore halo elements
        if (!elem_pt->is_halo())
        {
#endif

          // Find number of degrees of freedom in the element
          const unsigned nvar = assembly_handler_pt->ndof(elem_pt);

          // Resize the storage for elemental jacobian and residuals
          for (unsigned v = 0; v < n_vector; v++)
          {
            el_residuals[v].resize(nvar);
          }
          for (unsigned m = 0; m < n_matrix; m++)
          {
            el_jacobian[m].resize(nvar);
          }

          // Now get the residuals and jacobian for the element
          assembly_handler_pt->get_all_vectors_and_matrices(
            elem_pt, el_residuals, el_jacobian);

          //---------------Insert the values into the vectors--------------

          // Loop over the first index of local variables
          for (unsigned i = 0; i < nvar; i++)
          {
            // Get the local equation number
            unsigned eqn_number = assembly_handler_pt->eqn_number(elem_pt, i);

            // Add the contribution to the residuals
            for (unsigned v = 0; v < n_vector; v++)
            {
              // Fill in each residuals vector
              residuals[v][eqn_number] += el_residuals[v][i];
            }

            // Now loop over the other index
            for (unsigned j = 0; j < nvar; j++)
            {
              // Get the number of the unknown
              unsigned unknown = assembly_handler_pt->eqn_number(elem_pt, j);

              // Loop over the matrices
              // If it's compressed row storage, then our vector of maps
              // is indexed by row (equation number)
              for (unsigned m = 0; m < n_matrix; m++)
              {
                // Get the value of the matrix at this point
                double value = el_jacobian[m](i, j);
                // Only bother to add to the vector if it's non-zero
                if (std::fabs(value) > Numerical_zero_for_sparse_assembly)
                {
                  // If it's compressed row storage, then our vector of maps
                  // is indexed by row (equation number)
                  if (compressed_row_flag)
                  {
                    // Find the correct position and add the data into the
                    // vectors
                    const unsigned size =
                      matrix_row_or_col_indices[m][eqn_number].size();

                    for (unsigned k = 0; k <= size; k++)
                    {
                      if (k == size)
                      {
                        matrix_row_or_col_indices[m][eqn_number].push_back(
                          unknown);
                        matrix_values[m][eqn_number].push_back(value);
                        break;
                      }
                      else if (matrix_row_or_col_indices[m][eqn_number][k] ==
                               unknown)
                      {
                        matrix_values[m][eqn_number][k] += value;
                        break;
                      }
                    }
                  }
                  // Otherwise it's compressed column storage and our vector is
                  // indexed by column (the unknown)
                  else
                  {
                    // Add the data into the vectors in the correct position
                    const unsigned size =
                      matrix_row_or_col_indices[m][unknown].size();
                    for (unsigned k = 0; k <= size; k++)
                    {
                      if (k == size)
                      {
                        matrix_row_or_col_indices[m][unknown].push_back(
                          eqn_number);
                        matrix_values[m][unknown].push_back(value);
                        break;
                      }
                      else if (matrix_row_or_col_indices[m][unknown][k] ==
                               eqn_number)
                      {
                        matrix_values[m][unknown][k] += value;
                        break;
                      }
                    }
                  }
                }
              } // End of loop over matrices
            }
          }

#ifdef OOMPH_HAS_MPI
        } // endif halo element
#endif


#ifdef OOMPH_HAS_MPI
        // Time it?
        if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
        {
          Elemental_assembly_time[e] =
            TimingHelpers::timer() - t_assemble_start;
        }
#endif

      } // End of loop over the elements

    } // End of vector assembly


#ifdef OOMPH_HAS_MPI

    // Postprocess timing information and re-allocate distribution of
    // elements during subsequent assemblies.
    if ((!doing_residuals) && (!Problem_has_been_distributed) &&
        Must_recompute_load_balance_for_assembly)
    {
      recompute_load_balanced_assembly();
    }

    // We have determined load balancing for current setup.
    // This can remain the same until assign_eqn_numbers() is called
    // again -- the flag is re-set to true there.
    if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
    {
      Must_recompute_load_balance_for_assembly = false;
    }

#endif

    //-----------Finally we need to convert this lousy vector storage scheme
    //------------------------to the containers required by SuperLU

    // Loop over the number of matrices
    for (unsigned m = 0; m < n_matrix; m++)
    {
      // Set the number of rows or columns
      row_or_column_start[m] = new int[ndof + 1];

      // fill row_or_column_start and find the number of entries
      row_or_column_start[m][0] = 0;
      for (unsigned long i = 0; i < ndof; i++)
      {
        row_or_column_start[m][i + 1] =
          row_or_column_start[m][i] + matrix_values[m][i].size();
      }
      const unsigned entries = row_or_column_start[m][ndof];

      // resize vectors
      column_or_row_index[m] = new int[entries];
      value[m] = new double[entries];
      nnz[m] = entries;

      // Now we merely loop over the number of rows or columns
      for (unsigned long i_global = 0; i_global < ndof; i_global++)
      {
        // If there are no entries in the vector then skip the rest of the loop
        if (matrix_values[m][i_global].empty())
        {
          continue;
        }

        // Loop over all the entries in the vectors corresponding to the given
        // row or column. It will NOT be ordered
        unsigned p = 0;
        for (int j = row_or_column_start[m][i_global];
             j < row_or_column_start[m][i_global + 1];
             j++)
        {
          column_or_row_index[m][j] = matrix_row_or_col_indices[m][i_global][p];
          value[m][j] = matrix_values[m][i_global][p];
          ++p;
        }
      }
    } // End of the loop over the matrices

    if (Pause_at_end_of_sparse_assembly)
    {
      oomph_info << "Pausing at end of sparse assembly." << std::endl;
      pause("Check memory usage now.");
    }
  }


  //=====================================================================
  /// This is a (private) helper function that is used to assemble system
  /// matrices in compressed row or column format
  /// and compute residual vectors using two vectors.
  /// The default action is to assemble the jacobian matrix and
  /// residuals for the Newton method. The action can be
  /// overloaded at an elemental level by chaging the default
  /// behaviour of the function Element::get_all_vectors_and_matrices().
  /// column_or_row_index: Column [or row] index of given entry
  /// row_or_column_start: Index of first entry for given row [or column]
  /// value              : Vector of nonzero entries
  /// residuals          : Residual vector
  /// compressed_row_flag: Bool flag to indicate if storage format is
  ///                      compressed row [if false interpretation of
  ///                      arguments is as stated in square brackets].
  //=====================================================================
  void Problem::sparse_assemble_row_or_column_compressed_with_two_arrays(
    Vector<int*>& column_or_row_index,
    Vector<int*>& row_or_column_start,
    Vector<double*>& value,
    Vector<unsigned>& nnz,
    Vector<double*>& residuals,
    bool compressed_row_flag)
  {
    // Total number of elements
    const unsigned long n_elements = mesh_pt()->nelement();

    // Default range of elements for distributed problems
    unsigned long el_lo = 0;
    unsigned long el_hi = n_elements - 1;


#ifdef OOMPH_HAS_MPI
    // Otherwise just loop over a fraction of the elements
    // (This will either have been initialised in
    // Problem::set_default_first_and_last_element_for_assembly() or
    // will have been re-assigned during a previous assembly loop
    // Note that following the re-assignment only the entries
    // for the current processor are relevant.
    if (!Problem_has_been_distributed)
    {
      el_lo = First_el_for_assembly[Communicator_pt->my_rank()];
      el_hi = Last_el_plus_one_for_assembly[Communicator_pt->my_rank()] - 1;
    }
#endif

    // number of local eqns
    unsigned ndof = this->ndof();

    // Find the number of vectors to be assembled
    const unsigned n_vector = residuals.size();

    // Find the number of matrices to be assembled
    const unsigned n_matrix = column_or_row_index.size();

    // Locally cache pointer to assembly handler
    AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

#ifdef OOMPH_HAS_MPI
    bool doing_residuals = false;
    if (dynamic_cast<ParallelResidualsHandler*>(Assembly_handler_pt) != 0)
    {
      doing_residuals = true;
    }
#endif

// Error check dimensions
#ifdef PARANOID
    if (row_or_column_start.size() != n_matrix)
    {
      std::ostringstream error_stream;
      error_stream << "Error: " << std::endl
                   << "row_or_column_start.size() "
                   << row_or_column_start.size() << " does not equal "
                   << "column_or_row_index.size() "
                   << column_or_row_index.size() << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    if (value.size() != n_matrix)
    {
      std::ostringstream error_stream;
      error_stream << "Error: " << std::endl
                   << "value.size() " << value.size() << " does not equal "
                   << "column_or_row_index.size() "
                   << column_or_row_index.size() << std::endl
                   << std::endl
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // The idea behind this sparse assembly routine is to use Vectors of
    // Vectors for the entries in each complete matrix. And a second
    // Vector of Vectors stores the global row (or column) indeces. This
    // will not have the memory overheads associated with the methods using
    // lists or maps, but insertion will be more costly.

    // Set up two vector of vectors to store the entries of each  matrix,
    // indexed by either the column or row. The entries of the vector for
    // each matrix correspond to all the rows or columns of that matrix.
    Vector<unsigned**> matrix_row_or_col_indices(n_matrix);
    Vector<double**> matrix_values(n_matrix);

    // Loop over the number of matrices being assembled and resize
    // each vector of vectors to the number of rows or columns of the matrix
    for (unsigned m = 0; m < n_matrix; m++)
    {
      matrix_row_or_col_indices[m] = new unsigned*[ndof];
      matrix_values[m] = new double*[ndof];
    }

    // Resize the residuals vectors
    for (unsigned v = 0; v < n_vector; v++)
    {
      residuals[v] = new double[ndof];
      for (unsigned i = 0; i < ndof; i++)
      {
        residuals[v][i] = 0;
      }
    }

#ifdef OOMPH_HAS_MPI

    // Storage for assembly time for elements
    double t_assemble_start = 0.0;

    // Storage for assembly times
    if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
    {
      Elemental_assembly_time.resize(n_elements);
    }

#endif

    // number of coefficients in each row
    Vector<Vector<unsigned>> ncoef(n_matrix);
    for (unsigned m = 0; m < n_matrix; m++)
    {
      ncoef[m].resize(ndof, 0);
    }

    if (Sparse_assemble_with_arrays_previous_allocation.size() == 0)
    {
      Sparse_assemble_with_arrays_previous_allocation.resize(n_matrix);
      for (unsigned m = 0; m < n_matrix; m++)
      {
        Sparse_assemble_with_arrays_previous_allocation[m].resize(ndof, 0);
      }
    }

    //----------------Assemble and populate the vector storage scheme-------
    {
      // Allocate local storage for the element's contribution to the
      // residuals vectors and system matrices of the size of the maximum
      // number of dofs in any element
      // This means that the storage will only be allocated (and deleted) once
      Vector<Vector<double>> el_residuals(n_vector);
      Vector<DenseMatrix<double>> el_jacobian(n_matrix);

      // Loop over the elements
      for (unsigned long e = el_lo; e <= el_hi; e++)
      {
#ifdef OOMPH_HAS_MPI
        // Time it?
        if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
        {
          t_assemble_start = TimingHelpers::timer();
        }
#endif

        // Get the pointer to the element
        GeneralisedElement* elem_pt = mesh_pt()->element_pt(e);

#ifdef OOMPH_HAS_MPI
        // Ignore halo elements
        if (!elem_pt->is_halo())
        {
#endif

          // Find number of degrees of freedom in the element
          const unsigned nvar = assembly_handler_pt->ndof(elem_pt);

          // Resize the storage for elemental jacobian and residuals
          for (unsigned v = 0; v < n_vector; v++)
          {
            el_residuals[v].resize(nvar);
          }
          for (unsigned m = 0; m < n_matrix; m++)
          {
            el_jacobian[m].resize(nvar);
          }

          // Now get the residuals and jacobian for the element
          assembly_handler_pt->get_all_vectors_and_matrices(
            elem_pt, el_residuals, el_jacobian);

          //---------------Insert the values into the vectors--------------

          // Loop over the first index of local variables
          for (unsigned i = 0; i < nvar; i++)
          {
            // Get the local equation number
            unsigned eqn_number = assembly_handler_pt->eqn_number(elem_pt, i);

            // Add the contribution to the residuals
            for (unsigned v = 0; v < n_vector; v++)
            {
              // Fill in each residuals vector
              residuals[v][eqn_number] += el_residuals[v][i];
            }

            // Now loop over the other index
            for (unsigned j = 0; j < nvar; j++)
            {
              // Get the number of the unknown
              unsigned unknown = assembly_handler_pt->eqn_number(elem_pt, j);

              // Loop over the matrices
              // If it's compressed row storage, then our vector of maps
              // is indexed by row (equation number)
              for (unsigned m = 0; m < n_matrix; m++)
              {
                // Get the value of the matrix at this point
                double value = el_jacobian[m](i, j);
                // Only bother to add to the vector if it's non-zero
                if (std::fabs(value) > Numerical_zero_for_sparse_assembly)
                {
                  // number of entrys in this row
                  const unsigned size = ncoef[m][eqn_number];

                  // if no data has been allocated for this row then allocate
                  if (size == 0)
                  {
                    // do we have previous allocation data
                    if (Sparse_assemble_with_arrays_previous_allocation
                          [m][eqn_number] != 0)
                    {
                      matrix_row_or_col_indices[m][eqn_number] = new unsigned
                        [Sparse_assemble_with_arrays_previous_allocation
                           [m][eqn_number]];
                      matrix_values[m][eqn_number] = new double
                        [Sparse_assemble_with_arrays_previous_allocation
                           [m][eqn_number]];
                    }
                    else
                    {
                      matrix_row_or_col_indices[m][eqn_number] = new unsigned
                        [Sparse_assemble_with_arrays_initial_allocation];
                      matrix_values[m][eqn_number] = new double
                        [Sparse_assemble_with_arrays_initial_allocation];
                      Sparse_assemble_with_arrays_previous_allocation
                        [m][eqn_number] =
                          Sparse_assemble_with_arrays_initial_allocation;
                    }
                  }

                  // If it's compressed row storage, then our vector of maps
                  // is indexed by row (equation number)
                  if (compressed_row_flag)
                  {
                    // next add the data
                    for (unsigned k = 0; k <= size; k++)
                    {
                      if (k == size)
                      {
                        // do we need to allocate more storage
                        if (Sparse_assemble_with_arrays_previous_allocation
                              [m][eqn_number] == ncoef[m][eqn_number])
                        {
                          unsigned new_allocation =
                            ncoef[m][eqn_number] +
                            Sparse_assemble_with_arrays_allocation_increment;
                          double* new_values = new double[new_allocation];
                          unsigned* new_indices = new unsigned[new_allocation];
                          for (unsigned c = 0; c < ncoef[m][eqn_number]; c++)
                          {
                            new_values[c] = matrix_values[m][eqn_number][c];
                            new_indices[c] =
                              matrix_row_or_col_indices[m][eqn_number][c];
                          }
                          delete[] matrix_values[m][eqn_number];
                          delete[] matrix_row_or_col_indices[m][eqn_number];
                          matrix_values[m][eqn_number] = new_values;
                          matrix_row_or_col_indices[m][eqn_number] =
                            new_indices;
                          Sparse_assemble_with_arrays_previous_allocation
                            [m][eqn_number] = new_allocation;
                        }
                        // and now add the data
                        unsigned entry = ncoef[m][eqn_number];
                        ncoef[m][eqn_number]++;
                        matrix_row_or_col_indices[m][eqn_number][entry] =
                          unknown;
                        matrix_values[m][eqn_number][entry] = value;
                        break;
                      }
                      else if (matrix_row_or_col_indices[m][eqn_number][k] ==
                               unknown)
                      {
                        matrix_values[m][eqn_number][k] += value;
                        break;
                      }
                    }
                  }
                  // Otherwise it's compressed column storage and our vector is
                  // indexed by column (the unknown)
                  else
                  {
                    // Add the data into the vectors in the correct position
                    for (unsigned k = 0; k <= size; k++)
                    {
                      if (k == size)
                      {
                        // do we need to allocate more storage
                        if (Sparse_assemble_with_arrays_previous_allocation
                              [m][unknown] == ncoef[m][unknown])
                        {
                          unsigned new_allocation =
                            ncoef[m][unknown] +
                            Sparse_assemble_with_arrays_allocation_increment;
                          double* new_values = new double[new_allocation];
                          unsigned* new_indices = new unsigned[new_allocation];
                          for (unsigned c = 0; c < ncoef[m][unknown]; c++)
                          {
                            new_values[c] = matrix_values[m][unknown][c];
                            new_indices[c] =
                              matrix_row_or_col_indices[m][unknown][c];
                          }
                          delete[] matrix_values[m][unknown];
                          delete[] matrix_row_or_col_indices[m][unknown];
                          Sparse_assemble_with_arrays_previous_allocation
                            [m][unknown] = new_allocation;
                        }
                        // and now add the data
                        unsigned entry = ncoef[m][unknown];
                        ncoef[m][unknown]++;
                        matrix_row_or_col_indices[m][unknown][entry] =
                          eqn_number;
                        matrix_values[m][unknown][entry] = value;
                        break;
                      }
                      else if (matrix_row_or_col_indices[m][unknown][k] ==
                               eqn_number)
                      {
                        matrix_values[m][unknown][k] += value;
                        break;
                      }
                    }
                  }
                }
              } // End of loop over matrices
            }
          }

#ifdef OOMPH_HAS_MPI
        } // endif halo element
#endif


#ifdef OOMPH_HAS_MPI
        // Time it?
        if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
        {
          Elemental_assembly_time[e] =
            TimingHelpers::timer() - t_assemble_start;
        }
#endif

      } // End of loop over the elements

    } // End of vector assembly


#ifdef OOMPH_HAS_MPI

    // Postprocess timing information and re-allocate distribution of
    // elements during subsequent assemblies.
    if ((!doing_residuals) && (!Problem_has_been_distributed) &&
        Must_recompute_load_balance_for_assembly)
    {
      recompute_load_balanced_assembly();
    }

    // We have determined load balancing for current setup.
    // This can remain the same until assign_eqn_numbers() is called
    // again -- the flag is re-set to true there.
    if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
    {
      Must_recompute_load_balance_for_assembly = false;
    }

#endif

    //-----------Finally we need to convert this lousy vector storage scheme
    //------------------------to the containers required by SuperLU

    // Loop over the number of matrices
    for (unsigned m = 0; m < n_matrix; m++)
    {
      // Set the number of rows or columns
      row_or_column_start[m] = new int[ndof + 1];

      // fill row_or_column_start and find the number of entries
      row_or_column_start[m][0] = 0;
      for (unsigned long i = 0; i < ndof; i++)
      {
        row_or_column_start[m][i + 1] = row_or_column_start[m][i] + ncoef[m][i];
        Sparse_assemble_with_arrays_previous_allocation[m][i] = ncoef[m][i];
      }
      const unsigned entries = row_or_column_start[m][ndof];

      // resize vectors
      column_or_row_index[m] = new int[entries];
      value[m] = new double[entries];
      nnz[m] = entries;

      // Now we merely loop over the number of rows or columns
      for (unsigned long i_global = 0; i_global < ndof; i_global++)
      {
        // If there are no entries in the vector then skip the rest of the loop
        if (ncoef[m][i_global] == 0)
        {
          continue;
        }

        // Loop over all the entries in the vectors corresponding to the given
        // row or column. It will NOT be ordered
        unsigned p = 0;
        for (int j = row_or_column_start[m][i_global];
             j < row_or_column_start[m][i_global + 1];
             j++)
        {
          column_or_row_index[m][j] = matrix_row_or_col_indices[m][i_global][p];
          value[m][j] = matrix_values[m][i_global][p];
          ++p;
        }

        // and delete
        delete[] matrix_row_or_col_indices[m][i_global];
        delete[] matrix_values[m][i_global];
      }

      //
      delete[] matrix_row_or_col_indices[m];
      delete[] matrix_values[m];
    } // End of the loop over the matrices

    if (Pause_at_end_of_sparse_assembly)
    {
      oomph_info << "Pausing at end of sparse assembly." << std::endl;
      pause("Check memory usage now.");
    }
  }


#ifdef OOMPH_HAS_MPI
  //=======================================================================
  ///\short Helper method that returns the global equations to which
  /// the elements in the range el_lo to el_hi contribute on this
  /// processor
  //=======================================================================
  void Problem::get_my_eqns(AssemblyHandler* const& assembly_handler_pt,
                            const unsigned& el_lo,
                            const unsigned& el_hi,
                            Vector<unsigned>& my_eqns)
  {
    // Index to keep track of the equations counted
    unsigned my_eqns_index = 0;

    // Loop over the selection of elements
    for (unsigned long e = el_lo; e <= el_hi; e++)
    {
      // Get the pointer to the element
      GeneralisedElement* elem_pt = this->mesh_pt()->element_pt(e);

      // Ignore halo elements
      if (!elem_pt->is_halo())
      {
        // Find number of degrees of freedom in the element
        const unsigned nvar = assembly_handler_pt->ndof(elem_pt);
        // Add the number of dofs to the current size of my_eqns
        my_eqns.resize(my_eqns_index + nvar);

        // Loop over the first index of local variables
        for (unsigned i = 0; i < nvar; i++)
        {
          // Get the local equation number
          unsigned global_eqn_number =
            assembly_handler_pt->eqn_number(elem_pt, i);
          // Add into the vector
          my_eqns[my_eqns_index + i] = global_eqn_number;
        }
        // Update the number of elements in the vector
        my_eqns_index += nvar;
      }
    }

    //  now sort and remove duplicate entries in the vector
    std::sort(my_eqns.begin(), my_eqns.end());
    Vector<unsigned>::iterator it = std::unique(my_eqns.begin(), my_eqns.end());
    my_eqns.resize(it - my_eqns.begin());
  }


  //=============================================================================
  /// \short Helper method to assemble CRDoubleMatrices from distributed
  /// on multiple processors.
  //=============================================================================
  void Problem::parallel_sparse_assemble(
    const LinearAlgebraDistribution* const& target_dist_pt,
    Vector<int*>& column_indices,
    Vector<int*>& row_start,
    Vector<double*>& values,
    Vector<unsigned>& nnz,
    Vector<double*>& residuals)
  {
    // Time assembly
    double t_start = TimingHelpers::timer();

    // my rank and nproc
    unsigned my_rank = Communicator_pt->my_rank();
    unsigned nproc = Communicator_pt->nproc();

    // Total number of elements
    const unsigned long n_elements = mesh_pt()->nelement();

#ifdef PARANOID
    // No elements? This is usually a sign that the problem distribution has
    // led to one processor not having any elements. Either
    // a sign of something having gone wrong or a relatively small
    // problem on a huge number of processors
    if (n_elements == 0)
    {
      std::ostringstream error_stream;
      error_stream << "Processsor " << my_rank << " has no elements. \n"
                   << "This is usually a sign that the problem distribution \n"
                   << "or the load balancing have gone wrong.";
      OomphLibWarning(error_stream.str(),
                      "Problem::parallel_sparse_assemble()",
                      OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // Default range of elements for distributed problems.
    unsigned long el_lo = 0;
    unsigned long el_hi_plus_one = n_elements;

    // Otherwise just loop over a fraction of the elements
    // (This will either have been initialised in
    // Problem::set_default_first_and_last_element_for_assembly() or
    // will have been re-assigned during a previous assembly loop
    // Note that following the re-assignment only the entries
    // for the current processor are relevant.
    if (!Problem_has_been_distributed)
    {
      el_lo = First_el_for_assembly[my_rank];
      el_hi_plus_one = Last_el_plus_one_for_assembly[my_rank];
    }

    // Find the number of vectors to be assembled
    const unsigned n_vector = residuals.size();

    // Find the number of matrices to be assembled
    const unsigned n_matrix = column_indices.size();

    // Locally cache pointer to assembly handler
    AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

    bool doing_residuals = false;
    if (dynamic_cast<ParallelResidualsHandler*>(Assembly_handler_pt) != 0)
    {
      doing_residuals = true;
    }

// Error check dimensions
#ifdef PARANOID
    if (row_start.size() != n_matrix)
    {
      std::ostringstream error_stream;
      error_stream << "Error: " << std::endl
                   << "row_or_column_start.size() " << row_start.size()
                   << " does not equal "
                   << "column_or_row_index.size() " << column_indices.size()
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    if (values.size() != n_matrix)
    {
      std::ostringstream error_stream;
      error_stream << "Error: " << std::endl
                   << "value.size() " << values.size() << " does not equal "
                   << "column_or_row_index.size() " << column_indices.size()
                   << std::endl
                   << std::endl
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // start by assembling the sorted set of equations to which this processor
    // contributes. Essentially this is every global equation that features in
    // all the non-halo elements. This may not be the same as the locally-stored
    // dofs because some of the Nodes in non-halo elements may actually
    // be halos.
    //======================================================================
    Vector<unsigned> my_eqns;
    if (n_elements != 0)
    {
      this->get_my_eqns(
        assembly_handler_pt, el_lo, el_hi_plus_one - 1, my_eqns);
    }

    // number of equations
    unsigned my_n_eqn = my_eqns.size();

    // next we assemble the data into an array of arrays
    // =================================================
    // The idea behind this sparse assembly routine is to use an array of
    // arrays for the entries in each complete matrix. And a second
    // array of arrays stores the global row (or column) indeces.

    // Set up two vector of vectors to store the entries of each  matrix,
    // indexed by either the column or row. The entries of the vector for
    // each matrix correspond to all the rows or columns of that matrix.
    Vector<unsigned**> matrix_col_indices(n_matrix);
    Vector<double**> matrix_values(n_matrix);

    // Loop over the number of matrices being assembled and resize
    // each vector of vectors to the number of rows or columns of the matrix
    for (unsigned m = 0; m < n_matrix; m++)
    {
      matrix_col_indices[m] = new unsigned*[my_n_eqn];
      matrix_values[m] = new double*[my_n_eqn];
      for (unsigned i = 0; i < my_n_eqn; i++)
      {
        matrix_col_indices[m][i] = 0;
        matrix_values[m][i] = 0;
      }
    }

    // Resize the residuals vectors
    Vector<double*> residuals_data(n_vector);
    for (unsigned v = 0; v < n_vector; v++)
    {
      residuals_data[v] = new double[my_n_eqn];
      for (unsigned i = 0; i < my_n_eqn; i++)
      {
        residuals_data[v][i] = 0;
      }
    }

    // Storage for assembly time for elements
    double t_assemble_start = 0.0;

    // Storage for assembly times
    if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
    {
      Elemental_assembly_time.resize(n_elements);
    }

    // number of coefficients in each row
    Vector<Vector<unsigned>> ncoef(n_matrix);
    for (unsigned m = 0; m < n_matrix; m++)
    {
      ncoef[m].resize(my_n_eqn, 0);
    }

    // Sparse_assemble_with_arrays_previous_allocation stores the number of
    // coefs in each row.
    // if a matrix of this size has not been assembled before then resize this
    // storage
    if (Sparse_assemble_with_arrays_previous_allocation.size() == 0)
    {
      Sparse_assemble_with_arrays_previous_allocation.resize(n_matrix);
      for (unsigned m = 0; m < n_matrix; m++)
      {
        Sparse_assemble_with_arrays_previous_allocation[m].resize(my_n_eqn, 0);
      }
    }


    // assemble and populate an array based storage scheme
    {
      // Allocate local storage for the element's contribution to the
      // residuals vectors and system matrices of the size of the maximum
      // number of dofs in any element
      // This means that the storage will only be allocated (and deleted) once
      Vector<Vector<double>> el_residuals(n_vector);
      Vector<DenseMatrix<double>> el_jacobian(n_matrix);

      // Loop over the elements
      for (unsigned long e = el_lo; e < el_hi_plus_one; e++)
      {
        // Time it?
        if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
        {
          t_assemble_start = TimingHelpers::timer();
        }

        // Get the pointer to the element
        GeneralisedElement* elem_pt = mesh_pt()->element_pt(e);

        // Ignore halo elements
        if (!elem_pt->is_halo())
        {
          // Find number of degrees of freedom in the element
          const unsigned nvar = assembly_handler_pt->ndof(elem_pt);

          // Resize the storage for elemental jacobian and residuals
          for (unsigned v = 0; v < n_vector; v++)
          {
            el_residuals[v].resize(nvar);
          }
          for (unsigned m = 0; m < n_matrix; m++)
          {
            el_jacobian[m].resize(nvar);
          }

          // Now get the residuals and jacobian for the element
          assembly_handler_pt->get_all_vectors_and_matrices(
            elem_pt, el_residuals, el_jacobian);

          //---------------Insert the values into the vectors--------------

          // Loop over the first index of local variables
          for (unsigned i = 0; i < nvar; i++)
          {
            // Get the local equation number
            unsigned global_eqn_number =
              assembly_handler_pt->eqn_number(elem_pt, i);

            // determine the element number in my set of eqns using the
            // bisection method
            int left = 0;
            int right = my_n_eqn - 1;
            int eqn_number = right / 2;
            while (my_eqns[eqn_number] != global_eqn_number)
            {
              if (left == right)
              {
                // Check that the residuals associated with the
                // eqn number that can't be found are all zero
                bool broken = false;
                for (unsigned v = 0; v < n_vector; v++)
                {
                  if (el_residuals[v][i] != 0.0)
                  {
                    broken = true;
                    break;
                  }
                }

                // Now loop over the other index to check the entries
                // in the appropriate row of the Jacobians are zero too
                for (unsigned j = 0; j < nvar; j++)
                {
                  // Get the number of the unknown
                  // unsigned unknown =
                  // assembly_handler_pt->eqn_number(elem_pt,j);

                  // Loop over the matrices
                  // If it's compressed row storage, then our vector of maps
                  // is indexed by row (equation number)
                  for (unsigned m = 0; m < n_matrix; m++)
                  {
                    // Get the value of the matrix at this point
                    double value = el_jacobian[m](i, j);
                    if (value != 0.0)
                    {
                      broken = true;
                      break;
                    }
                    if (broken) break;
                  }
                }

                if (broken)
                {
                  std::ostringstream error_stream;
                  error_stream
                    << "Internal Error: " << std::endl
                    << "Could not find global equation number "
                    << global_eqn_number
                    << " in my_eqns vector of equation numbers but\n"
                    << "at least one entry in the residual vector is nonzero.";
                  throw OomphLibError(error_stream.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
                else
                {
                  break;
                }
              }
              if (my_eqns[eqn_number] > global_eqn_number)
              {
                right = std::max(eqn_number - 1, left);
              }
              else
              {
                left = std::min(eqn_number + 1, right);
              }
              eqn_number = (right + left) / 2;
            }

            // Add the contribution to the residuals
            for (unsigned v = 0; v < n_vector; v++)
            {
              // Fill in each residuals vector
              residuals_data[v][eqn_number] += el_residuals[v][i];
            }

            // Now loop over the other index
            for (unsigned j = 0; j < nvar; j++)
            {
              // Get the number of the unknown
              unsigned unknown = assembly_handler_pt->eqn_number(elem_pt, j);

              // Loop over the matrices
              // If it's compressed row storage, then our vector of maps
              // is indexed by row (equation number)
              for (unsigned m = 0; m < n_matrix; m++)
              {
                // Get the value of the matrix at this point
                double value = el_jacobian[m](i, j);
                // Only bother to add to the vector if it's non-zero
                if (std::fabs(value) > Numerical_zero_for_sparse_assembly)
                {
                  // number of entrys in this row
                  const unsigned size = ncoef[m][eqn_number];

                  // if no data has been allocated for this row then allocate
                  if (size == 0)
                  {
                    // do we have previous allocation data
                    if (Sparse_assemble_with_arrays_previous_allocation
                          [m][eqn_number] != 0)
                    {
                      matrix_col_indices[m][eqn_number] = new unsigned
                        [Sparse_assemble_with_arrays_previous_allocation
                           [m][eqn_number]];

                      matrix_values[m][eqn_number] = new double
                        [Sparse_assemble_with_arrays_previous_allocation
                           [m][eqn_number]];
                    }
                    else
                    {
                      matrix_col_indices[m][eqn_number] = new unsigned
                        [Sparse_assemble_with_arrays_initial_allocation];

                      matrix_values[m][eqn_number] = new double
                        [Sparse_assemble_with_arrays_initial_allocation];

                      Sparse_assemble_with_arrays_previous_allocation
                        [m][eqn_number] =
                          Sparse_assemble_with_arrays_initial_allocation;
                    }
                  }

                  // next add the data
                  for (unsigned k = 0; k <= size; k++)
                  {
                    if (k == size)
                    {
                      // do we need to allocate more storage
                      if (Sparse_assemble_with_arrays_previous_allocation
                            [m][eqn_number] == ncoef[m][eqn_number])
                      {
                        unsigned new_allocation =
                          ncoef[m][eqn_number] +
                          Sparse_assemble_with_arrays_allocation_increment;
                        double* new_values = new double[new_allocation];
                        unsigned* new_indices = new unsigned[new_allocation];
                        for (unsigned c = 0; c < ncoef[m][eqn_number]; c++)
                        {
                          new_values[c] = matrix_values[m][eqn_number][c];
                          new_indices[c] = matrix_col_indices[m][eqn_number][c];
                        }
                        delete[] matrix_values[m][eqn_number];
                        delete[] matrix_col_indices[m][eqn_number];

                        matrix_values[m][eqn_number] = new_values;
                        matrix_col_indices[m][eqn_number] = new_indices;

                        Sparse_assemble_with_arrays_previous_allocation
                          [m][eqn_number] = new_allocation;
                      }
                      // and now add the data
                      unsigned entry = ncoef[m][eqn_number];
                      ncoef[m][eqn_number]++;
                      matrix_col_indices[m][eqn_number][entry] = unknown;
                      matrix_values[m][eqn_number][entry] = value;
                      break;
                    }
                    else if (matrix_col_indices[m][eqn_number][k] == unknown)
                    {
                      matrix_values[m][eqn_number][k] += value;
                      break;
                    }
                  }
                } // numerical zero check
              } // End of loop over matrices
            }
          }
        } // endif halo element

        // Time it?
        if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
        {
          Elemental_assembly_time[e] =
            TimingHelpers::timer() - t_assemble_start;
        }
      } // End of loop over the elements
    } // End of vector assembly


    // Doc?
    double t_end = 0.0;
    double t_local = 0.0;
    double t_max = 0.0;
    double t_min = 0.0;
    double t_sum = 0.0;
    if (Doc_imbalance_in_parallel_assembly)
    {
      t_end = TimingHelpers::timer();
      t_local = t_end - t_start;
      t_max = 0.0;
      t_min = 0.0;
      t_sum = 0.0;
      MPI_Allreduce(&t_local,
                    &t_max,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX,
                    this->communicator_pt()->mpi_comm());
      MPI_Allreduce(&t_local,
                    &t_min,
                    1,
                    MPI_DOUBLE,
                    MPI_MIN,
                    this->communicator_pt()->mpi_comm());
      MPI_Allreduce(&t_local,
                    &t_sum,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    this->communicator_pt()->mpi_comm());
      double imbalance = (t_max - t_min) / (t_sum / double(nproc)) * 100.0;

      if (doing_residuals)
      {
        oomph_info << "\nCPU for residual computation (loc/max/min/imbal): ";
      }
      else
      {
        oomph_info << "\nCPU for Jacobian computation (loc/max/min/imbal): ";
      }
      oomph_info << t_local << " " << t_max << " " << t_min << " " << imbalance
                 << "%\n";

      t_start = TimingHelpers::timer();
    }


    // Adjust number of coefficients in each row
    for (unsigned m = 0; m < n_matrix; m++)
    {
      unsigned max = 0;
      unsigned min = INT_MAX;
      unsigned sum = 0;
      unsigned sum_total = 0;
      for (unsigned e = 0; e < my_n_eqn; e++)
      {
        sum += ncoef[m][e];
        sum_total += Sparse_assemble_with_arrays_previous_allocation[m][e];
        if (ncoef[m][e] > max) max = ncoef[m][e];
        if (ncoef[m][e] < min) min = ncoef[m][e];

        // Now shrink the storage to what we actually need
        unsigned new_allocation = ncoef[m][e];
        double* new_values = new double[new_allocation];
        unsigned* new_indices = new unsigned[new_allocation];
        for (unsigned c = 0; c < ncoef[m][e]; c++)
        {
          new_values[c] = matrix_values[m][e][c];
          new_indices[c] = matrix_col_indices[m][e][c];
        }
        delete[] matrix_values[m][e];
        delete[] matrix_col_indices[m][e];

        matrix_values[m][e] = new_values;
        matrix_col_indices[m][e] = new_indices;
      }
    }


    // Postprocess timing information and re-allocate distribution of
    // elements during subsequent assemblies.
    if ((!doing_residuals) && (!Problem_has_been_distributed) &&
        Must_recompute_load_balance_for_assembly)
    {
      recompute_load_balanced_assembly();
    }

    // We have determined load balancing for current setup.
    // This can remain the same until assign_eqn_numbers() is called
    // again -- the flag is re-set to true there.
    if ((!doing_residuals) && Must_recompute_load_balance_for_assembly)
    {
      Must_recompute_load_balance_for_assembly = false;
    }


    // next we compute the number of equations and number of non-zeros to be
    // sent to each processor, and send/recv that information
    // =====================================================================

    // determine the number of eqns to be sent to each processor
    Vector<unsigned> n_eqn_for_proc(nproc, 0);
    Vector<unsigned> first_eqn_element_for_proc(nproc, 0);
    // If no equations are assembled then we don't need to do any of this
    if (my_n_eqn > 0)
    {
      unsigned current_p = target_dist_pt->rank_of_global_row(my_eqns[0]);
      first_eqn_element_for_proc[current_p] = 0;
      n_eqn_for_proc[current_p] = 1;
      for (unsigned i = 1; i < my_n_eqn; i++)
      {
        unsigned next_p = target_dist_pt->rank_of_global_row(my_eqns[i]);
        if (next_p != current_p)
        {
          current_p = next_p;
          first_eqn_element_for_proc[current_p] = i;
        }
        n_eqn_for_proc[current_p]++;
      }
    }

    // determine the number of non-zeros to be sent to each processor for each
    // matrix (if n_eqn_for_proc[p]=0, then nothing will be assembled)
    DenseMatrix<unsigned> nnz_for_proc(nproc, n_matrix, 0);
    for (unsigned p = 0; p < nproc; p++)
    {
      int first_eqn_element = first_eqn_element_for_proc[p];
      int last_eqn_element = (int)(first_eqn_element + n_eqn_for_proc[p]) - 1;
      for (unsigned m = 0; m < n_matrix; m++)
      {
        for (int i = first_eqn_element; i <= last_eqn_element; i++)
        {
          nnz_for_proc(p, m) += ncoef[m][i];
        }
      }
    }

    // next post the sends and recvs to the corresponding processors
    Vector<unsigned*> temp_send_storage(nproc);
    Vector<unsigned*> temp_recv_storage(nproc);
    Vector<MPI_Request> send_nnz_reqs;
    Vector<MPI_Request> recv_nnz_reqs;
    for (unsigned p = 0; p < nproc; p++)
    {
      if (p != my_rank)
      {
        temp_send_storage[p] = new unsigned[n_matrix + 1];
        temp_send_storage[p][0] = n_eqn_for_proc[p];
        for (unsigned m = 0; m < n_matrix; m++)
        {
          temp_send_storage[p][m + 1] = nnz_for_proc(p, m);
        }
        MPI_Request sreq;
        MPI_Isend(temp_send_storage[p],
                  n_matrix + 1,
                  MPI_UNSIGNED,
                  p,
                  0,
                  Communicator_pt->mpi_comm(),
                  &sreq);
        send_nnz_reqs.push_back(sreq);
        temp_recv_storage[p] = new unsigned[n_matrix + 1];
        MPI_Request rreq;
        MPI_Irecv(temp_recv_storage[p],
                  n_matrix + 1,
                  MPI_UNSIGNED,
                  p,
                  0,
                  Communicator_pt->mpi_comm(),
                  &rreq);
        recv_nnz_reqs.push_back(rreq);
      }
    }

    // assemble the data to be sent to each processor
    // ==============================================

    // storage
    Vector<unsigned*> eqns_for_proc(nproc);
    DenseMatrix<double*> residuals_for_proc(nproc, n_vector);
    DenseMatrix<unsigned*> row_start_for_proc(nproc, n_matrix);
    DenseMatrix<unsigned*> column_indices_for_proc(nproc, n_matrix);
    DenseMatrix<double*> values_for_proc(nproc, n_matrix);

    // equation numbers
    for (unsigned p = 0; p < nproc; p++)
    {
      unsigned n_eqns_p = n_eqn_for_proc[p];
      if (n_eqns_p > 0)
      {
        unsigned first_eqn_element = first_eqn_element_for_proc[p];
        unsigned first_row = target_dist_pt->first_row(p);
        eqns_for_proc[p] = new unsigned[n_eqns_p];
        for (unsigned i = 0; i < n_eqns_p; i++)
        {
          eqns_for_proc[p][i] = my_eqns[i + first_eqn_element] - first_row;
        }
      }
    }

    // residuals for p
    for (unsigned v = 0; v < n_vector; v++)
    {
      for (unsigned p = 0; p < nproc; p++)
      {
        unsigned n_eqns_p = n_eqn_for_proc[p];
        if (n_eqns_p > 0)
        {
          unsigned first_eqn_element = first_eqn_element_for_proc[p];
          residuals_for_proc(p, v) = new double[n_eqns_p];
          for (unsigned i = 0; i < n_eqns_p; i++)
          {
            residuals_for_proc(p, v)[i] =
              residuals_data[v][first_eqn_element + i];
          }
        }
      }
      delete[] residuals_data[v];
    }

    // matrices for p
    for (unsigned m = 0; m < n_matrix; m++)
    {
      for (unsigned p = 0; p < nproc; p++)
      {
        unsigned n_eqns_p = n_eqn_for_proc[p];
        if (n_eqns_p > 0)
        {
          unsigned first_eqn_element = first_eqn_element_for_proc[p];
          row_start_for_proc(p, m) = new unsigned[n_eqns_p + 1];
          column_indices_for_proc(p, m) = new unsigned[nnz_for_proc(p, m)];
          values_for_proc(p, m) = new double[nnz_for_proc(p, m)];
          unsigned entry = 0;
          for (unsigned i = 0; i < n_eqns_p; i++)
          {
            row_start_for_proc(p, m)[i] = entry;
            unsigned n_coef_in_row = ncoef[m][first_eqn_element + i];
            for (unsigned j = 0; j < n_coef_in_row; j++)
            {
              column_indices_for_proc(p, m)[entry] =
                matrix_col_indices[m][i + first_eqn_element][j];
              values_for_proc(p, m)[entry] =
                matrix_values[m][i + first_eqn_element][j];
              entry++;
            }
          }
          row_start_for_proc(p, m)[n_eqns_p] = entry;
        }
      }
      for (unsigned i = 0; i < my_n_eqn; i++)
      {
        delete[] matrix_col_indices[m][i];
        delete[] matrix_values[m][i];
      }
      delete[] matrix_col_indices[m];
      delete[] matrix_values[m];
    }

    // need to wait for the recv nnzs to complete
    // before we can allocate storage for the matrix recvs
    // ===================================================

    // recv and copy the datafrom the recv storage to
    // + nnz_from_proc
    // + n_eqn_from_proc
    Vector<MPI_Status> recv_nnz_stat(nproc - 1);
    MPI_Waitall(nproc - 1, &recv_nnz_reqs[0], &recv_nnz_stat[0]);
    Vector<unsigned> n_eqn_from_proc(nproc);
    DenseMatrix<unsigned> nnz_from_proc(nproc, n_matrix);
    for (unsigned p = 0; p < nproc; p++)
    {
      if (p != my_rank)
      {
        n_eqn_from_proc[p] = temp_recv_storage[p][0];
        for (unsigned m = 0; m < n_matrix; m++)
        {
          nnz_from_proc(p, m) = temp_recv_storage[p][m + 1];
        }
        delete[] temp_recv_storage[p];
      }
      else
      {
        n_eqn_from_proc[p] = n_eqn_for_proc[p];
        for (unsigned m = 0; m < n_matrix; m++)
        {
          nnz_from_proc(p, m) = nnz_for_proc(p, m);
        }
      }
    }
    recv_nnz_stat.clear();
    recv_nnz_reqs.clear();

    // allocate the storage for the data to be recv and post the sends recvs
    // =====================================================================

    // storage
    Vector<unsigned*> eqns_from_proc(nproc);
    DenseMatrix<double*> residuals_from_proc(nproc, n_vector);
    DenseMatrix<unsigned*> row_start_from_proc(nproc, n_matrix);
    DenseMatrix<unsigned*> column_indices_from_proc(nproc, n_matrix);
    DenseMatrix<double*> values_from_proc(nproc, n_matrix);

    // allocate and post sends and recvs
    double base;
    MPI_Aint communication_base;
    MPI_Get_address(&base, &communication_base);
    unsigned n_comm_types = 1 + 1 * n_vector + 3 * n_matrix;
    Vector<MPI_Request> recv_reqs;
    Vector<MPI_Request> send_reqs;
    for (unsigned p = 0; p < nproc; p++)
    {
      if (p != my_rank)
      {
        // allocate
        if (n_eqn_from_proc[p] > 0)
        {
          eqns_from_proc[p] = new unsigned[n_eqn_from_proc[p]];
          for (unsigned v = 0; v < n_vector; v++)
          {
            residuals_from_proc(p, v) = new double[n_eqn_from_proc[p]];
          }
          for (unsigned m = 0; m < n_matrix; m++)
          {
            row_start_from_proc(p, m) = new unsigned[n_eqn_from_proc[p] + 1];
            column_indices_from_proc(p, m) = new unsigned[nnz_from_proc(p, m)];
            values_from_proc(p, m) = new double[nnz_from_proc(p, m)];
          }
        }

        // recv
        if (n_eqn_from_proc[p] > 0)
        {
          MPI_Datatype types[n_comm_types];
          MPI_Aint offsets[n_comm_types];
          int count[n_comm_types];
          int pt = 0;

          // equations
          count[pt] = 1;
          MPI_Get_address(eqns_from_proc[p], &offsets[pt]);
          offsets[pt] -= communication_base;
          MPI_Type_contiguous(n_eqn_from_proc[p], MPI_UNSIGNED, &types[pt]);
          MPI_Type_commit(&types[pt]);
          pt++;

          // vectors
          for (unsigned v = 0; v < n_vector; v++)
          {
            count[pt] = 1;
            MPI_Get_address(residuals_from_proc(p, v), &offsets[pt]);
            offsets[pt] -= communication_base;
            MPI_Type_contiguous(n_eqn_from_proc[p], MPI_DOUBLE, &types[pt]);
            MPI_Type_commit(&types[pt]);
            pt++;
          }

          // matrices
          for (unsigned m = 0; m < n_matrix; m++)
          {
            // row start
            count[pt] = 1;
            MPI_Get_address(row_start_from_proc(p, m), &offsets[pt]);
            offsets[pt] -= communication_base;
            MPI_Type_contiguous(
              n_eqn_from_proc[p] + 1, MPI_UNSIGNED, &types[pt]);
            MPI_Type_commit(&types[pt]);
            pt++;


            // column indices
            count[pt] = 1;
            MPI_Get_address(column_indices_from_proc(p, m), &offsets[pt]);
            offsets[pt] -= communication_base;
            MPI_Type_contiguous(nnz_from_proc(p, m), MPI_UNSIGNED, &types[pt]);
            MPI_Type_commit(&types[pt]);
            pt++;

            // values
            count[pt] = 1;
            MPI_Get_address(values_from_proc(p, m), &offsets[pt]);
            offsets[pt] -= communication_base;
            MPI_Type_contiguous(nnz_from_proc(p, m), MPI_DOUBLE, &types[pt]);
            MPI_Type_commit(&types[pt]);
            pt++;
          }

          // build the combined type
          MPI_Datatype recv_type;
          MPI_Type_create_struct(
            n_comm_types, count, offsets, types, &recv_type);
          MPI_Type_commit(&recv_type);
          for (unsigned t = 0; t < n_comm_types; t++)
          {
            MPI_Type_free(&types[t]);
          }
          MPI_Request req;
          MPI_Irecv(
            &base, 1, recv_type, p, 1, Communicator_pt->mpi_comm(), &req);
          MPI_Type_free(&recv_type);
          recv_reqs.push_back(req);
        }

        // send
        if (n_eqn_for_proc[p] > 0)
        {
          MPI_Datatype types[n_comm_types];
          MPI_Aint offsets[n_comm_types];
          int count[n_comm_types];
          int pt = 0;

          // equations
          count[pt] = 1;
          MPI_Get_address(eqns_for_proc[p], &offsets[pt]);
          offsets[pt] -= communication_base;
          MPI_Type_contiguous(n_eqn_for_proc[p], MPI_UNSIGNED, &types[pt]);
          MPI_Type_commit(&types[pt]);
          pt++;

          // vectors
          for (unsigned v = 0; v < n_vector; v++)
          {
            count[pt] = 1;
            MPI_Get_address(residuals_for_proc(p, v), &offsets[pt]);
            offsets[pt] -= communication_base;
            MPI_Type_contiguous(n_eqn_for_proc[p], MPI_DOUBLE, &types[pt]);
            MPI_Type_commit(&types[pt]);
            pt++;
          }

          // matrices
          for (unsigned m = 0; m < n_matrix; m++)
          {
            // row start
            count[pt] = 1;
            MPI_Get_address(row_start_for_proc(p, m), &offsets[pt]);
            offsets[pt] -= communication_base;
            MPI_Type_contiguous(
              n_eqn_for_proc[p] + 1, MPI_UNSIGNED, &types[pt]);
            MPI_Type_commit(&types[pt]);
            pt++;


            // column indices
            count[pt] = 1;
            MPI_Get_address(column_indices_for_proc(p, m), &offsets[pt]);
            offsets[pt] -= communication_base;
            MPI_Type_contiguous(nnz_for_proc(p, m), MPI_UNSIGNED, &types[pt]);
            MPI_Type_commit(&types[pt]);
            pt++;

            // values
            count[pt] = 1;
            MPI_Get_address(values_for_proc(p, m), &offsets[pt]);
            offsets[pt] -= communication_base;
            MPI_Type_contiguous(nnz_for_proc(p, m), MPI_DOUBLE, &types[pt]);
            MPI_Type_commit(&types[pt]);
            pt++;
          }

          // build the combined type
          MPI_Datatype send_type;
          MPI_Type_create_struct(
            n_comm_types, count, offsets, types, &send_type);
          MPI_Type_commit(&send_type);
          for (unsigned t = 0; t < n_comm_types; t++)
          {
            MPI_Type_free(&types[t]);
          }
          MPI_Request req;
          MPI_Isend(
            &base, 1, send_type, p, 1, Communicator_pt->mpi_comm(), &req);
          MPI_Type_free(&send_type);
          send_reqs.push_back(req);
        }
      }
      // otherwise send to self
      else
      {
        eqns_from_proc[p] = eqns_for_proc[p];
        for (unsigned v = 0; v < n_vector; v++)
        {
          residuals_from_proc(p, v) = residuals_for_proc(p, v);
        }
        for (unsigned m = 0; m < n_matrix; m++)
        {
          row_start_from_proc(p, m) = row_start_for_proc(p, m);
          column_indices_from_proc(p, m) = column_indices_for_proc(p, m);
          values_from_proc(p, m) = values_for_proc(p, m);
        }
      }
    }

    // wait for the recvs to complete
    unsigned n_recv_req = recv_reqs.size();
    if (n_recv_req > 0)
    {
      Vector<MPI_Status> recv_stat(n_recv_req);
      MPI_Waitall(n_recv_req, &recv_reqs[0], &recv_stat[0]);
    }

    // ==============================================
    unsigned target_nrow_local = target_dist_pt->nrow_local();

    // loop over the matrices
    for (unsigned m = 0; m < n_matrix; m++)
    {
      // allocate row_start
      row_start[m] = new int[target_nrow_local + 1];
      row_start[m][0] = 0;

      // initially allocate storage based on the maximum number of non-zeros
      // from any one processor
      unsigned nnz_allocation = Parallel_sparse_assemble_previous_allocation;
      for (unsigned p = 0; p < nproc; p++)
      {
        nnz_allocation = std::max(nnz_allocation, nnz_from_proc(p, m));
      }
      Vector<double*> values_chunk(1);
      values_chunk[0] = new double[nnz_allocation];
      Vector<int*> column_indices_chunk(1);
      column_indices_chunk[0] = new int[nnz_allocation];
      Vector<unsigned> ncoef_in_chunk(1, 0);
      Vector<unsigned> size_of_chunk(1, 0);
      size_of_chunk[0] = nnz_allocation;
      unsigned current_chunk = 0;

      // for each row on this processor
      for (unsigned i = 0; i < target_nrow_local; i++)
      {
        row_start[m][i] = 0;

        // determine the processors that this row is on
        Vector<int> row_on_proc(nproc, -1);
        for (unsigned p = 0; p < nproc; p++)
        {
          if (n_eqn_from_proc[p] == 0)
          {
            row_on_proc[p] = -1;
          }
          else
          {
            int left = 0;
            int right = n_eqn_from_proc[p] - 1;
            int midpoint = right / 2;
            bool complete = false;
            while (!complete)
            {
              midpoint = (right + left) / 2;
              if (midpoint > right)
              {
                midpoint = right;
              }
              if (midpoint < left)
              {
                midpoint = left;
              }
              if (left == right)
              {
                if (eqns_from_proc[p][midpoint] == i)
                {
                  midpoint = left;
                }
                else
                {
                  midpoint = -1;
                }
                complete = true;
              }
              else if (eqns_from_proc[p][midpoint] == i)
              {
                complete = true;
              }
              else if (eqns_from_proc[p][midpoint] > i)
              {
                right = std::max(midpoint - 1, left);
              }
              else
              {
                left = std::min(midpoint + 1, right);
              }
            }
            row_on_proc[p] = midpoint;
          }
        }

        // for each processor build this row of the matrix
        unsigned check_first = ncoef_in_chunk[current_chunk];
        unsigned check_last = check_first;
        for (unsigned p = 0; p < nproc; p++)
        {
          if (row_on_proc[p] != -1)
          {
            int row = row_on_proc[p];
            unsigned first = row_start_from_proc(p, m)[row];
            unsigned last = row_start_from_proc(p, m)[row + 1];
            for (unsigned l = first; l < last; l++)
            {
              bool done = false;
              for (unsigned j = check_first; j <= check_last && !done; j++)
              {
                if (j == check_last)
                {
                  // is this temp array full, do we need to allocate
                  // a new temp array
                  if (ncoef_in_chunk[current_chunk] ==
                      size_of_chunk[current_chunk])
                  {
                    // number of chunks allocated
                    unsigned n_chunk = values_chunk.size();

                    // determine the number of non-zeros added so far
                    // (excluding the current row)
                    unsigned nnz_so_far = 0;
                    for (unsigned c = 0; c < n_chunk; c++)
                    {
                      nnz_so_far += ncoef_in_chunk[c];
                    }
                    nnz_so_far -= row_start[m][i];

                    // average number of non-zeros per row
                    unsigned avg_nnz = nnz_so_far / (i + 1);

                    // number of rows left +1
                    unsigned nrows_left = target_nrow_local - i;

                    // allocation for next chunk
                    unsigned next_chunk_size =
                      avg_nnz * nrows_left + row_start[m][i];

                    // allocate storage in next chunk
                    current_chunk++;
                    n_chunk++;
                    values_chunk.resize(n_chunk);
                    values_chunk[current_chunk] = new double[next_chunk_size];
                    column_indices_chunk.resize(n_chunk);
                    column_indices_chunk[current_chunk] =
                      new int[next_chunk_size];
                    size_of_chunk.resize(n_chunk);
                    size_of_chunk[current_chunk] = next_chunk_size;
                    ncoef_in_chunk.resize(n_chunk);

                    // copy current row from previous chunk to new chunk
                    for (unsigned k = check_first; k < check_last; k++)
                    {
                      values_chunk[current_chunk][k - check_first] =
                        values_chunk[current_chunk - 1][k];
                      column_indices_chunk[current_chunk][k - check_first] =
                        column_indices_chunk[current_chunk - 1][k];
                    }
                    ncoef_in_chunk[current_chunk - 1] -= row_start[m][i];
                    ncoef_in_chunk[current_chunk] = row_start[m][i];

                    // update first_check and last_check
                    check_first = 0;
                    check_last = row_start[m][i];
                    j = check_last;
                  }

                  // add the coefficient
                  values_chunk[current_chunk][j] = values_from_proc(p, m)[l];
                  column_indices_chunk[current_chunk][j] =
                    column_indices_from_proc(p, m)[l];
                  ncoef_in_chunk[current_chunk]++;
                  row_start[m][i]++;
                  check_last++;
                  done = true;
                }
                else if (column_indices_chunk[current_chunk][j] ==
                         (int)column_indices_from_proc(p, m)[l])
                {
                  values_chunk[current_chunk][j] += values_from_proc(p, m)[l];
                  done = true;
                }
              }
            }
          }
        }
      }

      // delete recv data for this matrix
      for (unsigned p = 0; p < nproc; p++)
      {
        if (n_eqn_from_proc[p] > 0)
        {
          delete[] row_start_from_proc(p, m);
          delete[] column_indices_from_proc(p, m);
          delete[] values_from_proc(p, m);
        }
      }

      // next we take the chunk base storage of the column indices and values
      // and copy into a single contiguous block of memory
      // ====================================================================
      unsigned n_chunk = values_chunk.size();
      nnz[m] = 0;
      for (unsigned c = 0; c < n_chunk; c++)
      {
        nnz[m] += ncoef_in_chunk[c];
      }
      Parallel_sparse_assemble_previous_allocation = nnz[m];

      // allocate
      values[m] = new double[nnz[m]];
      column_indices[m] = new int[nnz[m]];

      // copy
      unsigned pt = 0;
      for (unsigned c = 0; c < n_chunk; c++)
      {
        unsigned nc = ncoef_in_chunk[c];
        for (unsigned i = 0; i < nc; i++)
        {
          values[m][pt + i] = values_chunk[c][i];
          column_indices[m][pt + i] = column_indices_chunk[c][i];
        }
        pt += nc;
        delete[] values_chunk[c];
        delete[] column_indices_chunk[c];
      }

      // the row_start vector currently contains the number of coefs in each
      // row. Update
      // ===================================================================
      unsigned g = row_start[m][0];
      row_start[m][0] = 0;
      for (unsigned i = 1; i < target_nrow_local; i++)
      {
        unsigned h = g + row_start[m][i];
        row_start[m][i] = g;
        g = h;
      }
      row_start[m][target_nrow_local] = g;
    }

    // next accumulate the residuals
    for (unsigned v = 0; v < n_vector; v++)
    {
      residuals[v] = new double[target_nrow_local];
      for (unsigned i = 0; i < target_nrow_local; i++)
      {
        residuals[v][i] = 0;
      }
      for (unsigned p = 0; p < nproc; p++)
      {
        if (n_eqn_from_proc[p] > 0)
        {
          unsigned n_eqn_p = n_eqn_from_proc[p];
          for (unsigned i = 0; i < n_eqn_p; i++)
          {
            residuals[v][eqns_from_proc[p][i]] += residuals_from_proc(p, v)[i];
          }
          delete[] residuals_from_proc(p, v);
        }
      }
    }

    // delete list of eqns from proc
    for (unsigned p = 0; p < nproc; p++)
    {
      if (n_eqn_from_proc[p] > 0)
      {
        delete[] eqns_from_proc[p];
      }
    }

    // and wait for sends to complete
    Vector<MPI_Status> send_nnz_stat(nproc - 1);
    MPI_Waitall(nproc - 1, &send_nnz_reqs[0], &send_nnz_stat[0]);
    for (unsigned p = 0; p < nproc; p++)
    {
      if (p != my_rank)
      {
        delete[] temp_send_storage[p];
      }
    }
    send_nnz_stat.clear();
    send_nnz_reqs.clear();

    // wait for the matrix data sends to complete and delete the data
    unsigned n_send_reqs = send_reqs.size();
    if (n_send_reqs > 0)
    {
      Vector<MPI_Status> send_stat(n_send_reqs);
      MPI_Waitall(n_send_reqs, &send_reqs[0], &send_stat[0]);
      for (unsigned p = 0; p < nproc; p++)
      {
        if (p != my_rank)
        {
          if (n_eqn_for_proc[p])
          {
            delete[] eqns_for_proc[p];
            for (unsigned m = 0; m < n_matrix; m++)
            {
              delete[] row_start_for_proc(p, m);
              delete[] column_indices_for_proc(p, m);
              delete[] values_for_proc(p, m);
            }
            for (unsigned v = 0; v < n_vector; v++)
            {
              delete[] residuals_for_proc(p, v);
            }
          }
        }
      }
    }

    // Doc?
    if (Doc_imbalance_in_parallel_assembly)
    {
      t_end = TimingHelpers::timer();
      t_local = t_end - t_start;
      t_max = 0.0;
      t_min = 0.0;
      t_sum = 0.0;
      MPI_Allreduce(&t_local,
                    &t_max,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX,
                    this->communicator_pt()->mpi_comm());
      MPI_Allreduce(&t_local,
                    &t_min,
                    1,
                    MPI_DOUBLE,
                    MPI_MIN,
                    this->communicator_pt()->mpi_comm());
      MPI_Allreduce(&t_local,
                    &t_sum,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    this->communicator_pt()->mpi_comm());
      double imbalance = (t_max - t_min) / (t_sum / double(nproc)) * 100.0;
      if (doing_residuals)
      {
        oomph_info << "CPU for residual distribut.  (loc/max/min/imbal): ";
      }
      else
      {
        oomph_info << "CPU for Jacobian distribut.  (loc/max/min/imbal): ";
      }
      oomph_info << t_local << " " << t_max << " " << t_min << " " << imbalance
                 << "%\n\n";
    }
  }

#endif


  //================================================================
  /// \short Get the full Jacobian by finite differencing
  //================================================================
  void Problem::get_fd_jacobian(DoubleVector& residuals,
                                DenseMatrix<double>& jacobian)
  {
#ifdef OOMPH_HAS_MPI

    if (Problem_has_been_distributed)
    {
      OomphLibWarning("This is unlikely to work with a distributed problem",
                      " Problem::get_fd_jacobian()",
                      OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // Find number of dofs
    const unsigned long n_dof = ndof();

    // Advanced residuals
    DoubleVector residuals_pls;

    // Get reference residuals
    get_residuals(residuals);

    const double FD_step = 1.0e-8;

    // Make sure the Jacobian is the right size (since we don't care about
    // speed).
    jacobian.resize(n_dof, n_dof);

    // Loop over all dofs
    for (unsigned long jdof = 0; jdof < n_dof; jdof++)
    {
      double backup = *Dof_pt[jdof];
      *Dof_pt[jdof] += FD_step;

      // We're checking if the new values for Dof_pt[] actually
      // solve the entire problem --> update as if problem had
      // been solved
      actions_before_newton_solve();
      actions_before_newton_convergence_check();
      actions_after_newton_solve();

      // Get advanced residuals
      get_residuals(residuals_pls);

      for (unsigned long ieqn = 0; ieqn < n_dof; ieqn++)
      {
        jacobian(ieqn, jdof) =
          (residuals_pls[ieqn] - residuals[ieqn]) / FD_step;
      }

      *Dof_pt[jdof] = backup;
    }

    // Reset problem to state it was in
    actions_before_newton_solve();
    actions_before_newton_convergence_check();
    actions_after_newton_solve();
  }

  //======================================================================
  /// \short Get derivative of the residuals vector wrt a global parameter
  /// This is required in continuation problems
  //=======================================================================
  void Problem::get_derivative_wrt_global_parameter(double* const& parameter_pt,
                                                    DoubleVector& result)
  {
    // If we are doing the calculation analytically then call the appropriate
    // handler and then calling get_residuals
    if (is_dparameter_calculated_analytically(parameter_pt))
    {
      // Locally cache pointer to assembly handler
      AssemblyHandler* const old_assembly_handler_pt = Assembly_handler_pt;
      // Create a new assembly handler that replaces get_residuals by
      // get_dresiduals_dparameter for each element
      Assembly_handler_pt =
        new ParameterDerivativeHandler(old_assembly_handler_pt, parameter_pt);
      // Get the residuals, which will be dresiduals by dparameter
      this->get_residuals(result);
      // Delete the parameter derivative handler
      delete Assembly_handler_pt;
      // Reset the assembly handler to the original handler
      Assembly_handler_pt = old_assembly_handler_pt;

      /*AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;
      //Loop over all the elements
      unsigned long Element_pt_range = Mesh_pt->nelement();
      for(unsigned long e=0;e<Element_pt_range;e++)
       {
        //Get the pointer to the element
        GeneralisedElement* elem_pt = Mesh_pt->element_pt(e);
        //Find number of dofs in the element
        unsigned n_element_dofs = assembly_handler_pt->ndof(elem_pt);
        //Set up an array
        Vector<double> element_residuals(n_element_dofs);
        //Fill the array
        assembly_handler_pt->get_dresiduals_dparameter(elem_pt,parameter_pt,
                                                       element_residuals);
        //Now loop over the dofs and assign values to global Vector
        for(unsigned l=0;l<n_element_dofs;l++)
         {
          result[assembly_handler_pt->eqn_number(elem_pt,l)]
           += element_residuals[l];
          }
          }*/

      // for(unsigned n=0;n<n_dof;n++)
      // {std::cout << "BLA " << n << " " <<  result[n] << "\n";}
    }
    // Otherwise use the finite difference default
    else
    {
      // Get the (global) residuals and store in the result vector
      get_residuals(result);

      // Storage for the new residuals
      DoubleVector newres;

      // Increase the global parameter
      const double FD_step = 1.0e-8;

      // Store the current value of the parameter
      double param_value = *parameter_pt;

      // Increase the parameter
      *parameter_pt += FD_step;

      // Do any possible updates
      actions_after_change_in_global_parameter(parameter_pt);

      // Get the new residuals
      get_residuals(newres);

      // Find the number of local rows
      //(I think it's a global vector, so that should be fine)
      const unsigned ndof_local = result.nrow_local();

      // Do the finite differencing in the local variables
      for (unsigned n = 0; n < ndof_local; ++n)
      {
        result[n] = (newres[n] - result[n]) / FD_step;
      }

      // Reset the value of the parameter
      *parameter_pt = param_value;

      // Do any possible updates
      actions_after_change_in_global_parameter(parameter_pt);
    }
  }


  //======================================================================
  /// Return the product of the global hessian (derivative of Jacobian
  /// matrix  with respect to all variables) with
  /// an eigenvector, Y, and any number of other specified vectors C
  /// (d(J_{ij})/d u_{k}) Y_{j} C_{k}.
  /// This function is used in assembling and solving the augmented systems
  /// associated with bifurcation tracking.
  /// The default implementation is to use finite differences at the global
  /// level.
  //========================================================================
  void Problem::get_hessian_vector_products(
    DoubleVectorWithHaloEntries const& Y,
    Vector<DoubleVectorWithHaloEntries> const& C,
    Vector<DoubleVectorWithHaloEntries>& product)
  {
    // How many vector products must we construct
    const unsigned n_vec = C.size();

    // currently only global (non-distributed) distributions are allowed
    // LinearAlgebraDistribution* dist_pt = new
    // LinearAlgebraDistribution(Communicator_pt,n_dof,false);

    // Cache the assembly hander
    AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

    // Rebuild the results vectors and initialise to zero
    // use the same distribution of the vector Y
    for (unsigned i = 0; i < n_vec; i++)
    {
      product[i].build(Y.distribution_pt(), 0.0);
      product[i].initialise(0.0);
    }

// Setup the halo schemes for the result
#ifdef OOMPH_HAS_MPI
    if (Problem_has_been_distributed)
    {
      for (unsigned i = 0; i < n_vec; i++)
      {
        product[i].build_halo_scheme(this->Halo_scheme_pt);
      }
    }
#endif

    // If we are doing the calculation analytically then call the appropriate
    // handler
    // A better way to do this is probably to hook into the get_residuals
    // framework but with a different member function of the assembly
    // handler
    if (this->are_hessian_products_calculated_analytically())
    {
      // Loop over all the elements
      unsigned long Element_pt_range = Mesh_pt->nelement();
      for (unsigned long e = 0; e < Element_pt_range; e++)
      {
        // Get the pointer to the element
        GeneralisedElement* elem_pt = Mesh_pt->element_pt(e);
// Do not loop over halo elements
#ifdef OOMPH_HAS_MPI
        if (!elem_pt->is_halo())
        {
#endif
          // Find number of dofs in the element
          unsigned n_var = assembly_handler_pt->ndof(elem_pt);
          // Set up a matrix for the input and output
          Vector<double> Y_local(n_var);
          DenseMatrix<double> C_local(n_vec, n_var);
          DenseMatrix<double> product_local(n_vec, n_var);

          // Translate the global input vectors into the local storage
          // Probably horribly inefficient, but otherwise things get really
          // messy at the elemental level
          for (unsigned l = 0; l < n_var; l++)
          {
            // Cache the global equation number
            const unsigned long eqn_number =
              assembly_handler_pt->eqn_number(elem_pt, l);

            Y_local[l] = Y.global_value(eqn_number);
            for (unsigned i = 0; i < n_vec; i++)
            {
              C_local(i, l) = C[i].global_value(eqn_number);
            }
          }

          // Fill the array
          assembly_handler_pt->get_hessian_vector_products(
            elem_pt, Y_local, C_local, product_local);

          // Assign the local results to the global vector
          for (unsigned l = 0; l < n_var; l++)
          {
            const unsigned long eqn_number =
              assembly_handler_pt->eqn_number(elem_pt, l);

            for (unsigned i = 0; i < n_vec; i++)
            {
              product[i].global_value(eqn_number) += product_local(i, l);
              // std::cout << "BLA " << e << " " << i << " "
              //          << l << " " << product_local(i,l) << "\n";
            }
          }
#ifdef OOMPH_HAS_MPI
        }
#endif
      }
    }
    // Otherwise calculate using finite differences by
    // perturbing the jacobian along a particular direction
    else
    {
      // Cache the finite difference step
      /// Alice: My bifurcation tracking converges better with this FD_step
      /// as 1.0e-5. The default value remains at 1.0e-8.
      const double FD_step = FD_step_used_in_get_hessian_vector_products;

      // We can now construct our multipliers
      const unsigned n_dof_local = this->Dof_distribution_pt->nrow_local();
      // Prepare to scale
      double dof_length = 0.0;
      Vector<double> C_length(n_vec, 0.0);

      for (unsigned n = 0; n < n_dof_local; n++)
      {
        if (std::fabs(this->dof(n)) > dof_length)
        {
          dof_length = std::fabs(this->dof(n));
        }
      }

      // C is assumed to have the same distribution as the dofs
      for (unsigned i = 0; i < n_vec; i++)
      {
        for (unsigned n = 0; n < n_dof_local; n++)
        {
          if (std::fabs(C[i][n]) > C_length[i])
          {
            C_length[i] = std::fabs(C[i][n]);
          }
        }
      }

      // Now broadcast the information, if distributed
#ifdef OOMPH_HAS_MPI
      if (Problem_has_been_distributed)
      {
        const unsigned n_length = n_vec + 1;
        double all_length[n_length];
        all_length[0] = dof_length;
        for (unsigned i = 0; i < n_vec; i++)
        {
          all_length[i + 1] = C_length[i];
        }

        // Do the MPI call
        double all_length_reduce[n_length];
        MPI_Allreduce(all_length,
                      all_length_reduce,
                      n_length,
                      MPI_DOUBLE,
                      MPI_MAX,
                      this->communicator_pt()->mpi_comm());

        // Read out the information
        dof_length = all_length_reduce[0];
        for (unsigned i = 0; i < n_vec; i++)
        {
          C_length[i] = all_length_reduce[i + 1];
        }
      }
#endif

      // Form the multipliers
      Vector<double> C_mult(n_vec, 0.0);
      for (unsigned i = 0; i < n_vec; i++)
      {
        C_mult[i] = dof_length / C_length[i];
        C_mult[i] += FD_step;
        C_mult[i] *= FD_step;
      }


      // Dummy vector to stand in the place of the residuals
      Vector<double> dummy_res;

      // Calculate the product of the jacobian matrices, etc by looping over the
      // elements
      const unsigned long n_element = this->mesh_pt()->nelement();
      for (unsigned long e = 0; e < n_element; e++)
      {
        GeneralisedElement* elem_pt = this->mesh_pt()->element_pt(e);
        // Ignore halo's of course
#ifdef OOMPH_HAS_MPI
        if (!elem_pt->is_halo())
        {
#endif
          // Loop over the ndofs in each element
          unsigned n_var = assembly_handler_pt->ndof(elem_pt);
          // Resize the dummy residuals vector
          dummy_res.resize(n_var);
          // Allocate storage for the unperturbed jacobian matrix
          DenseMatrix<double> jac(n_var);
          // Get unperturbed jacobian
          assembly_handler_pt->get_jacobian(elem_pt, dummy_res, jac);

          // Backup the dofs
          Vector<double> dof_bac(n_var);
          for (unsigned n = 0; n < n_var; n++)
          {
            unsigned eqn_number = assembly_handler_pt->eqn_number(elem_pt, n);
            dof_bac[n] = *this->global_dof_pt(eqn_number);
          }

          // Now loop over all vectors C
          for (unsigned i = 0; i < n_vec; i++)
          {
            // Perturb the dofs by the appropriate vector
            for (unsigned n = 0; n < n_var; n++)
            {
              unsigned eqn_number = assembly_handler_pt->eqn_number(elem_pt, n);
              // Perturb by vector C[i]
              *this->global_dof_pt(eqn_number) +=
                C_mult[i] * C[i].global_value(eqn_number);
            }
            actions_before_newton_convergence_check();

            // Allocate storage for the perturbed jacobian
            DenseMatrix<double> jac_C(n_var);

            // Now get the new jacobian
            assembly_handler_pt->get_jacobian(elem_pt, dummy_res, jac_C);

            // Reset the dofs
            for (unsigned n = 0; n < n_var; n++)
            {
              unsigned eqn_number = assembly_handler_pt->eqn_number(elem_pt, n);
              *this->global_dof_pt(eqn_number) = dof_bac[n];
            }
            actions_before_newton_convergence_check();

            // Now work out the products
            for (unsigned n = 0; n < n_var; n++)
            {
              unsigned eqn_number = assembly_handler_pt->eqn_number(elem_pt, n);
              double prod_c = 0.0;
              for (unsigned m = 0; m < n_var; m++)
              {
                unsigned unknown = assembly_handler_pt->eqn_number(elem_pt, m);
                prod_c += (jac_C(n, m) - jac(n, m)) * Y.global_value(unknown);
              }
              // std::cout << "FD   " << e << " " << i << " "
              //          << n << " " << prod_c/C_mult[i] << "\n";
              product[i].global_value(eqn_number) += prod_c / C_mult[i];
            }
          }
#ifdef OOMPH_HAS_MPI
        }
#endif
      } // End of loop over elements
    }

    // If we have a distributed problem then gather all
    // values
#ifdef OOMPH_HAS_MPI
    if (Problem_has_been_distributed)
    {
      // Sum all values if distributed
      for (unsigned i = 0; i < n_vec; i++)
      {
        product[i].sum_all_halo_and_haloed_values();
      }
    }
#endif
  }


  //================================================================
  /// \short Get derivative of an element in the problem wrt a global
  /// parameter, to be used in continuation problems
  //================================================================
  /*void Problem::get_derivative_wrt_global_parameter(
   double* const &parameter_pt,
   GeneralisedElement* const &elem_pt,
   Vector<double> &result)
  {

  #ifdef OOMPH_HAS_MPI

   if (Problem_has_been_distributed)
    {
     OomphLibWarning("This is unlikely to work with a distributed problem",
                     "Problem::get_derivative_wrt_global_parameter()",
                     OOMPH_EXCEPTION_LOCATION);
    }
  #endif

   //Locally cache pointer to assembly handler
   AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

   //Should definitely give this a more global scope
   double FD_Jstep = 1.0e-8;

   //Find the number of variables in the element, e
   unsigned nvar = assembly_handler_pt->ndof(elem_pt);
   //Create storage for residuals
   Vector<double> residuals(nvar), newres(nvar);

   //Get the "original" residuals
   assembly_handler_pt->get_residuals(elem_pt,residuals);

   //Save the old value of the global parameter
   double old_var = *parameter_pt;

   //Increment the value
   *parameter_pt += FD_Jstep;

   //Now do any possible updates
   actions_after_change_in_global_parameter();

   //Get the "new" residuals
   assembly_handler_pt->get_residuals(elem_pt,newres);

   //Do the finite differences
   for(unsigned m=0;m<nvar;m++)
    {
     result[m] = (newres[m] - residuals[m])/FD_Jstep;
    }

   //Reset value of the global parameter
   *parameter_pt = old_var;

   //Now do any possible updates
   actions_after_change_in_global_parameter();
  }*/

  //==================================================================
  /// Solve the eigenproblem
  //==================================================================
  void Problem::solve_eigenproblem(const unsigned& n_eval,
                                   Vector<std::complex<double>>& eigenvalue,
                                   Vector<DoubleVector>& eigenvector,
                                   const bool& steady)
  {
    // If the boolean flag is steady, then make all the timesteppers steady
    // before solving the eigenproblem. This will "switch off" the
    // time-derivative terms in the jacobian matrix
    if (steady)
    {
      // Find out how many timesteppers there are
      const unsigned n_time_steppers = ntime_stepper();

      // Vector of bools to store the is_steady status of the various
      // timesteppers when we came in here
      std::vector<bool> was_steady(n_time_steppers);

      // Loop over them all and make them (temporarily) static
      for (unsigned i = 0; i < n_time_steppers; i++)
      {
        was_steady[i] = time_stepper_pt(i)->is_steady();
        time_stepper_pt(i)->make_steady();
      }

      // Call the Eigenproblem for the eigensolver
      Eigen_solver_pt->solve_eigenproblem(
        this, n_eval, eigenvalue, eigenvector);

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
    }
    // Otherwise if we don't want to make the problem steady, just
    // assemble and solve the eigensystem
    else
    {
      // Call the Eigenproblem for the eigensolver
      Eigen_solver_pt->solve_eigenproblem(
        this, n_eval, eigenvalue, eigenvector);
    }
  }


  //===================================================================
  /// Get the matrices required to solve an eigenproblem
  /// WARNING: temporarily this method only works with non-distributed
  /// matrices
  //===================================================================
  void Problem::get_eigenproblem_matrices(CRDoubleMatrix& mass_matrix,
                                          CRDoubleMatrix& main_matrix,
                                          const double& shift)
  {
    // Three different cases again here:
    // 1) Compiled with MPI, but run in serial
    // 2) Compiled with MPI, but MPI not initialised in driver
    // 3) Serial version


#ifdef PARANOID
    if (mass_matrix.distribution_built() && main_matrix.distribution_built())
    {
      // Check that the distribution of the mass matrix and jacobian match
      if (!(*mass_matrix.distribution_pt() == *main_matrix.distribution_pt()))
      {
        std::ostringstream error_stream;
        error_stream
          << "The distributions of the jacobian and mass matrix are\n"
          << "not the same and they must be.\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      if (mass_matrix.nrow() != this->ndof())
      {
        std::ostringstream error_stream;
        error_stream
          << "mass_matrix has a distribution, but the number of rows is not "
          << "equal to the number of degrees of freedom in the problem.";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      if (main_matrix.nrow() != this->ndof())
      {
        std::ostringstream error_stream;
        error_stream
          << "main_matrix has a distribution, but the number of rows is not "
          << "equal to the number of degrees of freedom in the problem.";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
    // If the distributions are not the same, then complain
    else if (main_matrix.distribution_built() !=
             mass_matrix.distribution_built())
    {
      std::ostringstream error_stream;
      error_stream << "The distribution of the jacobian and mass matrix must "
                   << "both be setup or both not setup";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Store the old assembly handler
    AssemblyHandler* old_assembly_handler_pt = Assembly_handler_pt;
    // Now setup the eigenproblem handler, pass in the value of the shift
    Assembly_handler_pt = new EigenProblemHandler(shift);

    // Prepare the storage formats.
    Vector<int*> column_or_row_index(2);
    Vector<int*> row_or_column_start(2);
    Vector<double*> value(2);
    Vector<unsigned> nnz(2);
    // Allocate pointer to residuals, although not used in these problems
    Vector<double*> residuals_vectors(0);

    // total number of rows in each matrix
    unsigned nrow = this->ndof();

    // determine the distribution for the jacobian (main matrix)
    // IF the jacobian has distribution setup then use that
    // ELSE determine the distribution based on the
    // distributed_matrix_distribution enum
    LinearAlgebraDistribution* dist_pt = 0;
    if (main_matrix.distribution_built())
    {
      dist_pt = new LinearAlgebraDistribution(main_matrix.distribution_pt());
    }
    else
    {
#ifdef OOMPH_HAS_MPI
      // if problem is only one one processor
      if (Communicator_pt->nproc() == 1)
      {
        dist_pt = new LinearAlgebraDistribution(Communicator_pt, nrow, false);
      }
      // if the problem is not distributed then assemble the matrices with
      // a uniform distributed distribution
      else if (!Problem_has_been_distributed)
      {
        dist_pt = new LinearAlgebraDistribution(Communicator_pt, nrow, true);
      }
      // otherwise the problem is a distributed problem
      else
      {
        switch (Dist_problem_matrix_distribution)
        {
          case Uniform_matrix_distribution:
            dist_pt =
              new LinearAlgebraDistribution(Communicator_pt, nrow, true);
            break;
          case Problem_matrix_distribution:
            dist_pt = new LinearAlgebraDistribution(Dof_distribution_pt);
            break;
          case Default_matrix_distribution:
            LinearAlgebraDistribution* uniform_dist_pt =
              new LinearAlgebraDistribution(Communicator_pt, nrow, true);
            bool use_problem_dist = true;
            unsigned nproc = Communicator_pt->nproc();
            for (unsigned p = 0; p < nproc; p++)
            {
              if ((double)Dof_distribution_pt->nrow_local(p) >
                  ((double)uniform_dist_pt->nrow_local(p)) * 1.1)
              {
                use_problem_dist = false;
              }
            }
            if (use_problem_dist)
            {
              dist_pt = new LinearAlgebraDistribution(Dof_distribution_pt);
            }
            else
            {
              dist_pt = new LinearAlgebraDistribution(uniform_dist_pt);
            }
            delete uniform_dist_pt;
            break;
        }
      }
#else
      dist_pt = new LinearAlgebraDistribution(Communicator_pt, nrow, false);
#endif
    }


    // The matrix is in compressed row format
    bool compressed_row_flag = true;

#ifdef OOMPH_HAS_MPI
    //
    if (Communicator_pt->nproc() == 1)
    {
#endif

      sparse_assemble_row_or_column_compressed(column_or_row_index,
                                               row_or_column_start,
                                               value,
                                               nnz,
                                               residuals_vectors,
                                               compressed_row_flag);

      // The main matrix is the first entry
      main_matrix.build(dist_pt);
      main_matrix.build_without_copy(dist_pt->nrow(),
                                     nnz[0],
                                     value[0],
                                     column_or_row_index[0],
                                     row_or_column_start[0]);
      // The mass matrix is the second entry
      mass_matrix.build(dist_pt);
      mass_matrix.build_without_copy(dist_pt->nrow(),
                                     nnz[1],
                                     value[1],
                                     column_or_row_index[1],
                                     row_or_column_start[1]);
#ifdef OOMPH_HAS_MPI
    }
    else
    {
      if (dist_pt->distributed())
      {
        parallel_sparse_assemble(dist_pt,
                                 column_or_row_index,
                                 row_or_column_start,
                                 value,
                                 nnz,
                                 residuals_vectors);
        // The main matrix is the first entry
        main_matrix.build(dist_pt);
        main_matrix.build_without_copy(dist_pt->nrow(),
                                       nnz[0],
                                       value[0],
                                       column_or_row_index[0],
                                       row_or_column_start[0]);
        // The mass matrix is the second entry
        mass_matrix.build(dist_pt);
        mass_matrix.build_without_copy(dist_pt->nrow(),
                                       nnz[1],
                                       value[1],
                                       column_or_row_index[1],
                                       row_or_column_start[1]);
      }
      else
      {
        LinearAlgebraDistribution* temp_dist_pt =
          new LinearAlgebraDistribution(Communicator_pt, dist_pt->nrow(), true);
        parallel_sparse_assemble(temp_dist_pt,
                                 column_or_row_index,
                                 row_or_column_start,
                                 value,
                                 nnz,
                                 residuals_vectors);
        // The main matrix is the first entry
        main_matrix.build(temp_dist_pt);
        main_matrix.build_without_copy(dist_pt->nrow(),
                                       nnz[0],
                                       value[0],
                                       column_or_row_index[0],
                                       row_or_column_start[0]);
        main_matrix.redistribute(dist_pt);
        // The mass matrix is the second entry
        mass_matrix.build(temp_dist_pt);
        mass_matrix.build_without_copy(dist_pt->nrow(),
                                       nnz[1],
                                       value[1],
                                       column_or_row_index[1],
                                       row_or_column_start[1]);
        mass_matrix.redistribute(dist_pt);
        delete temp_dist_pt;
      }
    }
#endif

    // clean up dist_pt and residuals_vector pt
    delete dist_pt;

    // Delete the eigenproblem handler
    delete Assembly_handler_pt;
    // Reset the assembly handler to the original handler
    Assembly_handler_pt = old_assembly_handler_pt;
  }


  //=======================================================================
  /// Stored the current values of the dofs
  //=======================================================================
  void Problem::store_current_dof_values()
  {
    // If memory has not been allocated, then allocated memory for the saved
    // dofs
    if (Saved_dof_pt == 0)
    {
      Saved_dof_pt = new Vector<double>;
    }

#ifdef OOMPH_HAS_MPI
    // If the problem is distributed I have to do something different
    if (Problem_has_been_distributed)
    {
      // How many entries do we store locally?
      const unsigned n_row_local = Dof_distribution_pt->nrow_local();

      // Resize the vector
      Saved_dof_pt->resize(n_row_local);

      // Back 'em up
      for (unsigned i = 0; i < n_row_local; i++)
      {
        (*Saved_dof_pt)[i] = *(this->Dof_pt[i]);
      }
    }
    // Otherwise just store all the dofs
    else
#endif
    {
      // Find the number of dofs
      unsigned long n_dof = ndof();

      // Resize the vector
      Saved_dof_pt->resize(n_dof);

      // Transfer the values over
      for (unsigned long n = 0; n < n_dof; n++)
      {
        (*Saved_dof_pt)[n] = dof(n);
      }
    }
  }

  //====================================================================
  /// Restore the saved dofs
  //====================================================================
  void Problem::restore_dof_values()
  {
    // Check that we can do this
    if (Saved_dof_pt == 0)
    {
      throw OomphLibError(
        "There are no stored values, use store_current_dof_values()\n",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }


#ifdef OOMPH_HAS_MPI
    // If the problem is distributed I have to do something different
    if (Problem_has_been_distributed)
    {
      // How many entries do we store locally?
      const unsigned n_row_local = Dof_distribution_pt->nrow_local();

      if (Saved_dof_pt->size() != n_row_local)
      {
        throw OomphLibError("The number of stored values is not equal to the "
                            "current number of dofs\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Transfer the values over
      for (unsigned long n = 0; n < n_row_local; n++)
      {
        *(this->Dof_pt[n]) = (*Saved_dof_pt)[n];
      }
    }
    // Otherwise just restore all the dofs
    else
#endif
    {
      // Find the number of dofs
      unsigned long n_dof = ndof();

      if (Saved_dof_pt->size() != n_dof)
      {
        throw OomphLibError("The number of stored values is not equal to the "
                            "current number of dofs\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Transfer the values over
      for (unsigned long n = 0; n < n_dof; n++)
      {
        dof(n) = (*Saved_dof_pt)[n];
      }
    }

    // Delete the memory
    delete Saved_dof_pt;
    Saved_dof_pt = 0;
  }

  //======================================================================
  /// Assign the eigenvector passed to the function to the dofs
  //======================================================================
  void Problem::assign_eigenvector_to_dofs(DoubleVector& eigenvector)
  {
    unsigned long n_dof = ndof();
    // Check that the eigenvector has the correct size
    if (eigenvector.nrow() != n_dof)
    {
      std::ostringstream error_message;
      error_message << "Eigenvector has size " << eigenvector.nrow()
                    << ", not equal to the number of dofs in the problem,"
                    << n_dof << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Loop over the dofs and assign the eigenvector
    for (unsigned long n = 0; n < eigenvector.nrow_local(); n++)
    {
      dof(n) = eigenvector[n];
    }
// Of course we now need to synchronise
#ifdef OOMPH_HAS_MPI
    this->synchronise_all_dofs();
#endif
  }


  //======================================================================
  /// Add the eigenvector passed to the function to the dofs with
  /// magnitude epsilon
  //======================================================================
  void Problem::add_eigenvector_to_dofs(const double& epsilon,
                                        const DoubleVector& eigenvector)
  {
    unsigned long n_dof = ndof();
    // Check that the eigenvector has the correct size
    if (eigenvector.nrow() != n_dof)
    {
      std::ostringstream error_message;
      error_message << "Eigenvector has size " << eigenvector.nrow()
                    << ", not equal to the number of dofs in the problem,"
                    << n_dof << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Loop over the dofs and add the eigenvector
    // Only use local values
    for (unsigned long n = 0; n < eigenvector.nrow_local(); n++)
    {
      dof(n) += epsilon * eigenvector[n];
    }
// Of course we now need to synchronise
#ifdef OOMPH_HAS_MPI
    this->synchronise_all_dofs();
#endif
  }


  //================================================================
  /// General Newton solver. Requires only a convergence tolerance.
  /// The linear solver takes a pointer to the problem (which defines
  /// the Jacobian \b J and the residual Vector \b r) and returns
  /// the solution \b x of the system
  /// \f[ {\bf J} {\bf x} = - \bf{r} \f].
  //================================================================
  void Problem::newton_solve()
  {
    // Initialise timers
    double total_linear_solver_time = 0.0;
    double t_start = TimingHelpers::timer();
    Max_res.clear();

    // Find total number of dofs
    unsigned long n_dofs = ndof();

    // Set up the Vector to hold the solution
    DoubleVector dx;

    //-----Variables for the globally convergent Newton method------

    // Set up the vector to hold the gradient
    DoubleVector gradient;

    // Other variables
    double half_residual_squared = 0.0;
    double max_step = 0.0;

    //--------------------------------------------------------------

    // Set the counter
    unsigned count = 0;
    // Set the loop flag
    unsigned LOOP_FLAG = 1;

    if (Use_globally_convergent_newton_method)
    {
#ifdef OOMPH_HAS_MPI
      // Break if running in parallel
      if (MPI_Helpers::mpi_has_been_initialised())
      {
        std::ostringstream error_stream;
        error_stream << "Globally convergent Newton method has not been "
                     << "implemented in parallel yet!" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Get gradient
      Linear_solver_pt->enable_computation_of_gradient();
      // Reset the gradient (clear it), since the number of dofs and
      // hence the size of the DoubleVector might have changed
      Linear_solver_pt->reset_gradient();
    }

    // Update anything that needs updating
    actions_before_newton_solve();

    // Reset number of Newton iterations taken
    Nnewton_iter_taken = 0;

    // Now do the Newton loop
    do
    {
      count++;

      // Do any updates that are required
      actions_before_newton_step();


      // No degrees of freedom? What are you solving for?
      if (n_dofs == 0)
      {
        oomph_info << std::endl << std::endl << std::endl;
        oomph_info << "This is a bit bizarre: The problem has no dofs."
                   << std::endl;
        oomph_info
          << "I'll just return from the Newton solver without doing anything."
          << std::endl;

        // Do any updates that would have been performed
        actions_before_newton_convergence_check();
        actions_after_newton_step();
        actions_before_newton_convergence_check();
        actions_after_newton_solve();

        oomph_info << "I hope this is what you intended me to do..."
                   << std::endl;
        oomph_info << std::endl
                   << "Note: All actions_...() functions were called"
                   << std::endl;
        oomph_info << std::endl << "      before returning." << std::endl;
        oomph_info << std::endl << std::endl << std::endl;
        return;
      }

      // Calculate initial residuals
      if (count == 1)
      {
        // Is the problem nonlinear? If not ignore the pre-iteration
        // convergence check.
        if (Problem_is_nonlinear)
        {
#ifdef OOMPH_HAS_MPI
          // Synchronise the solution on different processors (on each submesh)
          this->synchronise_all_dofs();
#endif

          actions_before_newton_convergence_check();
          dx.clear();
          get_residuals(dx);

          // Get half of squared residual and find maximum step length
          // for step length control
          if (Use_globally_convergent_newton_method)
          {
            half_residual_squared = 0.0;
            double sum = 0.0;
            for (unsigned i = 0; i < n_dofs; i++)
            {
              sum += (*Dof_pt[i]) * (*Dof_pt[i]);
              half_residual_squared += dx[i] * dx[i];
            }
            half_residual_squared *= 0.5;
            max_step = 100.0 * std::max(sqrt(sum), double(n_dofs));
          }

          // Get maximum residuals
          double maxres = dx.max();
          Max_res.push_back(maxres);

          if (!Shut_up_in_newton_solve)
          {
            // Let's output the residuals
            // unsigned n_row_local = dx.distribution_pt()->nrow_local();
            // unsigned first_row = dx.distribution_pt()->first_row();
            // for(unsigned n=0;n<n_row_local;n++)
            //{
            // oomph_info << "residual: " << n + first_row << " " << dx[n] <<
            // "\n";
            //}

            oomph_info << "\nInitial Maximum residuals " << maxres << std::endl;
          }

          if ((maxres < Newton_solver_tolerance) &&
              (!Always_take_one_newton_step))
          {
            LOOP_FLAG = 0;
            continue;
          }
        }
        else
        {
          if (!Shut_up_in_newton_solve)
          {
            oomph_info
              << "Linear problem -- convergence in one iteration assumed."
              << std::endl;
          }
        }
      }


      // Increment number of Newton iterations taken
      Nnewton_iter_taken++;

      // Initialise timer for linear solver
      double t_solver_start = TimingHelpers::timer();

      // Now do the linear solve -- recycling Jacobian if requested
      if (Jacobian_reuse_is_enabled && Jacobian_has_been_computed)
      {
        if (!Shut_up_in_newton_solve)
        {
          oomph_info << "Not recomputing Jacobian! " << std::endl;
        }

        // If we're doing the first iteration and the problem is nonlinear,
        // the residuals have already been computed above during the
        // initial convergence check. Otherwise compute them here.
        if ((count != 1) || (!Problem_is_nonlinear)) get_residuals(dx);

        // Backup residuals
        DoubleVector resid(dx);

        // Resolve
        Linear_solver_pt->resolve(resid, dx);
      }
      else
      {
        if (Jacobian_reuse_is_enabled)
        {
          if (!Shut_up_in_newton_solve)
          {
            oomph_info << "Enabling resolve" << std::endl;
          }
          Linear_solver_pt->enable_resolve();
        }
        Linear_solver_pt->solve(this, dx);
        Jacobian_has_been_computed = true;
      }

      // End of linear solver
      double t_solver_end = TimingHelpers::timer();
      total_linear_solver_time += t_solver_end - t_solver_start;

      if (!Shut_up_in_newton_solve)
      {
        oomph_info << std::endl;
        oomph_info << "Time for linear solver ( ndof = " << n_dofs
                   << " ) [sec]: " << t_solver_end - t_solver_start << std::endl
                   << std::endl;
      }

      // Subtract the new values from the true dofs
      dx.redistribute(Dof_distribution_pt);
      double* dx_pt = dx.values_pt();
      unsigned ndof_local = Dof_distribution_pt->nrow_local();

      if (Use_globally_convergent_newton_method)
      {
        // Get the gradient
        Linear_solver_pt->get_gradient(gradient);

        for (unsigned i = 0; i < ndof_local; i++)
        {
          dx_pt[i] *= -1.0;
        }

        // Update with steplength control
        Vector<double> unknowns_old(ndof_local);

        for (unsigned i = 0; i < ndof_local; i++)
        {
          unknowns_old[i] = *Dof_pt[i];
        }

        double half_residual_squared_old = half_residual_squared;
        globally_convergent_line_search(unknowns_old,
                                        half_residual_squared_old,
                                        gradient,
                                        dx,
                                        half_residual_squared,
                                        max_step);
      }
      // direct Newton update
      else
      {
        for (unsigned l = 0; l < ndof_local; l++)
        {
          *Dof_pt[l] -= Relaxation_factor * dx_pt[l];
        }
      }
#ifdef OOMPH_HAS_MPI
      // Synchronise the solution on different processors (on each submesh)
      this->synchronise_all_dofs();
#endif

      // Do any updates that are required
      actions_after_newton_step();
      actions_before_newton_convergence_check();

      // Maximum residuals
      double maxres = 0.0;
      // If the user has declared that the Problem is linear
      // we ignore the convergence check
      if (Problem_is_nonlinear)
      {
        // Get the maximum residuals
        // maxres = std::fabs(*std::max_element(dx.begin(),dx.end(),
        //                                    AbsCmp<double>()));
        // oomph_info << "Maxres correction " << maxres << "\n";

        // Calculate the new residuals
        dx.clear();
        get_residuals(dx);

        // Get the maximum residuals
        maxres = dx.max();
        Max_res.push_back(maxres);

        if (!Shut_up_in_newton_solve)
        {
          oomph_info << "Newton Step " << count << ": Maximum residuals "
                     << maxres << std::endl
                     << std::endl;
        }
      }

      // If we have converged jump straight to the test at the end of the loop
      if (maxres < Newton_solver_tolerance)
      {
        LOOP_FLAG = 0;
        continue;
      }

      // This section will not be reached if we have converged already
      // If the maximum number of residuals is too high or the maximum number
      // of iterations has been reached
      if ((maxres > Max_residuals) || (count == Max_newton_iterations))
      {
        // Print a warning -- regardless of what the throw does
        if (maxres > Max_residuals)
        {
          oomph_info << "Max. residual (" << Max_residuals
                     << ") has been exceeded in Newton solver." << std::endl;
        }
        if (count == Max_newton_iterations)
        {
          oomph_info << "Reached max. number of iterations ("
                     << Max_newton_iterations << ") in Newton solver."
                     << std::endl;
        }
        // Now throw...
        throw NewtonSolverError(count, maxres);
      }

    } while (LOOP_FLAG);

    // Now update anything that needs updating
    actions_after_newton_solve();

    // Finalise/doc timings
    if (!Shut_up_in_newton_solve)
    {
      oomph_info << std::endl;
      oomph_info << "Total time for linear solver (ndof=" << n_dofs
                 << ") [sec]: " << total_linear_solver_time << std::endl;
    }

    double t_end = TimingHelpers::timer();
    double total_time = t_end - t_start;

    if (!Shut_up_in_newton_solve)
    {
      oomph_info << "Total time for Newton solver (ndof=" << n_dofs
                 << ") [sec]: " << total_time << std::endl;
    }
    if (total_time > 0.0)
    {
      if (!Shut_up_in_newton_solve)
      {
        oomph_info << "Time outside linear solver        : "
                   << (total_time - total_linear_solver_time) / total_time *
                        100.0
                   << " %" << std::endl;
      }
    }
    else
    {
      if (!Shut_up_in_newton_solve)
      {
        oomph_info << "Time outside linear solver        : "
                   << "[too fast]" << std::endl;
      }
    }
    if (!Shut_up_in_newton_solve) oomph_info << std::endl;
  }

  //========================================================================
  /// Helper function for the globally convergent Newton solver
  //========================================================================
  void Problem::globally_convergent_line_search(
    const Vector<double>& x_old,
    const double& half_residual_squared_old,
    DoubleVector& gradient,
    DoubleVector& newton_dir,
    double& half_residual_squared,
    const double& stpmax)
  {
    const double min_fct_decrease = 1.0e-4;
    double convergence_tol_on_x = 1.0e-16;
    double f_aux = 0.0;
    double lambda_aux = 0.0;
    double proposed_lambda;
    unsigned long n_dof = ndof();
    double sum = 0.0;
    for (unsigned i = 0; i < n_dof; i++)
    {
      sum += newton_dir[i] * newton_dir[i];
    }
    sum = sqrt(sum);
    if (sum > stpmax)
    {
      for (unsigned i = 0; i < n_dof; i++)
      {
        newton_dir[i] *= stpmax / sum;
      }
    }
    double slope = 0.0;
    for (unsigned i = 0; i < n_dof; i++)
    {
      slope += gradient[i] * newton_dir[i];
    }
    if (slope >= 0.0)
    {
      std::ostringstream warn_message;
      warn_message << "WARNING: Non-negative slope, probably due to a "
                   << " roundoff \nproblem in the linesearch: slope=" << slope
                   << "\n";
      OomphLibWarning(warn_message.str(),
                      "Problem::globally_convergent_line_search()",
                      OOMPH_EXCEPTION_LOCATION);
    }
    double test = 0.0;
    for (unsigned i = 0; i < n_dof; i++)
    {
      double temp =
        std::fabs(newton_dir[i]) / std::max(std::fabs(x_old[i]), 1.0);
      if (temp > test) test = temp;
    }
    double lambda_min = convergence_tol_on_x / test;
    double lambda = 1.0;
    while (true)
    {
      for (unsigned i = 0; i < n_dof; i++)
      {
        *Dof_pt[i] = x_old[i] + lambda * newton_dir[i];
      }

      // Evaluate current residuals
      DoubleVector residuals;
      get_residuals(residuals);
      half_residual_squared = 0.0;
      for (unsigned i = 0; i < n_dof; i++)
      {
        half_residual_squared += residuals[i] * residuals[i];
      }
      half_residual_squared *= 0.5;

      if (lambda < lambda_min)
      {
        for (unsigned i = 0; i < n_dof; i++) *Dof_pt[i] = x_old[i];

        std::ostringstream warn_message;
        warn_message << "WARNING: Line search converged on x only!\n";
        OomphLibWarning(warn_message.str(),
                        "Problem::globally_convergent_line_search()",
                        OOMPH_EXCEPTION_LOCATION);
        return;
      }
      else if (half_residual_squared <=
               half_residual_squared_old + min_fct_decrease * lambda * slope)
      {
        oomph_info << "Returning from linesearch with lambda=" << lambda
                   << std::endl;
        return;
      }
      else
      {
        if (lambda == 1.0)
        {
          proposed_lambda =
            -slope /
            (2.0 * (half_residual_squared - half_residual_squared_old - slope));
        }
        else
        {
          double r1 =
            half_residual_squared - half_residual_squared_old - lambda * slope;
          double r2 = f_aux - half_residual_squared_old - lambda_aux * slope;
          double a_poly =
            (r1 / (lambda * lambda) - r2 / (lambda_aux * lambda_aux)) /
            (lambda - lambda_aux);
          double b_poly = (-lambda_aux * r1 / (lambda * lambda) +
                           lambda * r2 / (lambda_aux * lambda_aux)) /
                          (lambda - lambda_aux);
          if (a_poly == 0.0)
          {
            proposed_lambda = -slope / (2.0 * b_poly);
          }
          else
          {
            double discriminant = b_poly * b_poly - 3.0 * a_poly * slope;
            if (discriminant < 0.0)
            {
              proposed_lambda = 0.5 * lambda;
            }
            else if (b_poly <= 0.0)
            {
              proposed_lambda = (-b_poly + sqrt(discriminant)) / (3.0 * a_poly);
            }
            else
            {
              proposed_lambda = -slope / (b_poly + sqrt(discriminant));
            }
          }
          if (proposed_lambda > 0.5 * lambda)
          {
            proposed_lambda = 0.5 * lambda;
          }
        }
      }
      lambda_aux = lambda;
      f_aux = half_residual_squared;
      lambda = std::max(proposed_lambda, 0.1 * lambda);
    }
  }


  //========================================================================
  /// Solve a steady problem, in the context of an overall unsteady problem.
  /// This is achieved by setting the weights in the timesteppers to be zero
  /// which has the effect of rendering them steady timesteppers
  /// The optional argument max_adapt specifies the max. number of
  /// adaptations of all refineable submeshes are performed to
  /// achieve the the error targets specified in the refineable submeshes.
  //========================================================================
  void Problem::steady_newton_solve(unsigned const& max_adapt)
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

    try
    {
      // Solve the non-linear problem with Newton's method
      if (max_adapt == 0)
      {
        newton_solve();
      }
      else
      {
        newton_solve(max_adapt);
      }
    }
    // Catch any exceptions thrown in the Newton solver
    catch (NewtonSolverError& error)
    {
      oomph_info << std::endl
                 << "USER-DEFINED ERROR IN NEWTON SOLVER " << std::endl;
      // Check whether it's the linear solver
      if (error.linear_solver_error)
      {
        oomph_info << "ERROR IN THE LINEAR SOLVER" << std::endl;
      }
      // Check to see whether we have reached Max_iterations
      else if (error.iterations == Max_newton_iterations)
      {
        oomph_info << "MAXIMUM NUMBER OF ITERATIONS (" << error.iterations
                   << ") REACHED WITHOUT CONVERGENCE " << std::endl;
      }
      // If not, it must be that we have exceeded the maximum residuals
      else
      {
        oomph_info << "MAXIMUM RESIDUALS: " << error.maxres
                   << " EXCEEDS PREDEFINED MAXIMUM " << Max_residuals
                   << std::endl;
      }

      // Die horribly!!
      std::ostringstream error_stream;
      error_stream << "Error occured in Newton solver. " << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
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
  }

  //===========================================================================
  /// Perform a basic continuation step using Newton's method. The governing
  /// parameter of the problem is passed as a pointer to the routine. The
  /// number of Newton steps taken is returned
  //==========================================================================
  unsigned Problem::newton_solve_continuation(double* const& parameter_pt)
  {
    // Set up memory for z
    // unsigned long n_dofs = ndof();
    // LinearAlgebraDistribution dist(Communicator_pt,n_dofs,false);
    // DoubleVector z(&dist,0.0);
    DoubleVector z;
    // Call the solver
    return newton_solve_continuation(parameter_pt, z);
  }


  //===================================================================
  /// This function performs a basic continuation step using the Newton method.
  /// The number of Newton steps taken is returned, to be used in any
  /// external step-size control routines.
  /// The governing parameter of the problem is passed as a pointer to the
  /// routine, as is the sign of the Jacobian and a Vector in which
  /// to store the derivatives wrt the parameter, if required.
  //==================================================================
  unsigned Problem::newton_solve_continuation(double* const& parameter_pt,
                                              DoubleVector& z)
  {
    // Find the total number of dofs
    // unsigned long n_dofs = ndof();

    // Find the local number of dofs
    unsigned ndof_local = Dof_distribution_pt->nrow_local();

    // create the distribution (not distributed)
    // LinearAlgebraDistribution dist(this->communicator_pt(),n_dofs,false);

    // Assign memory for solutions of the equations
    // DoubleVector y(&dist,0.0);
    DoubleVector y;

    // Assign memory for the dot products of the uderivatives and y and z
    double uderiv_dot_y = 0.0, uderiv_dot_z = 0.0;
    // Set and initialise the counter
    unsigned count = 0;
    // Set the loop flag
    unsigned LOOP_FLAG = 1;

    // Update anything that needs updating
    actions_before_newton_solve();

    // Check the arc-length constraint
    double arc_length_constraint_residual = 0.0;

    // Are we storing the matrix in the linear solve
    bool enable_resolve = Linear_solver_pt->is_resolve_enabled();

    // For this problem, we must store the residuals
    Linear_solver_pt->enable_resolve();

    // Now do the Newton loop
    do
    {
      count++;

      // Do any updates that are required
      actions_before_newton_step();

      // Calculate initial residuals
      if (count == 1)
      {
#ifdef OOMPH_HAS_MPI
        // Synchronise the solution on different processors (on each submesh)
        this->synchronise_all_dofs();
#endif

        actions_before_newton_convergence_check();
        y.clear();
        get_residuals(y);
        // Get maximum residuals, using our own abscmp function
        double maxres = y.max();

        // Assemble the residuals for the arc-length step
        arc_length_constraint_residual = 0.0;
        // Add the variables
        for (unsigned long l = 0; l < ndof_local; l++)
        {
          arc_length_constraint_residual +=
            dof_derivative(l) * (*Dof_pt[l] - dof_current(l));
        }

        // Now reduce if we have been distributed
#ifdef OOMPH_HAS_MPI
        double arc_length_cons_res2 = arc_length_constraint_residual;
        if ((Dof_distribution_pt->distributed()) &&
            (Dof_distribution_pt->communicator_pt()->nproc() > 1))
        {
          MPI_Allreduce(&arc_length_constraint_residual,
                        &arc_length_cons_res2,
                        1,
                        MPI_DOUBLE,
                        MPI_SUM,
                        Dof_distribution_pt->communicator_pt()->mpi_comm());
        }
        arc_length_constraint_residual = arc_length_cons_res2;
#endif

        arc_length_constraint_residual *= Theta_squared;
        arc_length_constraint_residual +=
          Parameter_derivative * (*parameter_pt - Parameter_current) -
          Ds_current;

        // Is it the max
        if (std::fabs(arc_length_constraint_residual) > maxres)
        {
          maxres = std::fabs(arc_length_constraint_residual);
        }

        // Find the max
        if (!Shut_up_in_newton_solve)
        {
          oomph_info << "Initial Maximum residuals " << maxres << std::endl;
        }

        // If we are below the Tolerance, then return immediately
        if ((maxres < Newton_solver_tolerance) &&
            (!Always_take_one_newton_step))
        {
          LOOP_FLAG = 0;
          count = 0;
          continue;
        }
      }

      // If it's the block hopf solver we need to solve for both rhs's
      // simultaneously. This is because the block decomposition involves
      // solves with two different matrices and storing both at once to
      // allow general resolves would be more expensive than necessary.
      if (dynamic_cast<BlockHopfLinearSolver*>(Linear_solver_pt))
      {
        // Get the vector dresiduals/dparameter
        z.clear();
        get_derivative_wrt_global_parameter(parameter_pt, z);

        // Copy rhs vector into local storage so it doesn't get overwritten
        // if the linear solver decides to initialise the solution vector, say,
        // which it's quite entitled to do!
        DoubleVector input_z(z);

        // Solve the system for the two right-hand sides.
        dynamic_cast<BlockHopfLinearSolver*>(Linear_solver_pt)
          ->solve_for_two_rhs(this, y, input_z, z);
      }
      // Otherwise
      else
      {
        // Solve the standard problem
        Linear_solver_pt->solve(this, y);

        // Get the vector dresiduals/dparameter
        z.clear();
        get_derivative_wrt_global_parameter(parameter_pt, z);

        // Copy rhs vector into local storage so it doesn't get overwritten
        // if the linear solver decides to initialise the solution vector, say,
        // which it's quite entitled to do!
        DoubleVector input_z(z);

        // Redistribute the RHS to match the linear solver
        // input_z.redistribute(Linear_solver_pt->distribution_pt());
        // Do not clear z because we assume that it has dR/dparam
        z.clear();
        // Now resolve the system with the new RHS
        Linear_solver_pt->resolve(input_z, z);
      }

      // Redistribute the results into the natural distribution
      y.redistribute(Dof_distribution_pt);
      z.redistribute(Dof_distribution_pt);

      // Now we need to calculate dparam, for which we must calculate the
      // dot product of the derivatives and y and z
      // Reset these values to zero
      uderiv_dot_y = 0.0;
      uderiv_dot_z = 0.0;
      // Now calculate the dot products of the derivative and the solutions
      // of the linear system
      // Cache pointers to the data in the distributed vectors
      double* const y_pt = y.values_pt();
      double* const z_pt = z.values_pt();
      for (unsigned long l = 0; l < ndof_local; l++)
      {
        uderiv_dot_y += dof_derivative(l) * y_pt[l];
        uderiv_dot_z += dof_derivative(l) * z_pt[l];
      }

      // Now reduce if we have been distributed
#ifdef OOMPH_HAS_MPI
      // Create send and receive arrays of size two
      double uderiv_dot[2];
      double uderiv_dot2[2];
      uderiv_dot[0] = uderiv_dot_y;
      uderiv_dot[1] = uderiv_dot_z;
      uderiv_dot2[0] = uderiv_dot_y;
      uderiv_dot2[1] = uderiv_dot_z;
      // Now reduce both together
      if ((Dof_distribution_pt->distributed()) &&
          (Dof_distribution_pt->communicator_pt()->nproc() > 1))
      {
        MPI_Allreduce(uderiv_dot,
                      uderiv_dot2,
                      2,
                      MPI_DOUBLE,
                      MPI_SUM,
                      Dof_distribution_pt->communicator_pt()->mpi_comm());
      }
      uderiv_dot_y = uderiv_dot2[0];
      uderiv_dot_z = uderiv_dot2[1];
#endif

      // Now scale the results
      uderiv_dot_y *= Theta_squared;
      uderiv_dot_z *= Theta_squared;

      // The set the change in the parameter, given by the pseudo-arclength
      // equation. Note that here we are assuming that the arc-length
      // equation is always exactly zero,
      // which seems to work OK, and saves on some storage.
      // In fact, it's more subtle than that. If we include this
      // proper residual then we will have to solve the eigenproblem.
      // This will make the solver more robust and *should* be done
      // ... at some point.
      double dparam = (arc_length_constraint_residual - uderiv_dot_y) /
                      (Parameter_derivative - uderiv_dot_z);

      // Set the new value of the parameter
      *parameter_pt -= dparam;

      // Update the values of the other degrees of freedom
      for (unsigned long l = 0; l < ndof_local; l++)
      {
        *Dof_pt[l] -= y_pt[l] - dparam * z_pt[l];
      }

      // Calculate the new residuals
#ifdef OOMPH_HAS_MPI
      // Synchronise the solution on different processors (on each submesh)
      this->synchronise_all_dofs();
#endif

      // Do any updates that are required
      actions_after_newton_step();
      actions_before_newton_convergence_check();

      y.clear();
      get_residuals(y);

      // Get the maximum residuals
      double maxres = y.max();

      // Assemble the residuals for the arc-length step
      arc_length_constraint_residual = 0.0;
      // Add the variables
      for (unsigned long l = 0; l < ndof_local; l++)
      {
        arc_length_constraint_residual +=
          dof_derivative(l) * (*Dof_pt[l] - dof_current(l));
      }

      // Now reduce if we have been distributed
#ifdef OOMPH_HAS_MPI
      double arc_length_cons_res2 = arc_length_constraint_residual;
      if ((Dof_distribution_pt->distributed()) &&
          (Dof_distribution_pt->communicator_pt()->nproc() > 1))
      {
        MPI_Allreduce(&arc_length_constraint_residual,
                      &arc_length_cons_res2,
                      1,
                      MPI_DOUBLE,
                      MPI_SUM,
                      Dof_distribution_pt->communicator_pt()->mpi_comm());
      }
      arc_length_constraint_residual = arc_length_cons_res2;
#endif

      arc_length_constraint_residual *= Theta_squared;
      arc_length_constraint_residual +=
        Parameter_derivative * (*parameter_pt - Parameter_current) - Ds_current;

      // Is it the max
      if (std::fabs(arc_length_constraint_residual) > maxres)
      {
        maxres = std::fabs(arc_length_constraint_residual);
      }

      if (!Shut_up_in_newton_solve)
      {
        oomph_info << "Continuation Step " << count << ":  Maximum residuals "
                   << maxres << "\n";
      }

      // If we have converged jump straight to the test at the end of the loop
      if (maxres < Newton_solver_tolerance)
      {
        LOOP_FLAG = 0;
        continue;
      }

      // This section will not be reached if we have converged already
      // If the maximum number of residuals is too high or the maximum number
      // of iterations has been reached
      if ((maxres > Max_residuals) || (count == Max_newton_iterations))
      {
        throw NewtonSolverError(count, maxres);
      }

    } while (LOOP_FLAG);

    // Now update anything that needs updating
    actions_after_newton_solve();

    // Reset the storage of the matrix on the linear solver to what it was
    // on entry to this routine
    if (enable_resolve)
    {
      Linear_solver_pt->enable_resolve();
    }
    else
    {
      Linear_solver_pt->disable_resolve();
    }

    // Return the number of Newton Steps taken
    return count;
  }

  //=========================================================================
  /// A function to calculate the derivatives wrt the arc-length. This version
  /// of the function actually does a linear solve so that the derivatives
  /// are calculated "exactly" rather than using the values at the Newton
  /// step just before convergence. This is only necessary in spatially adaptive
  /// problems, in which the number of degrees of freedom changes and so
  /// the appropriate derivatives must be calculated for the new variables.
  //=========================================================================
  void Problem::calculate_continuation_derivatives(double* const& parameter_pt)
  {
    // Find the number of degrees of freedom in the problem
    const unsigned long n_dofs = ndof();

    // create a non-distributed z vector
    LinearAlgebraDistribution dist(Communicator_pt, n_dofs, false);

    // Assign memory for solutions of the equations
    DoubleVector z(&dist, 0.0);

    // If it's the block hopf solver need to solve for both RHS
    // at once, but this would all be alleviated if we have the solve
    // for the non-residuals RHS.
    if (dynamic_cast<BlockHopfLinearSolver*>(Linear_solver_pt))
    {
      // Get the vector dresiduals/dparameter
      get_derivative_wrt_global_parameter(parameter_pt, z);

      // Copy rhs vector into local storage so it doesn't get overwritten
      // if the linear solver decides to initialise the solution vector, say,
      // which it's quite entitled to do!
      DoubleVector dummy(&dist, 0.0), input_z(z);

      // Solve for the two RHSs
      dynamic_cast<BlockHopfLinearSolver*>(Linear_solver_pt)
        ->solve_for_two_rhs(this, dummy, input_z, z);
    }
    // Otherwise we can use the normal resolve
    else
    {
      // Save the status before entry to this routine
      bool enable_resolve = Linear_solver_pt->is_resolve_enabled();

      // We need to do resolves
      Linear_solver_pt->enable_resolve();

      // Solve the standard problem, we only want to make sure that
      // we factorise the matrix, if it has not been factorised. We shall
      // ignore the return value of z.
      Linear_solver_pt->solve(this, z);

      // Get the vector dresiduals/dparameter
      get_derivative_wrt_global_parameter(parameter_pt, z);


      // Copy rhs vector into local storage so it doesn't get overwritten
      // if the linear solver decides to initialise the solution vector, say,
      // which it's quite entitled to do!
      DoubleVector input_z(z);

      // Now resolve the system with the new RHS and overwrite the solution
      Linear_solver_pt->resolve(input_z, z);

      // Restore the storage status of the linear solver
      if (enable_resolve)
      {
        Linear_solver_pt->enable_resolve();
      }
      else
      {
        Linear_solver_pt->disable_resolve();
      }
    }

    // Now, we can calculate the derivatives, etc
    calculate_continuation_derivatives(z);
  }

  //=======================================================================
  /// A function to calculate the derivatives with respect to the arc-length
  /// required for continuation. The arguments is the solution of the
  /// linear system,
  /// Jz = dR/dparameter, that gives du/dparameter and the direction
  /// output from the newton_solve_continuation function. The derivatives
  /// are stored in the ContinuationParameters namespace.
  //===================================================================
  void Problem::calculate_continuation_derivatives(const DoubleVector& z)
  {
    // Calculate the continuation derivatives
    calculate_continuation_derivatives_helper(z);

    // Scale the value of theta if the control flag is set
    if (Scale_arc_length)
    {
      // Don't divide by zero!
      if (Parameter_derivative != 1.0)
      {
        Theta_squared *= (Parameter_derivative * Parameter_derivative /
                          Desired_proportion_of_arc_length) *
                         ((1.0 - Desired_proportion_of_arc_length) /
                          (1.0 - Parameter_derivative * Parameter_derivative));

        // Recalculate the continuation derivatives with the new scaled values
        calculate_continuation_derivatives_helper(z);
      }
    }
  }

  //=======================================================================
  /// A function to calculate the derivatives with respect to the arc-length
  /// required for continuation using finite differences.
  //===================================================================
  void Problem::calculate_continuation_derivatives_fd(
    double* const& parameter_pt)
  {
    // Calculate the continuation derivatives
    calculate_continuation_derivatives_fd_helper(parameter_pt);

    // Scale the value of theta if the control flag is set
    if (Scale_arc_length)
    {
      // Don't divide by zero!
      if (Parameter_derivative != 1.0)
      {
        Theta_squared *= (Parameter_derivative * Parameter_derivative /
                          Desired_proportion_of_arc_length) *
                         ((1.0 - Desired_proportion_of_arc_length) /
                          (1.0 - Parameter_derivative * Parameter_derivative));

        // Recalculate the continuation derivatives with the new scaled values
        calculate_continuation_derivatives_fd_helper(parameter_pt);
      }
    }
  }

  //======================================================================
  /// Function that returns a boolean flag to indicate whether the pointer
  /// parameter_pt refers to memory that is a value in a Data object used
  /// within the problem
  //======================================================================
  bool Problem::does_pointer_correspond_to_problem_data(
    double* const& parameter_pt)
  {
    // Firstly check the global data
    const unsigned n_global = Global_data_pt.size();
    for (unsigned i = 0; i < n_global; ++i)
    {
      // If we find it then return true
      if (Global_data_pt[i]->does_pointer_correspond_to_value(parameter_pt))
      {
        return true;
      }
    }

    // If we find the pointer in the mesh data return true
    if (Mesh_pt->does_pointer_correspond_to_mesh_data(parameter_pt))
    {
      return true;
    }

    // Loop over the submeshes to handle the case of spine data
    const unsigned n_sub_mesh = this->nsub_mesh();
    // If there is only one mesh
    if (n_sub_mesh == 0)
    {
      if (SpineMesh* const spine_mesh_pt = dynamic_cast<SpineMesh*>(Mesh_pt))
      {
        if (spine_mesh_pt->does_pointer_correspond_to_spine_data(parameter_pt))
        {
          return true;
        }
      }
    }
    // Otherwise loop over the sub meshes
    else
    {
      // Assign global equation numbers first
      for (unsigned i = 0; i < n_sub_mesh; i++)
      {
        if (SpineMesh* const spine_mesh_pt =
              dynamic_cast<SpineMesh*>(Sub_mesh_pt[i]))
        {
          if (spine_mesh_pt->does_pointer_correspond_to_spine_data(
                parameter_pt))
          {
            return true;
          }
        }
      }
    }

    // If we have got here then the data is not stored in the problem, so return
    // false
    return false;
  }


  //=======================================================================
  /// A private helper function to
  /// calculate the derivatives with respect to the arc-length
  /// required for continuation. The arguments is the solution of the
  /// linear system,
  /// Jz = dR/dparameter, that gives du/dparameter and the direction
  /// output from the newton_solve_continuation function. The derivatives
  /// are stored in the ContinuationParameters namespace.
  //===================================================================
  void Problem::calculate_continuation_derivatives_helper(const DoubleVector& z)
  {
    // Find the number of degrees of freedom in the problem
    // unsigned long n_dofs = ndof();
    // Find the number of local dofs in the problem
    const unsigned long ndof_local = Dof_distribution_pt->nrow_local();

    // Work out the continuation direction
    // The idea is that (du/ds)_{old} . (du/ds)_{new} >= 0
    // if the direction is to remain the same.
    // du/ds_{new} = [dlambda/ds; du/ds] = [dlambda/ds ; - dlambda/ds z]
    // so (du/ds)_{new} . (du/ds)_{old}
    // = dlambda/ds [1 ; - z] . [ Parameter_derivative ; Dof_derivatives]
    // = dlambda/ds (Parameter_derivative - Dof_derivative . z)

    // Create a local copy of z that can be redistributed without breaking
    // the constness of z
    DoubleVector local_z(z);

    // Redistribute z so that it has the (natural) dof distribution
    local_z.redistribute(Dof_distribution_pt);

    // Calculate the local contribution to the Continuation direction
    Continuation_direction = 0.0;
    // Cache the pointer to z
    double* const local_z_pt = local_z.values_pt();
    for (unsigned long l = 0; l < ndof_local; l++)
    {
      Continuation_direction -= dof_derivative(l) * local_z_pt[l];
    }

    // Now reduce if we have been distributed
#ifdef OOMPH_HAS_MPI
    double cont_dir2 = Continuation_direction;
    if ((Dof_distribution_pt->distributed()) &&
        (Dof_distribution_pt->communicator_pt()->nproc() > 1))
    {
      MPI_Allreduce(&Continuation_direction,
                    &cont_dir2,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    Dof_distribution_pt->communicator_pt()->mpi_comm());
    }
    Continuation_direction = cont_dir2;
#endif

    // Add parameter derivative
    Continuation_direction += Parameter_derivative;

    // Calculate the magnitude of the du/ds Vector

    // Note that actually, we are usually approximating by using the value at
    // newton step just before convergence, which saves one additional
    // Newton solve.

    // First calculate the magnitude of du/dparameter, chi
    double chi = local_z.dot(local_z);

    // Calculate the current derivative of the parameter wrt the arc-length
    Parameter_derivative = 1.0 / sqrt(1.0 + Theta_squared * chi);

    // If the dot product of the current derivative wrt the Direction
    // is less than zero, switch the sign of the derivative to ensure
    // smooth continuation
    if (Parameter_derivative * Continuation_direction < 0.0)
    {
      Parameter_derivative *= -1.0;
    }

    // Resize the derivatives array, if necessary
    if (!Use_continuation_timestepper)
    {
      if (Dof_derivative.size() != ndof_local)
      {
        Dof_derivative.resize(ndof_local, 0.0);
      }
    }
    // Calculate the new derivatives wrt the arc-length
    for (unsigned long l = 0; l < ndof_local; l++)
    {
      // This comes from the formulation J u_dot + dr/dlambda  lambda_dot = 0
      // on the curve and then it follows that.
      dof_derivative(l) = -Parameter_derivative * local_z_pt[l];
    }
  }

  //=======================================================================
  /// A private helper function to
  /// calculate the derivatives with respect to the arc-length
  /// required for continuation using finite differences.
  //===================================================================
  void Problem::calculate_continuation_derivatives_fd_helper(
    double* const& parameter_pt)
  {
    // Find the number of values
    // const unsigned long n_dofs = this->ndof();
    // Find the number of local dofs in the problem
    const unsigned long ndof_local = Dof_distribution_pt->nrow_local();

    // Temporary storage for the finite-difference approximation to the helper
    Vector<double> z(ndof_local);
    double length = 0.0;
    // Calculate the change in values and contribution to total length
    for (unsigned long l = 0; l < ndof_local; l++)
    {
      z[l] = (*Dof_pt[l] - Dof_current[l]) / Ds_current;
      length += Theta_squared * z[l] * z[l];
    }

    // Reduce if parallel
#ifdef OOMPH_HAS_MPI
    double length2 = length;
    if ((Dof_distribution_pt->distributed()) &&
        (Dof_distribution_pt->communicator_pt()->nproc() > 1))
    {
      MPI_Allreduce(&length,
                    &length2,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    Dof_distribution_pt->communicator_pt()->mpi_comm());
    }
    length = length2;
#endif

    // Calculate change in parameter
    double Z = (*parameter_pt - Parameter_current) / Ds_current;
    length += Z * Z;

    // Scale the approximations to the derivatives
    length = sqrt(length);
    for (unsigned long l = 0; l < ndof_local; l++)
    {
      dof_derivative(l) = z[l] / length;
    }
    Parameter_derivative = Z / length;
  }


  /// \short Virtual function that is used to symmetrise the problem so that
  /// the current solution exactly satisfies any symmetries within the system.
  /// Used when adpativly solving pitchfork detection problems when small
  /// asymmetries in the coarse solution can be magnified
  /// leading to very inaccurate answers on the fine mesh.
  /// This is always problem-specific and must be filled in by the user
  /// The default issues a warning
  void Problem::symmetrise_eigenfunction_for_adaptive_pitchfork_tracking()
  {
    std::ostringstream warn_message;
    warn_message
      << "Warning: This function is called after spatially adapting the\n"
      << "eigenfunction associated with a pitchfork bifurcation and should\n"
      << "ensure that the exact (anti-)symmetries of problem are enforced\n"
      << "within that eigenfunction. It is problem specific and must be\n"
      << "filled in by the user if required.\n"
      << "A sign of problems is if the slack paramter gets too large and\n"
      << "if the solution at the Pitchfork is not symmetric.\n";
    OomphLibWarning(
      warn_message.str(),
      "Problem::symmetrise_eigenfunction_for_adaptive_pitchfork_tracking()",
      OOMPH_EXCEPTION_LOCATION);
  }

  //====================================================================
  /// Return pointer to the parameter that is used in the
  /// bifurcation detection. If we are not tracking a bifurcation then
  /// an error will be thrown by the AssemblyHandler
  //====================================================================
  double* Problem::bifurcation_parameter_pt() const
  {
    return Assembly_handler_pt->bifurcation_parameter_pt();
  }

  //====================================================================
  /// Return the eigenfunction calculated as part of a
  /// bifurcation tracking process. If we are not tracking a bifurcation
  /// then an error will be thrown by the AssemblyHandler
  //======================================================================
  void Problem::get_bifurcation_eigenfunction(
    Vector<DoubleVector>& eigenfunction)
  {
    // Simply call the appropriate assembly handler function
    Assembly_handler_pt->get_eigenfunction(eigenfunction);
  }

  //============================================================
  /// Activate the fold tracking system by changing the assembly
  /// handler and initialising it using the parameter addressed
  /// by parameter_pt.
  //============================================================
  void Problem::activate_fold_tracking(double* const& parameter_pt,
                                       const bool& block_solve)
  {
    // Reset the assembly handler to default
    reset_assembly_handler_to_default();
    // Set the new assembly handler. Note that the constructor actually
    // solves the original problem to get some initial conditions, but
    // this is OK because the RHS is always evaluated before assignment.
    Assembly_handler_pt = new FoldHandler(this, parameter_pt);

    // If we are using a block solver, we must set the linear solver pointer
    // to the block fold solver. The present linear solver is
    // used by the block solver and so must be passed as an argument.
    // The destructor of the Fold handler returns the linear
    // solver to the original non-block version.
    if (block_solve)
    {
      Linear_solver_pt = new AugmentedBlockFoldLinearSolver(Linear_solver_pt);
    }
  }

  //===============================================================
  /// Activate the generic bifurcation ///tracking system by changing the
  /// assembly handler and initialising it using the parameter addressed by
  /// parameter_pt.
  //============================================================
  void Problem::activate_bifurcation_tracking(double* const& parameter_pt,
                                              const DoubleVector& eigenvector,
                                              const bool& block_solve)
  {
    // Reset the assembly handler to default
    reset_assembly_handler_to_default();
    // Set the new assembly handler. Note that the constructor actually
    // solves the original problem to get some initial conditions, but
    // this is OK because the RHS is always evaluated before assignment.
    Assembly_handler_pt = new FoldHandler(this, parameter_pt, eigenvector);

    // If we are using a block solver, we must set the linear solver pointer
    // to the block fold solver. The present linear solver is
    // used by the block solver and so must be passed as an argument.
    // The destructor of the Fold handler returns the linear
    // solver to the original non-block version.
    if (block_solve)
    {
      Linear_solver_pt = new AugmentedBlockFoldLinearSolver(Linear_solver_pt);
    }
  }


  //===============================================================
  /// Activate the generic bifurcation ///tracking system by changing the
  /// assembly handler and initialising it using the parameter addressed by
  /// parameter_pt.
  //============================================================
  void Problem::activate_bifurcation_tracking(double* const& parameter_pt,
                                              const DoubleVector& eigenvector,
                                              const DoubleVector& normalisation,
                                              const bool& block_solve)
  {
    // Reset the assembly handler to default
    reset_assembly_handler_to_default();
    // Set the new assembly handler. Note that the constructor actually
    // solves the original problem to get some initial conditions, but
    // this is OK because the RHS is always evaluated before assignment.
    Assembly_handler_pt =
      new FoldHandler(this, parameter_pt, eigenvector, normalisation);

    // If we are using a block solver, we must set the linear solver pointer
    // to the block fold solver. The present linear solver is
    // used by the block solver and so must be passed as an argument.
    // The destructor of the Fold handler returns the linear
    // solver to the original non-block version.
    if (block_solve)
    {
      Linear_solver_pt = new AugmentedBlockFoldLinearSolver(Linear_solver_pt);
    }
  }


  //==================================================================
  /// Activate the pitchfork tracking system by changing the assembly
  /// handler and initialising it using the parameter addressed
  /// by parameter_pt and a symmetry vector. The boolean flag is
  /// used to specify whether a block solver is used, default is true.
  //===================================================================
  void Problem::activate_pitchfork_tracking(double* const& parameter_pt,
                                            const DoubleVector& symmetry_vector,
                                            const bool& block_solve)
  {
    // Reset the assembly handler to default
    reset_assembly_handler_to_default();

    // Set the new assembly handler. Note that the constructor actually
    // solves the original problem to get some initial conditions, but
    // this is OK because the RHS is always evaluated before assignment.
    Assembly_handler_pt = new PitchForkHandler(
      this, this->assembly_handler_pt(), parameter_pt, symmetry_vector);

    // If we are using a block solver, we must set the linear solver pointer
    // to the block pitchfork solver. The present linear solver is
    // used by the block solver and so must be passed as an argument.
    // The destructor of the PitchFork handler returns the linear
    // solver to the original non-block version.
    if (block_solve)
    {
      Linear_solver_pt = new BlockPitchForkLinearSolver(Linear_solver_pt);
    }
  }


  //============================================================
  /// Activate the hopf tracking system by changing the assembly
  /// handler and initialising it using the parameter addressed
  /// by parameter_pt.
  //============================================================
  void Problem::activate_hopf_tracking(double* const& parameter_pt,
                                       const bool& block_solve)
  {
    // Reset the assembly handler to default
    reset_assembly_handler_to_default();
    // Set the new assembly handler. Note that the constructor actually
    // solves the original problem to get some initial conditions, but
    // this is OK because the RHS is always evaluated before assignment.
    Assembly_handler_pt = new HopfHandler(this, parameter_pt);

    // If we are using a block solver, we must set the linear solver pointer
    // to the block hopf solver. The present linear solver is
    // used by the block solver and so must be passed as an argument.
    // The destructor of the Hopf handler returns the linear
    // solver to the original non-block version.
    if (block_solve)
    {
      Linear_solver_pt = new BlockHopfLinearSolver(Linear_solver_pt);
    }
  }


  //============================================================
  /// Activate the hopf tracking system by changing the assembly
  /// handler and initialising it using the parameter addressed
  /// by parameter_pt and the frequency and null vectors
  /// specified.
  //============================================================
  void Problem::activate_hopf_tracking(double* const& parameter_pt,
                                       const double& omega,
                                       const DoubleVector& null_real,
                                       const DoubleVector& null_imag,
                                       const bool& block_solve)
  {
    // Reset the assembly handler to default
    reset_assembly_handler_to_default();
    // Set the new assembly handler. Note that the constructor actually
    // solves the original problem to get some initial conditions, but
    // this is OK because the RHS is always evaluated before assignment.
    Assembly_handler_pt =
      new HopfHandler(this, parameter_pt, omega, null_real, null_imag);

    // If we are using a block solver, we must set the linear solver pointer
    // to the block hopf solver. The present linear solver is
    // used by the block solver and so must be passed as an argument.
    // The destructor of the Hopf handler returns the linear
    // solver to the original non-block version.
    if (block_solve)
    {
      Linear_solver_pt = new BlockHopfLinearSolver(Linear_solver_pt);
    }
  }


  //===============================================================
  /// Reset the assembly handler to default
  //===============================================================
  void Problem::reset_assembly_handler_to_default()
  {
    // If we have a non-default handler
    if (Assembly_handler_pt != Default_assembly_handler_pt)
    {
      // Delete the current assembly handler
      delete Assembly_handler_pt;
      // Reset the assembly handler
      Assembly_handler_pt = Default_assembly_handler_pt;
    }
  }

  //===================================================================
  /// This function takes one step of length ds in pseudo-arclength.The
  /// argument parameter_pt is a pointer to the parameter (global variable)
  /// that is being traded for arc-length. The function returns the next desired
  /// arc-length according to criteria based upon the desired number of Newton
  /// Iterations per solve.
  //=====================================================================
  double Problem::arc_length_step_solve(double* const& parameter_pt,
                                        const double& ds,
                                        const unsigned& max_adapt)
  {
    // First check that we shouldn't use the other interface
    // by checking that the parameter isn't already stored as data
    if (does_pointer_correspond_to_problem_data(parameter_pt))
    {
      std::ostringstream error_message;
      error_message
        << "The parameter addressed by " << parameter_pt << " with the value "
        << *parameter_pt
        << "\n is supposed to be used for arc-length contiunation,\n"
        << " but it is stored in a Data object used by the problem.\n\n"
        << "This is bad for two reasons:\n"
        << "1. If it's a variable in the problem, it must already have an\n"
           "associated equation, so it can't be used for continuation;\n"
        << "2. The problem data will be reorganised in memory during "
           "continuation,\n"
        << "   which means that the pointer will become invalid.\n\n"
        << "If you are sure that this is what you want to do you must:\n"
        << "A. Ensure that the value is pinned (don't worry we'll shout again "
           "if not)\n"
        << "B. Use the alternative interface\n"
        << "   Problem::arc_length_step_solve(Data*,unsigned,...)\n"
        << "   which uses a pointer to the data object and not the raw double "
           "pointer."
        << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


    // If we are using the continuation timestepper
    if (Use_continuation_timestepper)
    {
      // Has the timestepper already been added to the problem
      bool continuation_time_stepper_added = false;
      const unsigned n_time_steppers = this->ntime_stepper();
      for (unsigned i = 0; i < n_time_steppers; i++)
      {
        if (this->time_stepper_pt(i) == &Continuation_time_stepper)
        {
          continuation_time_stepper_added = true;
          break;
        }
      }

      // If not add it
      if (!continuation_time_stepper_added)
      {
        oomph_info << "Adding the continuation time stepper\n";
        this->add_time_stepper_pt(&Continuation_time_stepper);
      }

      // Need to treat case of eigenproblems and bifurcation detection/tracking
      // here

      // Backup the current timesteppers for each mesh!


      // If an arc length step has not been taken then set the timestepper
      if (!Arc_length_step_taken)
      {
        // Set the continuation timestepper for all data in the problem
        oomph_info << this->set_timestepper_for_all_data(
                        &Continuation_time_stepper)
                   << " equation numbers allocated for continuation\n";
      }

    } // End of continuation time stepper case


    // Just call the helper function (parameter is not from data)
    return arc_length_step_solve_helper(parameter_pt, ds, max_adapt);
  }


  //===================================================================
  /// This function takes one step of length ds in pseudo-arclength.The
  /// argument data_pt is a pointer to the data that holds the
  /// parameter (global variable)
  /// that is being traded for arc-length. The exact value is located at
  /// the location given by data_index.
  /// The function returns the next desired
  /// arc-length according to criteria based upon the desired number of Newton
  /// Iterations per solve.
  //=====================================================================
  double Problem::arc_length_step_solve(Data* const& data_pt,
                                        const unsigned& data_index,
                                        const double& ds,
                                        const unsigned& max_adapt)
  {
    // Firstly check that the data is pinned
    if (!data_pt->is_pinned(data_index))
    {
      std::ostringstream error_stream;
      error_stream << "The value at index " << data_index
                   << " in the data object to be used for continuation\n"
                   << "is not pinned, which means that it is already a\n"
                   << "variable in the problem "
                   << "and cannot be used for continuation.\n\n"
                   << "Please correct your formulation by either:\n"
                   << "A. Pinning the value"
                   << "\n or \n"
                   << "B. Using a different parameter for continuation"
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


    // If we are using the continuation timestepper
    if (Use_continuation_timestepper)
    {
      // Has the timestepper already been added to the problem
      bool continuation_time_stepper_added = false;
      const unsigned n_time_steppers = this->ntime_stepper();
      for (unsigned i = 0; i < n_time_steppers; i++)
      {
        if (this->time_stepper_pt(i) == &Continuation_time_stepper)
        {
          continuation_time_stepper_added = true;
          break;
        }
      }

      // If not add it
      if (!continuation_time_stepper_added)
      {
        oomph_info << "Adding the continuation time stepper\n";
        this->add_time_stepper_pt(&Continuation_time_stepper);
      }

      // Need to treat case of eigenproblems and bifurcation detection/tracking
      // here


      // Backup the current timesteppers for each mesh!


      // If an arc length step has not been taken then set the timestepper
      if (!Arc_length_step_taken)
      {
        // Set the continuation timestepper for all data in the problem
        oomph_info << this->set_timestepper_for_all_data(
                        &Continuation_time_stepper)
                   << " equation numbers allocated for continuation\n";
      }


    } // End of continuation time stepper case


    // Now make a pointer to the (newly allocated) data object
    double* parameter_pt = data_pt->value_pt(data_index);
    // Call the helper function, this will change the parameter_pt if
    // the data storage is changed (if the timestepper has to be changed,
    // which happens if this is the first time that a continuation step is
    // taken)
    // ALH: Don't think this is true because it has happened above....
    return arc_length_step_solve_helper(parameter_pt, ds, max_adapt);
  }

  //======================================================================
  /// \short Private helper function that is used to set the appropriate
  /// pinned values for continuation. If the data is pinned, the its
  /// current value is always the same as the original value and
  /// the derivative is always zero. If these are not set properly
  /// then interpolation and projection in spatial adaptivity will
  /// not give the best answers.
  //=====================================================================
  void Problem::set_consistent_pinned_values_for_continuation()
  {
    // Set the consistent values for the global mesh
    Mesh_pt->set_consistent_pinned_values_for_continuation(
      &Continuation_time_stepper);

    // Deal with the spine meshes additional numbering separately
    const unsigned n_sub_mesh = this->nsub_mesh();
    // If there is only one mesh
    if (n_sub_mesh == 0)
    {
      if (SpineMesh* const spine_mesh_pt = dynamic_cast<SpineMesh*>(Mesh_pt))
      {
        spine_mesh_pt->set_consistent_pinned_spine_values_for_continuation(
          &Continuation_time_stepper);
      }
      // If it's a triangle mesh the we need to set the
    }
    // Otherwise loop over the sub meshes
    else
    {
      // Assign global equation numbers first
      for (unsigned i = 0; i < n_sub_mesh; i++)
      {
        if (SpineMesh* const spine_mesh_pt =
              dynamic_cast<SpineMesh*>(Sub_mesh_pt[i]))
        {
          spine_mesh_pt->set_consistent_pinned_spine_values_for_continuation(
            &Continuation_time_stepper);
        }
      }
    }

    // Also set time stepper for global data
    const unsigned n_global = Global_data_pt.size();
    for (unsigned i = 0; i < n_global; ++i)
    {
      Continuation_time_stepper.set_consistent_pinned_values(Global_data_pt[i]);
    }
  }


  //===================================================================
  /// This function takes one step of length ds in pseudo-arclength.The
  /// argument parameter_pt is a pointer to the parameter (global variable)
  /// that is being traded for arc-length. The function returns the next desired
  /// arc-length according to criteria based upon the desired number of Newton
  /// Iterations per solve.
  //=====================================================================
  double Problem::arc_length_step_solve_helper(double* const& parameter_pt,
                                               const double& ds,
                                               const unsigned& max_adapt)
  {
    //----------------------MAKE THE PROBLEM STEADY-----------------------
    // Loop over the timesteppers and make them (temporarily) steady.
    // We can only do continuation for steady problems!
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


    // Max number of solves
    unsigned max_solve = max_adapt + 1;
    // Storage for newton steps in each adaptation
    unsigned max_count_in_adapt_loop = 0;


    //----SET UP MEMORY FOR QUANTITIES THAT ARE REQUIRED OUTSIDE THE LOOP----

    // Assign memory for solutions of the equations Jz = du/dparameter
    // This is needed here (outside the loop), so that we can save on
    // one linear solve when calculating the derivatives wrt the arc-length
    DoubleVector z;


    // Store sign of the Jacobian, used for bifurcation detection
    // If this is the first time that we are calling the arc-length solver,
    // this should not be used.
    int previous_sign = Sign_of_jacobian;

    // Flag to indicate a sign change
    bool SIGN_CHANGE = false;


    // Adaptation loop
    for (unsigned isolve = 0; isolve < max_solve; ++isolve)
    {
      // Only adapt after the first solve has been done
      if (isolve > 0)
      {
        unsigned n_refined;
        unsigned n_unrefined;

        // Adapt problem
        adapt(n_refined, n_unrefined);

#ifdef OOMPH_HAS_MPI
        // Adaptation only converges if ALL the processes have no
        // refinement or unrefinement to perform
        unsigned total_refined = 0;
        unsigned total_unrefined = 0;
        if (Problem_has_been_distributed)
        {
          MPI_Allreduce(&n_refined,
                        &total_refined,
                        1,
                        MPI_UNSIGNED,
                        MPI_SUM,
                        this->communicator_pt()->mpi_comm());
          n_refined = total_refined;
          MPI_Allreduce(&n_unrefined,
                        &total_unrefined,
                        1,
                        MPI_UNSIGNED,
                        MPI_SUM,
                        this->communicator_pt()->mpi_comm());
          n_unrefined = total_unrefined;
        }
#endif

        oomph_info << "---> " << n_refined << " elements were refined, and "
                   << n_unrefined << " were unrefined"
#ifdef OOMPH_HAS_MPI
                   << ", in total (over all processors).\n";
#else
                   << ".\n";
#endif


        // Check convergence of adaptation cycle
        if ((n_refined == 0) && (n_unrefined == 0))
        {
          oomph_info << "\n \n Solution is fully converged in "
                     << "Problem::newton_solver(). \n \n ";
          break;
        }
      }

      //----------SAVE THE INITIAL VALUES, IN CASE THE STEP FAILS-----------

      // Find the number of local dofs
      unsigned ndof_local = Dof_distribution_pt->nrow_local();

      // Only need to do this in the first loop
      if (isolve == 0)
      {
        if (!Use_continuation_timestepper)
        {
          // Safety check, set up the array of dof derivatives, if necessary
          // The distribution is the same as the (natural) distribution of the
          // dofs
          if (Dof_derivative.size() != ndof_local)
          {
            Dof_derivative.resize(ndof_local, 0.0);
          }

          // Safety check, set up the array of curren values, if necessary
          // Again the distribution reflects the (natural) distribution of the
          // dofs
          if (Dof_current.size() != ndof_local)
          {
            Dof_current.resize(ndof_local);
          }
        }

        // Save the current value of the parameter
        Parameter_current = *parameter_pt;

        // Save the current values of the degrees of freedom
        for (unsigned long l = 0; l < ndof_local; l++)
        {
          dof_current(l) = *Dof_pt[l];
        }

        // Set the value of ds_current
        Ds_current = ds;
      }

      // Counter for the number of newton steps
      unsigned count = 0;

      // Flag to indicate a successful step
      bool STEP_REJECTED = false;


      // Set the appropriate initial conditions for the pinned data
      if (Use_continuation_timestepper)
      {
        this->set_consistent_pinned_values_for_continuation();
      }

      // Loop around the step in arc-length
      do
      {
        // Check that the step has not fallen below the minimum tolerance
        if (std::fabs(Ds_current) < Minimum_ds)
        {
          std::ostringstream error_message;
          error_message << "DESIRED ARC-LENGTH STEP " << Ds_current
                        << " HAS FALLEN BELOW MINIMUM TOLERANCE, " << Minimum_ds
                        << std::endl;

          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        // Assume that we shall accept the step
        STEP_REJECTED = false;

        // Set initial value of the parameter
        *parameter_pt = Parameter_current + Parameter_derivative * Ds_current;

        // Perform any actions...
        actions_after_parameter_increase(parameter_pt);

        Ds_current = (*parameter_pt - Parameter_current) / Parameter_derivative;

        // Loop over the (local) variables and set their initial values
        for (unsigned long l = 0; l < ndof_local; l++)
        {
          *Dof_pt[l] = dof_current(l) + dof_derivative(l) * Ds_current;
        }

        // Actually do the newton solve stage for the continuation problem
        try
        {
          count = newton_solve_continuation(parameter_pt, z);
        }
        // Catch any exceptions thrown in the Newton solver
        catch (NewtonSolverError& error)
        {
          // Check whether it's the linear solver
          if (error.linear_solver_error)
          {
            std::ostringstream error_stream;
            error_stream << std::endl
                         << "USER-DEFINED ERROR IN NEWTON SOLVER " << std::endl;
            oomph_info << "ERROR IN THE LINEAR SOLVER" << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
          // Otherwise mark the step as having failed
          else
          {
            oomph_info << "STEP REJECTED --- TRYING AGAIN" << std::endl;
            STEP_REJECTED = true;
            // Let's take a smaller step
            Ds_current *= (2.0 / 3.0);
          }
        }
      } while (STEP_REJECTED); // continue until a step is accepted

      // Set the maximum count
      if (count > max_count_in_adapt_loop)
      {
        max_count_in_adapt_loop = count;
      }
    } /// end of adaptation loop

    // Only recalculate the derivatives if there has been a Newton solve
    // If not, the previous values should be close enough
    if (max_count_in_adapt_loop > 0)
    {
      //--------------------CHECK FOR POTENTIAL BIFURCATIONS-------------
      if (Bifurcation_detection)
      {
        // If the sign of the jacobian is zero issue a warning
        if (Sign_of_jacobian == 0)
        {
          std::string error_message =
            "The sign of the jacobian is zero after a linear solve\n";
          error_message += "Either the matrix is singular (unlikely),\n";
          error_message += "or the linear solver cannot compute the "
                           "determinant of the matrix;\n";
          error_message += "e.g. an iterative linear solver.\n";
          error_message +=
            "If the latter, bifurcation detection must be via an eigensolver\n";
          OomphLibWarning(error_message,
                          "Problem::arc_length_step_solve",
                          OOMPH_EXCEPTION_LOCATION);
        }
        // If this is the first step, we cannot rely on the previous value
        // of the jacobian so set the previous sign to the present sign
        if (!Arc_length_step_taken)
        {
          previous_sign = Sign_of_jacobian;
        }
        // If we have detected a sign change in the last converged Jacobian,
        // it must be a turning point or bifurcation
        if (Sign_of_jacobian != previous_sign)
        {
          // There has been, at least, one sign change
          First_jacobian_sign_change = true;

          // The sign has changed this time
          SIGN_CHANGE = true;

          // Calculate the dot product of the approximate null vector
          // of the Jacobian matrix ((badly) approximated by z)
          // and the vectors of derivatives of the residuals wrt the
          // global parameter
          // If this is small it is a bifurcation rather than a turning point.
          // Get the derivative wrt global parameter
          // DoubleVector dparam;
          // get_derivative_wrt_global_parameter(parameter_pt,dparam);
          // Calculate the dot product
          // double dot=0.0;
          // for(unsigned long n=0;n<n_dofs;++n) {dot += dparam[n]*z[n];}
          // z.dot(dparam);

          // Write the output message
          std::ostringstream message;
          message
            << "-----------------------------------------------------------";
          message << std::endl
                  << "SIGN CHANGE IN DETERMINANT OF JACOBIAN: " << std::endl;
          message << "BIFURCATION OR TURNING POINT DETECTED BETWEEN "
                  << Parameter_current << " AND " << *parameter_pt << std::endl;
          // message << "APPROXIMATE DOT PRODUCT : " << dot << "," << std::endl;
          // message << "IF CLOSE TO ZERO WE HAVE A BIFURCATION; ";
          // message << "OTHERWISE A TURNING POINT" << std::endl;
          message
            << "-----------------------------------------------------------"
            << std::endl;

          // Write the message to standard output
          oomph_info << message.str();

          // Open the information file for appending
          std::ofstream bifurcation_info("bifurcation_info",
                                         std::ios_base::app);
          // Write the message to the file
          bifurcation_info << message.str();
          bifurcation_info.close();
        }
      }

      // Calculate the derivatives required for the next stage of continuation
      // In this we pass the last value of z (i.e. approximation)
      if (!Use_finite_differences_for_continuation_derivatives)
      {
        calculate_continuation_derivatives(z);
      }
      // Or use finite differences
      else
      {
        calculate_continuation_derivatives_fd(parameter_pt);
      }

      // If it's the first step then the value of the next step should
      // be the change in parameter divided by the parameter derivative
      // to obtain approximately the same parameter change
      if (!Arc_length_step_taken)
      {
        Ds_current = (*parameter_pt - Parameter_current) / Parameter_derivative;
      }

      // We have taken our first step
      Arc_length_step_taken = true;
    }
    // If there has not been a newton step then we still need to estimate
    // the derivatives in the arc length direction
    else
    {
      // Default is to calculate the continuation derivatives by solving the
      // linear system. We must do this to ensure that the derivatives are in
      // sync It could lead to problems near turning points when we should
      // really be solving an eigenproblem, but seems OK so far!

      // Save the current sign of the jacobian
      int temp_sign = Sign_of_jacobian;

      // Calculate the continuation derivatives, which includes a solve
      // of the linear system if not using finite differences
      if (!Use_finite_differences_for_continuation_derivatives)
      {
        calculate_continuation_derivatives(parameter_pt);
      }
      // Otherwise use finite differences
      else
      {
        calculate_continuation_derivatives_fd(parameter_pt);
      }

      // Reset the sign of the jacobian, just in case the sign has changed when
      // solving the continuation derivatives. The sign change will be picked
      // up on the next continuation step.
      Sign_of_jacobian = temp_sign;
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

    // If we are trying to find a bifurcation and the first sign change
    // has occured, use bisection
    if ((Bifurcation_detection) && (Bisect_to_find_bifurcation) &&
        (First_jacobian_sign_change))
    {
      // If there has been a sign change we need to half the step size
      // and reverse the direction
      if (SIGN_CHANGE)
      {
        Ds_current *= -0.5;
      }
      // Otherwise
      else
      {
        // The size of the bracketed interval is always
        // 2ds - Ds_current (this will work even if the original step failed)
        // We want our new step size to be half this
        Ds_current = ds - 0.5 * Ds_current;
      }
      // Return the desired value of the step
      return Ds_current;
    }

    // If fewer than the desired number of Newton Iterations, increase the step
    if (max_count_in_adapt_loop < Desired_newton_iterations_ds)
    {
      return Ds_current * 1.5;
    }
    // If more than the desired number of Newton Iterations, reduce the step
    if (max_count_in_adapt_loop > Desired_newton_iterations_ds)
    {
      return Ds_current * (2.0 / 3.0);
    }
    // Otherwise return the step just taken
    return Ds_current;
  }


  //=======================================================================
  /// Take an explicit timestep of size dt
  //======================================================================
  void Problem::explicit_timestep(const double& dt, const bool& shift_values)
  {
#ifdef PARANOID
    if (this->explicit_time_stepper_pt() == 0)
    {
      throw OomphLibError("Explicit time stepper pointer is null in problem.",
                          OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
    }
#endif

    // Firstly we shift the time values
    if (shift_values)
    {
      shift_time_values();
    }
    // Set the current value of dt, if we can
    if (time_pt()->ndt() > 0)
    {
      time_pt()->dt() = dt;
    }

    // Take the explicit step
    this->explicit_time_stepper_pt()->timestep(this, dt);
  }


  //========================================================================
  /// Do one timestep of size dt using Newton's method with the specified
  /// tolerance and linear solver defined as member data of the Problem class.
  /// This will be the most commonly used version
  /// of  unsteady_newton_solve, in which the time values are always shifted
  /// This does not include any kind of adaptativity. If the solution fails to
  /// converge the program will end.
  //========================================================================
  void Problem::unsteady_newton_solve(const double& dt)
  {
    // We shift the values, so shift_values is true
    unsteady_newton_solve(dt, true);
  }

  //========================================================================
  /// Do one timestep forward of size dt using Newton's method with the
  /// specified tolerance and linear solver defined via member data of the
  /// Problem class.
  /// The boolean flag shift_values is used to control whether the time values
  /// should be shifted or not.
  //========================================================================
  void Problem::unsteady_newton_solve(const double& dt,
                                      const bool& shift_values)
  {
    // Shift the time values and the dts, according to the control flag
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

    // Run the individual timesteppers actions before timestep. These need to
    // be before the problem's actions_before_implicit_timestep so that the
    // boundary conditions are set consistently.
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      time_stepper_pt(i)->actions_before_timestep(this);
    }

    // Now update anything that needs updating before the timestep
    // This could be time-dependent boundary conditions, for example.
    actions_before_implicit_timestep();

    try
    {
      // Solve the non-linear problem for this timestep with Newton's method
      newton_solve();
    }
    // Catch any exceptions thrown in the Newton solver
    catch (NewtonSolverError& error)
    {
      oomph_info << std::endl
                 << "USER-DEFINED ERROR IN NEWTON SOLVER " << std::endl;
      // Check whether it's the linear solver
      if (error.linear_solver_error)
      {
        oomph_info << "ERROR IN THE LINEAR SOLVER" << std::endl;
      }
      // Check to see whether we have reached Max_iterations
      else if (error.iterations == Max_newton_iterations)
      {
        oomph_info << "MAXIMUM NUMBER OF ITERATIONS (" << error.iterations
                   << ") REACHED WITHOUT CONVERGENCE " << std::endl;
      }
      // If not, it must be that we have exceeded the maximum residuals
      else
      {
        oomph_info << "MAXIMUM RESIDUALS: " << error.maxres
                   << " EXCEEDS PREDEFINED MAXIMUM " << Max_residuals
                   << std::endl;
      }
      // Die horribly!!
      std::ostringstream error_stream;
      error_stream << "Error occured in unsteady Newton solver. " << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Run the individual timesteppers actions, these need to be before the
    // problem's actions_after_implicit_timestep so that the time step is
    // finished before the problem does any auxiliary calculations (e.g. in
    // semi-implicit micromagnetics the calculation of magnetostatic field).
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      time_stepper_pt(i)->actions_after_timestep(this);
    }


    // Now update anything that needs updating after the timestep
    actions_after_implicit_timestep();
    actions_after_implicit_timestep_and_error_estimation();
  }

  //=======================================================================
  /// Attempt to take one timestep forward using dt_desired. The error control
  /// parameter, epsilon, is used to specify the desired approximate value of
  /// the global error norm per timestep. The routine returns the value an
  /// estimate of the next value of dt that should be taken.
  //=======================================================================
  double Problem::adaptive_unsteady_newton_solve(const double& dt_desired,
                                                 const double& epsilon)
  {
    // We always want to shift the time values
    return adaptive_unsteady_newton_solve(dt_desired, epsilon, true);
  }


  //=======================================================================
  /// Attempt to take  one timestep forward using the dt_desired.
  /// This is the driver for a number of adaptive solvers. If the solution
  /// fails to converge at a given timestep, the routine will automatically
  /// halve the time step and try again, until the time step falls below the
  /// specified minimum value. The routine returns the value an estimate
  /// of the next value of dt that should be taken.
  /// Timestep is also rejected if the  error estimate post-solve
  /// (computed by global_temporal_error_norm()) exceeds epsilon.
  /// This behaviour can be over-ruled by setting the protected
  /// boolean Problem::Keep_temporal_error_below_tolerance to false.
  //========================================================================
  double Problem::adaptive_unsteady_newton_solve(const double& dt_desired,
                                                 const double& epsilon,
                                                 const bool& shift_values)
  {
    // First, we need to backup the existing dofs, in case the timestep is
    // rejected

    // Find total number of dofs on current processor
    unsigned n_dof_local = dof_distribution_pt()->nrow_local();

    // Now set up a Vector to hold current values
    Vector<double> dofs_current(n_dof_local);

    // Load values into dofs_current
    for (unsigned i = 0; i < n_dof_local; i++) dofs_current[i] = dof(i);

    // Store the time
    double time_current = time_pt()->time();

    // Flag to detect whether the timestep has been rejected or not
    bool reject_timestep = 0;

    // Flag to detect whether any of the timesteppers are adaptive
    unsigned adaptive_flag = 0;

    // The value of the actual timestep, by default the same as desired timestep
    double dt_actual = dt_desired;

    // Find out whether any of the timesteppers are adaptive
    unsigned n_time_steppers = ntime_stepper();
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      if (time_stepper_pt(i)->adaptive_flag())
      {
        adaptive_flag = 1;
        break;
      }
    }

    // Shift the time_values according to the control flag
    if (shift_values)
    {
      shift_time_values();
    }

    // This loop surrounds the adaptive time-stepping and will not be broken
    // until a timestep is accepted
    do
    {
      // Initially we assume that this step will succeed and that this dt
      // value is ok.
      reject_timestep = 0;
      double dt_rescaling_factor = 1.0;

      // Set the new time and value of dt
      time_pt()->time() += dt_actual;
      time_pt()->dt() = dt_actual;

      // Loop over all timesteppers and set the weights and predictor weights
      for (unsigned i = 0; i < n_time_steppers; i++)
      {
        // If the time_stepper is non-adaptive, this will be zero
        time_stepper_pt(i)->set_predictor_weights();
        time_stepper_pt(i)->set_weights();
      }

      // Now calculate the predicted values for the all data and all positions
      calculate_predictions();

      // Run the individual timesteppers actions before timestep. These need to
      // be before the problem's actions_before_implicit_timestep so that the
      // boundary conditions are set consistently.
      for (unsigned i = 0; i < n_time_steppers; i++)
      {
        time_stepper_pt(i)->actions_before_timestep(this);
      }

      // Do any updates/boundary conditions changes here
      actions_before_implicit_timestep();

      // Attempt to solve the non-linear system
      try
      {
        // Solve the non-linear problem at this timestep
        newton_solve();
      }
      // Catch any exceptions thrown
      catch (NewtonSolverError& error)
      {
        // If it's a solver error then die
        if (error.linear_solver_error ||
            Time_adaptive_newton_crash_on_solve_fail)
        {
          std::string error_message = "USER-DEFINED ERROR IN NEWTON SOLVER\n";
          error_message += "ERROR IN THE LINEAR SOLVER\n";

          // Die
          throw OomphLibError(
            error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
        else
        {
          // Reject the timestep, if we have an exception
          oomph_info << "TIMESTEP REJECTED" << std::endl;
          reject_timestep = 1;

          // Half the time step
          dt_rescaling_factor = Timestep_reduction_factor_after_nonconvergence;
        }
      }

      // Run the individual timesteppers actions, these need to be before the
      // problem's actions_after_implicit_timestep so that the time step is
      // finished before the problem does any auxiliary calculations (e.g. in
      // semi-implicit micromagnetics the calculation of magnetostatic field).
      for (unsigned i = 0; i < n_time_steppers; i++)
      {
        time_stepper_pt(i)->actions_after_timestep(this);
      }

      // Update anything that needs updating after the timestep
      actions_after_implicit_timestep();

      // If we have an adapative timestepper (and we haven't already failed)
      // then calculate the error estimate and rescaling factor.
      if (adaptive_flag && !reject_timestep)
      {
        // Once timestep has been accepted can do fancy error processing
        // Set the error weights
        for (unsigned i = 0; i < n_time_steppers; i++)
        {
          time_stepper_pt(i)->set_error_weights();
        }

        // Get a global error norm to use in adaptivity (as specified by the
        // problem sub-class writer). Prevent a divide by zero if the solution
        // gives very close to zero error. Error norm should never be negative
        // but use absolute value just in case.
        double error = std::max(std::abs(global_temporal_error_norm()), 1e-12);

        // Calculate the scaling  factor
        dt_rescaling_factor = std::pow(
          (epsilon / error), (1.0 / (1.0 + time_stepper_pt()->order())));

        oomph_info << "Timestep scaling factor is  " << dt_rescaling_factor
                   << std::endl;
        oomph_info << "Estimated timestepping error is " << error << std::endl;


        // Do we have to do it again?
        if (error > epsilon)
        {
          oomph_info << "Estimated timestepping error " << error
                     << " exceeds tolerance " << epsilon << " ";
          if (Keep_temporal_error_below_tolerance)
          {
            oomph_info << " --> rejecting timestep." << std::endl;
            reject_timestep = 1;
          }
          else
          {
            oomph_info << " ...but we're not rejecting the timestep"
                       << std::endl;
          }
          oomph_info
            << "Note: This behaviour can be adjusted by changing the protected "
            << "boolean" << std::endl
            << std::endl
            << "    Problem::Keep_temporal_error_below_tolerance" << std::endl;
        }


      } // End of if adaptive flag


      // Calculate the next time step size and check it's ok
      // ============================================================

      // Calculate the possible next time step, if no error conditions
      // trigger.
      double new_dt_candidate = dt_rescaling_factor * dt_actual;

      // Check that the scaling factor is within the allowed range
      if (dt_rescaling_factor > DTSF_max_increase)
      {
        oomph_info << "Tried to increase dt by the ratio "
                   << dt_rescaling_factor << " which is above the maximum ("
                   << DTSF_max_increase
                   << "). Attempting to increase by the maximum ratio instead."
                   << std::endl;
        new_dt_candidate = DTSF_max_increase * dt_actual;
      }
      // If we have already rejected the timestep then don't do this check
      // because DTSF will definitely be too small.
      else if ((!reject_timestep) && (dt_rescaling_factor <= DTSF_min_decrease))
      {
        // Handle this special case where we want to continue anyway (usually
        // Minimum_dt_but_still_proceed = -1 so this has no effect).
        if (new_dt_candidate < Minimum_dt_but_still_proceed)
        {
          oomph_info
            << "Warning: Adaptation of timestep to ensure satisfaction\n"
            << "         of error bounds during adaptive timestepping\n"
            << "         would lower dt below \n"
            << "         Problem::Minimum_dt_but_still_proceed="
            << Minimum_dt_but_still_proceed << "\n"
            << "         ---> We're continuing with present timestep.\n"
            << std::endl;
          dt_rescaling_factor = 1.0;
          // ??ds shouldn't we set new_dt_candidate =
          // Minimum_dt_but_still_proceed here, rather than not changing dt at
          // all?
        }
        else
        {
          // Otherwise reject
          oomph_info << "Timestep would decrease by " << dt_rescaling_factor
                     << " which is less than the minimum scaling factor "
                     << DTSF_min_decrease << std::endl;
          oomph_info << "TIMESTEP REJECTED" << std::endl;
          reject_timestep = 1;
        }
      }

      // Now check that the new dt is within the allowed range
      if (new_dt_candidate > Maximum_dt)
      {
        oomph_info << "Tried to increase dt to " << new_dt_candidate
                   << " which is above the maximum (" << Maximum_dt
                   << "). I increased it to the maximum value instead.";
        dt_actual = Maximum_dt;
      }
      else if (new_dt_candidate < Minimum_dt)
      {
        std::ostringstream err;
        err << "Tried to reduce dt to " << new_dt_candidate
            << " which is less than the minimum dt (" << Minimum_dt << ")."
            << std::endl;
        throw OomphLibError(
          err.str(), OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
      }
      else
      {
        dt_actual = new_dt_candidate;
      }


      actions_after_implicit_timestep_and_error_estimation();


      // If we are rejecting this attempt then revert the dofs etc.
      if (reject_timestep)
      {
        // Reset the time
        time_pt()->time() = time_current;

        // Reload the dofs
        unsigned ni = dofs_current.size();
        for (unsigned i = 0; i < ni; i++)
        {
          dof(i) = dofs_current[i];
        }

#ifdef OOMPH_HAS_MPI
        // Synchronise the solution on different processors (on each submesh)
        this->synchronise_all_dofs();
#endif

        // Call all "after" actions, e.g. to handle mesh updates
        actions_after_newton_step();
        actions_before_newton_convergence_check();
        actions_after_newton_solve();
        actions_after_implicit_timestep();
        actions_after_implicit_timestep_and_error_estimation();
      }

    }
    // Keep this loop going until we accept the timestep
    while (reject_timestep);

    // Once the timestep has been accepted, return the time step that should be
    // used next time.
    return dt_actual;
  }


  //=======================================================================
  /// Private helper function to perform
  /// unsteady "doubly" adaptive Newton solve: Does temporal
  /// adaptation first, i.e. we try to do a timestep with an increment
  /// of dt, and adjusting dt until the solution on the given mesh satisfies
  /// the temporal error measure with tolerance epsilon. Following
  /// this, we do up to max_adapt spatial adaptions (without
  /// re-examining the temporal error). If first==true, the initial conditions
  /// are re-assigned after the mesh adaptations.
  /// Shifting of time can be suppressed by overwriting the
  /// default value of shift (true). [Shifting must be done
  /// if first_timestep==true because we're constantly re-assigning
  /// the initial conditions; if first_timestep==true and shift==false
  /// shifting is performed anyway and a warning is issued.
  /// Pseudo-Boolean flag suppress_resolve_after_spatial_adapt [0: false;
  /// 1: true] does what it says.]
  //========================================================================
  double Problem::doubly_adaptive_unsteady_newton_solve_helper(
    const double& dt_desired,
    const double& epsilon,
    const unsigned& max_adapt,
    const unsigned& suppress_resolve_after_spatial_adapt_flag,
    const bool& first,
    const bool& shift_values)
  {
    // Store the initial time
    double initial_time = time_pt()->time();

    // Take adaptive timestep, adjusting dt until tolerance is satisfied
    double new_dt =
      adaptive_unsteady_newton_solve(dt_desired, epsilon, shift_values);
    double dt_taken = time_pt()->dt();
    oomph_info << "Accepted solution taken with timestep: " << dt_taken
               << std::endl;


    // Bail out straightaway if no spatial adaptation allowed
    if (max_adapt == 0)
    {
      oomph_info << "No spatial refinement allowed; max_adapt=0\n";
      return new_dt;
    }

    // Adapt problem/mesh
    unsigned n_refined = 0;
    unsigned n_unrefined = 0;
    adapt(n_refined, n_unrefined);

    // Check if mesh has been adapted on other processors
    Vector<int> total_ref_count(2);
    total_ref_count[0] = n_refined;
    total_ref_count[1] = n_unrefined;


#ifdef OOMPH_HAS_MPI
    if (Problem_has_been_distributed)
    {
      // Sum n_refine across all processors
      Vector<int> ref_count(2);
      ref_count[0] = n_refined;
      ref_count[1] = n_unrefined;
      MPI_Allreduce(&ref_count[0],
                    &total_ref_count[0],
                    2,
                    MPI_INT,
                    MPI_SUM,
                    communicator_pt()->mpi_comm());
    }
#endif


    // Re-solve the problem if the adaptation has changed anything
    if ((total_ref_count[0] != 0) || (total_ref_count[1] != 0))
    {
      if (suppress_resolve_after_spatial_adapt_flag == 1)
      {
        oomph_info << "Mesh was adapted but re-solve has been suppressed."
                   << std::endl;
      }
      else
      {
        oomph_info
          << "Mesh was adapted --> we'll re-solve for current timestep."
          << std::endl;

        // Reset time to what it was when we entered here
        // because it will be incremented again by dt_taken.
        time_pt()->time() = initial_time;

        // Shift the timesteps? No! They've been shifted already when we
        // called the solve with pure temporal adaptivity...
        bool shift = false;

        // Reset the inital condition on refined meshes
        if (first)
        {
          // Reset default set_initial_condition has been called flag to false
          Default_set_initial_condition_called = false;

          // Reset the initial conditions
          oomph_info << "Re-assigning initial condition at time="
                     << time_pt()->time() << std::endl;
          set_initial_condition();

          // This is the first timestep so shifting
          // has to be done following the assignment of initial conditions,
          // providing the default set_initial_condition function has not
          // been called.
          // In fact, unsteady_newton_solve(...) does that automatically.
          // We're changing the flag here to avoid warning messages.
          if (!Default_set_initial_condition_called)
          {
            shift = true;
          }
        }

        // Now take the step again on the refined mesh, using the same
        // timestep as used before.
        unsteady_newton_solve(dt_taken, max_adapt, first, shift);
      }
    }
    else
    {
      oomph_info << "Mesh wasn't adapted --> we'll accept spatial refinement."
                 << std::endl;
    }

    return new_dt;
  }


  //========================================================================
  /// \short Initialise the previous values of the variables for time stepping
  /// corresponding to an impulsive start. Previous history for all data
  /// is generated by the appropriate timesteppers. Previous nodal
  /// positions are simply copied backwards.
  //========================================================================
  void Problem::assign_initial_values_impulsive()
  {
    // Assign the impulsive values in the "master" mesh
    Mesh_pt->assign_initial_values_impulsive();

    // Loop over global data
    unsigned Nglobal = Global_data_pt.size();
    for (unsigned iglobal = 0; iglobal < Nglobal; iglobal++)
    {
      Global_data_pt[iglobal]
        ->time_stepper_pt()
        ->assign_initial_values_impulsive(Global_data_pt[iglobal]);
    }
  }


  //=======================================================================
  /// Assign the values for an impulsive start and also set the initial
  /// values of the previous dts to both be dt
  //======================================================================
  void Problem::assign_initial_values_impulsive(const double& dt)
  {
    // First initialise the dts and set the weights
    initialise_dt(dt);
    // Now call assign_initial_values_impulsive
    assign_initial_values_impulsive();
  }

  //=======================================================================
  /// Return the current value of continuous time. If not Time object
  /// has been assigned, then throw an error
  //======================================================================
  double& Problem::time()
  {
    if (Time_pt == 0)
    {
      throw OomphLibError("Time object has not been set",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      return Time_pt->time();
    }
  }

  //=======================================================================
  /// Return the current value of continuous time. If not Time object
  /// has been assigned, then throw an error. Const version.
  //======================================================================
  double Problem::time() const
  {
    if (Time_pt == 0)
    {
      throw OomphLibError("Time object has not been set",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      return Time_pt->time();
    }
  }


  //=======================================================================
  /// Set all problem data to have the same timestepper (timestepper_pt).
  /// This is mainly used in continuation and bifurcation detection problems
  /// in which case the total number of unknowns may change and the changes
  /// to the underlying memory layout means that the Dof_pt must be
  /// reallocated. Thus, the function calls assign_eqn_numbers() and returns
  /// the number of new equation numbers.
  //=========================================================================
  unsigned long Problem::set_timestepper_for_all_data(
    TimeStepper* const& time_stepper_pt, const bool& preserve_existing_data)
  {
    // Set the timestepper for the master mesh's nodal and elemental data
    // to be the
    // continuation time stepper. This will wipe all storage other than
    // the 0th (present time) value at all the data objects
    Mesh_pt->set_nodal_and_elemental_time_stepper(time_stepper_pt,
                                                  preserve_existing_data);

    // Deal with the any additional mesh level timestepper data separately
    const unsigned n_sub_mesh = this->nsub_mesh();
    // If there is only one mesh
    if (n_sub_mesh == 0)
    {
      Mesh_pt->set_mesh_level_time_stepper(time_stepper_pt,
                                           preserve_existing_data);
    }
    // Otherwise loop over the sub meshes
    else
    {
      // Assign global equation numbers first
      for (unsigned i = 0; i < n_sub_mesh; i++)
      {
        this->Sub_mesh_pt[i]->set_mesh_level_time_stepper(
          time_stepper_pt, preserve_existing_data);
      }
    }

    // Also set time stepper for global data
    const unsigned n_global = Global_data_pt.size();
    for (unsigned i = 0; i < n_global; ++i)
    {
      Global_data_pt[i]->set_time_stepper(time_stepper_pt,
                                          preserve_existing_data);
    }

    // We now need to reassign equations numbers because the Dof pointer
    // will be inappropriate  because memory has been reallocated

#ifdef OOMPH_HAS_MPI
    if (Problem_has_been_distributed)
    {
      std::ostringstream warning_stream;
      warning_stream << "This has not been comprehensively tested for "
                        "distributed problems.\n"
                     << "I'm sure that I need to worry about external halo and "
                        "external elements."
                     << std::endl;
      OomphLibWarning(
        warning_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

#endif

    return (this->assign_eqn_numbers());
  }


  //========================================================================
  /// Shift all time-dependent data along for next timestep.
  //========================================================================
  void Problem::shift_time_values()
  {
    // Move the values of dt in the Time object
    Time_pt->shift_dt();

    // Only shift time values in the "master" mesh, otherwise things will
    // get shifted twice in complex problems
    Mesh_pt->shift_time_values();

    // Shift global data with their own timesteppers
    unsigned Nglobal = Global_data_pt.size();
    for (unsigned iglobal = 0; iglobal < Nglobal; iglobal++)
    {
      Global_data_pt[iglobal]->time_stepper_pt()->shift_time_values(
        Global_data_pt[iglobal]);
    }
  }


  //========================================================================
  /// Calculate the predictions of all variables in problem
  //========================================================================
  void Problem::calculate_predictions()
  {
// Check that if we have multiple time steppers none of them want to
// predict by calling an explicit timestepper (as opposed to doing
// something like an explicit step by combining known history values, as
// done in BDF).
#ifdef PARANOID
    if (Time_stepper_pt.size() != 1)
    {
      for (unsigned j = 0; j < Time_stepper_pt.size(); j++)
      {
        if (time_stepper_pt()->predict_by_explicit_step())
        {
          std::string err = "Prediction by explicit step only works for "
                            "problems with a simple time";
          err += "stepper. I think implementing anything more general will";
          err += "require a rewrite of explicit time steppers. - David";
          throw OomphLibError(
            err, OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
        }
      }
    }
#endif


    // Predict using an explicit timestepper (don't do it if adaptive = false
    // because pointers probably aren't set up).
    if (time_stepper_pt()->predict_by_explicit_step() &&
        time_stepper_pt()->adaptive_flag())
    {
      // Copy the time stepper's predictor pt into problem's explicit time
      // stepper pt (unless problem already has its own explicit time
      // stepper).
      ExplicitTimeStepper* ets_pt = time_stepper_pt()->explicit_predictor_pt();
#ifdef PARANOID
      if (ets_pt == 0)
      {
        std::string err = "Requested predictions by explicit step but explicit";
        err += " predictor pt is null.";
        throw OomphLibError(
          err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      if ((explicit_time_stepper_pt() != ets_pt) &&
          (explicit_time_stepper_pt() != 0))
      {
        throw OomphLibError("Problem has explicit time stepper other than "
                            "predictor, not sure how to handle this yet ??ds",
                            OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif
      explicit_time_stepper_pt() = ets_pt;

      // Backup dofs and time
      store_current_dof_values();

#ifdef PARANOID
      double backup_time = time();
#endif

      // Move time back so that we are at the start of the timestep (as
      // explicit_timestep functions expect). This is needed because the
      // predictor calculations are done after unsteady newton solve has
      // started, and it has already moved time forwards.
      double dt = time_pt()->dt();
      time() -= dt;

      // Explicit step
      this->explicit_timestep(dt, false);

      // Copy predicted dofs and time to their storage slots.
      set_dofs(time_stepper_pt()->predictor_storage_index(), Dof_pt);
      time_stepper_pt()->update_predicted_time(time());

      // Check we got the times right
#ifdef PARANOID
      if (std::abs(time() - backup_time) > 1e-12)
      {
        using namespace StringConversion;
        std::string err = "Predictor landed at the wrong time!";
        err += " Expected time " + to_string(backup_time, 14) + " but got ";
        err += to_string(time(), 14);
        throw OomphLibError(
          err, OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
      }
#endif

      // Restore dofs and time
      restore_dof_values();
    }

    // Otherwise we can do predictions in a more object oriented way using
    // whatever timestepper the data provides (this is the normal case).
    else
    {
      // Calculate all predictions in the "master" mesh
      Mesh_pt->calculate_predictions();

      // Calculate predictions for global data with their own timesteppers
      unsigned Nglobal = Global_data_pt.size();
      for (unsigned iglobal = 0; iglobal < Nglobal; iglobal++)
      {
        Global_data_pt[iglobal]->time_stepper_pt()->calculate_predicted_values(
          Global_data_pt[iglobal]);
      }
    }

    // If requested then copy the predicted value into the current time data
    // slots, ready for the newton solver to use as an initial guess.
    if (use_predictor_values_as_initial_guess())
    {
      // Not sure I know enough about distributed problems to implement
      // this. Probably you just need to loop over ndof_local or something,
      // but I can't really test it...
#ifdef OOMPH_HAS_MPI
      if (distributed())
      {
        throw OomphLibError("Not yet implemented for distributed problems",
                            OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

      // With multiple time steppers this is much more complex becuase you
      // need to check the time stepper for each data to get the
      // predictor_storage_index(). Do-able if you need it though.
      if (Time_stepper_pt.size() != 1)
      {
        std::string err = "Not implemented for multiple time steppers";
        throw OomphLibError(
          err, OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
      }

      // Get predicted values
      DoubleVector predicted_dofs;
      get_dofs(time_stepper_pt()->predictor_storage_index(), predicted_dofs);

      // Update dofs at current step
      for (unsigned i = 0; i < ndof(); i++)
      {
        dof(i) = predicted_dofs[i];
      }
    }
  }

  //======================================================================
  ///\short Enable recycling of the mass matrix in explicit timestepping
  /// schemes. Useful for timestepping on fixed meshes when you want
  /// to avoid the linear solve phase.
  //=====================================================================
  void Problem::enable_mass_matrix_reuse()
  {
    Mass_matrix_reuse_is_enabled = true;
    Mass_matrix_has_been_computed = false;

    // If we have a discontinuous formulation set the elements to reuse
    // their own mass matrices
    if (Discontinuous_element_formulation)
    {
      const unsigned n_element = Problem::mesh_pt()->nelement();
      // Loop over the other elements
      for (unsigned e = 0; e < n_element; e++)
      {
        // Cache the element
        DGElement* const elem_pt =
          dynamic_cast<DGElement*>(Problem::mesh_pt()->element_pt(e));
        elem_pt->enable_mass_matrix_reuse();
      }
    }
  }

  //======================================================================
  /// \short Turn off the recyling of the mass matrix in explicit
  /// time-stepping schemes
  //======================================================================
  void Problem::disable_mass_matrix_reuse()
  {
    Mass_matrix_reuse_is_enabled = false;
    Mass_matrix_has_been_computed = false;

    // If we have a discontinuous formulation set the element-level
    // function
    if (Discontinuous_element_formulation)
    {
      const unsigned n_element = Problem::mesh_pt()->nelement();
      // Loop over the other elements
      for (unsigned e = 0; e < n_element; e++)
      {
        // Cache the element
        DGElement* const elem_pt =
          dynamic_cast<DGElement*>(Problem::mesh_pt()->element_pt(e));
        elem_pt->disable_mass_matrix_reuse();
      }
    }
  }


  //=========================================================================
  /// Copy Data values, nodal positions etc from specified problem.
  /// Note: This is not a copy constructor. We assume that the current
  /// and the "original" problem have both been created by calling
  /// the same problem constructor so that all Data objects,
  /// time steppers etc. in the two problems are completely independent.
  /// This function copies the nodal, internal and global values
  /// and the time parameters from the original problem into "this"
  /// one. This functionality is required, e.g. for
  /// multigrid computations.
  //=========================================================================
  void Problem::copy(Problem* orig_problem_pt)
  {
    // Copy time
    //----------

    // Flag to indicate that orig problem is unsteady problem
    bool unsteady_flag = (orig_problem_pt->time_pt() != 0);

    // Copy current time and previous time increments for proper unsteady run
    if (unsteady_flag)
    {
      oomph_info << "Copying an unsteady problem." << std::endl;
      // Current time
      this->time_pt()->time() = orig_problem_pt->time_pt()->time();
      // Timesteps
      unsigned n_dt = orig_problem_pt->time_pt()->ndt();
      time_pt()->resize(n_dt);
      for (unsigned i = 0; i < n_dt; i++)
      {
        time_pt()->dt(i) = orig_problem_pt->time_pt()->dt(i);
      }

      // Find out how many timesteppers there are
      unsigned n_time_steppers = ntime_stepper();

      // Loop over them all and set the weights
      for (unsigned i = 0; i < n_time_steppers; i++)
      {
        time_stepper_pt(i)->set_weights();
      }
    }

    // Copy nodes
    //-----------

    // Loop over submeshes:
    unsigned nmesh = nsub_mesh();
    if (nmesh == 0) nmesh = 1;
    for (unsigned m = 0; m < nmesh; m++)
    {
      // Find number of nodes in present mesh
      unsigned long n_node = mesh_pt(m)->nnode();

      // Check # of nodes:
      unsigned long n_node_orig = orig_problem_pt->mesh_pt(m)->nnode();
      if (n_node != n_node_orig)
      {
        std::ostringstream error_message;
        error_message << "Number of nodes in copy " << n_node
                      << " not equal to the number in the original "
                      << n_node_orig << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Loop over the nodes
      for (unsigned long i = 0; i < n_node; i++)
      {
        // Try to cast to elastic node
        SolidNode* el_node_pt =
          dynamic_cast<SolidNode*>(mesh_pt(m)->node_pt(i));
        if (el_node_pt != 0)
        {
          SolidNode* el_node_orig_pt =
            dynamic_cast<SolidNode*>(orig_problem_pt->mesh_pt(m)->node_pt(i));
          el_node_pt->copy(el_node_orig_pt);
        }
        else
        {
          mesh_pt(m)->node_pt(i)->copy(orig_problem_pt->mesh_pt(m)->node_pt(i));
        }
      }
    }


    // Copy global data:
    //------------------

    // Number of global data
    unsigned n_global = Global_data_pt.size();

    // Check # of nodes in orig problem
    unsigned long n_global_orig = orig_problem_pt->nglobal_data();
    if (n_global != n_global_orig)
    {
      std::ostringstream error_message;
      error_message << "Number of global data in copy " << n_global
                    << " not equal to the number in the original "
                    << n_global_orig << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    for (unsigned iglobal = 0; iglobal < n_global; iglobal++)
    {
      Global_data_pt[iglobal]->copy(orig_problem_pt->global_data_pt(iglobal));
    }


    // Copy internal data of elements:
    //--------------------------------

    // Loop over submeshes:
    for (unsigned m = 0; m < nmesh; m++)
    {
      // Loop over elements and deal with internal data
      unsigned n_element = mesh_pt(m)->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        GeneralisedElement* el_pt = mesh_pt(m)->element_pt(e);
        unsigned n_internal = el_pt->ninternal_data();
        if (n_internal > 0)
        {
          // Check # of internals :
          unsigned long n_internal_orig =
            orig_problem_pt->mesh_pt(m)->element_pt(e)->ninternal_data();
          if (n_internal != n_internal_orig)
          {
            std::ostringstream error_message;
            error_message << "Number of internal data in copy " << n_internal
                          << " not equal to the number in the original "
                          << n_internal_orig << std::endl;

            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
          for (unsigned i = 0; i < n_internal; i++)
          {
            el_pt->internal_data_pt(i)->copy(
              orig_problem_pt->mesh_pt(m)->element_pt(e)->internal_data_pt(i));
          }
        }
      }
    }
  }

  //=========================================================================
  /// Make and return a pointer to the copy of the problem. A virtual
  /// function that must be filled in by the user is they wish to perform
  /// adaptive refinement in bifurcation tracking or in multigrid problems.
  /// ALH: WILL NOT BE NECESSARY IN BIFURCATION TRACKING IN LONG RUN...
  //=========================================================================
  Problem* Problem::make_copy()
  {
    std::ostringstream error_stream;
    error_stream
      << "This function must be overloaded in your specific problem, and must\n"
      << "create an exact copy of your problem. Usually this will be achieved\n"
      << "by a call to the constructor with exactly the same arguments as "
         "used\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


  //=========================================================================
  /// Dump refinement pattern of all refineable meshes and all  generic
  /// Problem data to file for restart.
  //=========================================================================
  void Problem::dump(std::ofstream& dump_file) const
  {
    // Number of submeshes?
    unsigned n_mesh = nsub_mesh();

    dump_file << std::max(unsigned(1), n_mesh) << " # number of (sub)meshes "
              << std::endl;

    // Single mesh:
    //------------
    if (n_mesh == 0)
    {
      // Dump level of refinement before pruning
      if (TreeBasedRefineableMeshBase* mmesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(0)))
      {
        dump_file << mmesh_pt->uniform_refinement_level_when_pruned()
                  << " # uniform refinement when pruned " << std::endl;
      }
      else
      {
        dump_file << 0 << " # (fake) uniform refinement when pruned "
                  << std::endl;
      }
      dump_file << 9999 << " # test flag for end of sub-meshes " << std::endl;
    }

    // Multiple submeshes
    //------------------
    else
    {
      // Loop over submeshes to dump level of refinement before pruning
      for (unsigned imesh = 0; imesh < n_mesh; imesh++)
      {
        if (TreeBasedRefineableMeshBase* mmesh_pt =
              dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(imesh)))
        {
          dump_file << mmesh_pt->uniform_refinement_level_when_pruned()
                    << " # uniform refinement when pruned " << std::endl;
        }
        else
        {
          dump_file << 0 << " # (fake) uniform refinement when pruned "
                    << std::endl;
        }
      }
      dump_file << 9999 << " # test flag for end of sub-meshes " << std::endl;
    }

#ifdef OOMPH_HAS_MPI

    const int my_rank = this->communicator_pt()->my_rank();

    // Record destination of all base elements
    unsigned n = Base_mesh_element_pt.size();
    Vector<int> local_base_element_processor(n, -1);
    Vector<int> base_element_processor(n, -1);
    for (unsigned e = 0; e < n; e++)
    {
      GeneralisedElement* el_pt = Base_mesh_element_pt[e];
      if (el_pt != 0)
      {
        if (!el_pt->is_halo())
        {
          local_base_element_processor[e] = my_rank;
        }
      }
    }


    // Get target for all base elements by reduction
    if (Problem_has_been_distributed)
    {
      // Check that the base elements have been associated to a processor
      // (the Base_mesh_elemen_pt is only used for structured meshes,
      // therefore, if there are no ustructured meshes as part of the
      // problem this container will be empty)
      if (n > 0)
      {
        MPI_Allreduce(&local_base_element_processor[0],
                      &base_element_processor[0],
                      n,
                      MPI_INT,
                      MPI_MAX,
                      this->communicator_pt()->mpi_comm());
      }
    }
    else
    {
      // All the same...
      base_element_processor = local_base_element_processor;
    }


    dump_file << n << " # Number of base elements; partitioning follows.\n";
    for (unsigned e = 0; e < n; e++)
    {
      dump_file << base_element_processor[e] << "\n";
    }
    dump_file << "8888 #test flag for end of base element distribution\n";

#endif

    // Single mesh:
    //------------
    if (n_mesh == 0)
    {
      // Dump single mesh refinement pattern (if mesh is refineable)
      if (TreeBasedRefineableMeshBase* mmesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(0)))
      {
        mmesh_pt->dump_refinement(dump_file);
      }
#ifdef OOMPH_HAS_TRIANGLE_LIB
      // Dump triangle mesh TriangulateIO which represents mesh topology
      TriangleMeshBase* mmesh_pt = dynamic_cast<TriangleMeshBase*>(mesh_pt(0));
      if (mmesh_pt != 0 && mmesh_pt->use_triangulateio_restart())
      {
#ifdef OOMPH_HAS_MPI
        // Check if the mesh is distributed, if that is the case then
        // additional info. needs to be saved
        if (mmesh_pt->is_mesh_distributed())
        {
          // Dump the info. related with the distribution of the mesh
          mmesh_pt->dump_distributed_info_for_restart(dump_file);
        }
#endif
        mmesh_pt->dump_triangulateio(dump_file);
      }
#endif
    }

    // Multiple submeshes
    //------------------
    else
    {
      // Loop over submeshes
      for (unsigned imesh = 0; imesh < n_mesh; imesh++)
      {
        // Dump single mesh refinement pattern (if mesh is refineable)
        if (TreeBasedRefineableMeshBase* mmesh_pt =
              dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(imesh)))
        {
          mmesh_pt->dump_refinement(dump_file);
        }
#ifdef OOMPH_HAS_TRIANGLE_LIB
        // Dump triangle mesh TriangulateIO which represents mesh topology
        TriangleMeshBase* mmesh_pt =
          dynamic_cast<TriangleMeshBase*>(mesh_pt(imesh));
        if (mmesh_pt != 0 && mmesh_pt->use_triangulateio_restart())
        {
#ifdef OOMPH_HAS_MPI
          // Check if the mesh is distributed, if that is the case then
          // additional info. needs to be saved
          if (mmesh_pt->is_mesh_distributed())
          {
            // Dump the info. related with the distribution of the mesh
            mmesh_pt->dump_distributed_info_for_restart(dump_file);
          }
#endif
          mmesh_pt->dump_triangulateio(dump_file);
        }
#endif
      } // End of loop over submeshes
    }


    // Dump time
    // ---------

    // Flag to indicate unsteady run
    bool unsteady_flag = (time_pt() != 0);
    dump_file << unsteady_flag << " # bool flag for unsteady" << std::endl;

    // Current time and previous time increments for proper unsteady run
    if (unsteady_flag)
    {
      // Current time
      dump_file << time_pt()->time() << " # Time " << std::endl;
      // Timesteps
      unsigned n_dt = time_pt()->ndt();
      dump_file << n_dt << " # Number of timesteps " << std::endl;
      for (unsigned i = 0; i < n_dt; i++)
      {
        dump_file << time_pt()->dt(i) << " # dt " << std::endl;
      }
    }
    // Dummy time and previous time increments for steady run
    else
    {
      // Current time
      dump_file << "0.0 # Dummy time from steady run " << std::endl;
      // Timesteps
      dump_file << "0 # Dummy number of timesteps from steady run" << std::endl;
    }

    // Loop over submeshes and dump their data
    unsigned nmesh = nsub_mesh();
    if (nmesh == 0) nmesh = 1;
    for (unsigned m = 0; m < nmesh; m++)
    {
      mesh_pt(m)->dump(dump_file);
    }

    // Dump global data

    // Loop over global data
    unsigned Nglobal = Global_data_pt.size();
    dump_file << Nglobal << " # number of global Data items " << std::endl;
    for (unsigned iglobal = 0; iglobal < Nglobal; iglobal++)
    {
      Global_data_pt[iglobal]->dump(dump_file);
      dump_file << std::endl;
    }
  }

  //=========================================================================
  /// Read refinement pattern of all refineable meshes and refine them
  /// accordingly, then read all Data and nodal position info from
  /// file for restart. Return flag to indicate if the restart was from
  /// steady or unsteady solution.
  //=========================================================================
  void Problem::read(std::ifstream& restart_file, bool& unsteady_restart)
  {
    // Check if the file is actually open as it won't be if it doesn't
    // exist! In that case we're almost certainly restarting the run on
    // a larger number of processors than the restart data was produced.
    // Say so and return
    bool restart_file_is_open = true;
    if (!restart_file.is_open())
    {
      std::ostringstream warn_message;
      warn_message << "Restart file isn't open -- I'm assuming that this is\n";
      warn_message << "because we're restarting on a larger number of\n";
      warn_message << "processor than were in use when the restart data was \n";
      warn_message << "dumped.\n";
      OomphLibWarning(
        warn_message.str(), "Problem::read()", OOMPH_EXCEPTION_LOCATION);
      restart_file_is_open = false;
    }

    // Number of (sub)meshes?
    unsigned n_mesh = std::max(unsigned(1), nsub_mesh());

    std::string input_string;

    // Read line up to termination sign
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Read in number of sub-meshes
    unsigned n_submesh_read;
    n_submesh_read = std::atoi(input_string.c_str());

#ifdef PARANOID
    if (restart_file_is_open)
    {
      if (n_submesh_read != n_mesh)
      {
        std::ostringstream error_message;
        error_message
          << "Number of sub-meshes specified in restart file, "
          << n_submesh_read << " doesn't \n match the my number of sub-meshes,"
          << n_mesh << std::endl
          << "Make sure all sub-meshes have been added to the global mesh\n"
          << "when calling the Problem::dump() function.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#else
    // Suppress comiler warnings about non-used variable
    n_submesh_read++;
    n_submesh_read--;
#endif


    // Read levels of refinement before pruning
#ifdef OOMPH_HAS_MPI
    bool refine_and_prune_required = false;
#endif
    Vector<unsigned> nrefinement_for_mesh(n_mesh);
    for (unsigned i = 0; i < n_mesh; i++)
    {
      // Read line up to termination sign
      getline(restart_file, input_string, '#');

      // Ignore rest of line
      restart_file.ignore(80, '\n');

      // Convert
      nrefinement_for_mesh[i] = std::atoi(input_string.c_str());

      // Get pointer to sub-mesh in incarnation as tree-based refineable mesh
      TreeBasedRefineableMeshBase* ref_mesh_pt =
        dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(i));

      // If it's not a tree-based refineable mesh, ignore the following
      if (ref_mesh_pt == 0)
      {
        if (nrefinement_for_mesh[i] != 0)
        {
          std::ostringstream error_stream;
          error_stream << "Nonzero uniform-refinement-when-pruned specified\n"
                       << "even though mesh is not tree-based. Odd. May want\n"
                       << "to check this carefully before disabling this \n"
                       << "warning/error -- most likely if/when we start to\n"
                       << "prune unstructured meshes [though I can't see why\n"
                       << "we would want to do this, given that they are \n"
                       << "currently totally re-generated...]\n";
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
      else
      {
        // Get min and max refinement level
        unsigned local_min_ref = 0;
        unsigned local_max_ref = 0;
        ref_mesh_pt->get_refinement_levels(local_min_ref, local_max_ref);

        // Overall min refinement level over all meshes
        unsigned min_ref = local_min_ref;

#ifdef OOMPH_HAS_MPI
        if (Problem_has_been_distributed)
        {
          // Reconcile between processors: If (e.g. following
          // distribution/pruning) the mesh has no elements on this
          // processor) then ignore its contribution to the poll of
          // max/min refinement levels
          int int_local_min_ref = local_min_ref;
          if (ref_mesh_pt->nelement() == 0)
          {
            int_local_min_ref = INT_MAX;
          }
          int int_min_ref = 0;
          MPI_Allreduce(&int_local_min_ref,
                        &int_min_ref,
                        1,
                        MPI_INT,
                        MPI_MIN,
                        Communicator_pt->mpi_comm());

          // Overall min refinement level over all meshes
          min_ref = unsigned(int_min_ref);
        }
#endif

        // Need to refine less
        if (nrefinement_for_mesh[i] >= min_ref)
        {
          nrefinement_for_mesh[i] -= min_ref;
        }
      }

#ifdef OOMPH_HAS_MPI
      if (nrefinement_for_mesh[i] > 0)
      {
        refine_and_prune_required = true;
      }
#endif
    }


    // Reconcile overall need to refine and prune (even empty
    // processors have to participate in some communication!)
#ifdef OOMPH_HAS_MPI
    if (Problem_has_been_distributed)
    {
      unsigned local_req_flag = 0;
      unsigned req_flag = 0;
      if (refine_and_prune_required)
      {
        local_req_flag = 1;
      }
      MPI_Allreduce(&local_req_flag,
                    &req_flag,
                    1,
                    MPI_UNSIGNED,
                    MPI_MAX,
                    Communicator_pt->mpi_comm());
      refine_and_prune_required = false;
      if (req_flag == 1)
      {
        refine_and_prune_required = true;
      }

      // If refine and prune is required make number of uniform
      // refinements for each mesh consistent otherwise code
      // hangs on "empty" processors for which no restart file exists
      if (refine_and_prune_required)
      {
        // This is what we have locally
        Vector<unsigned> local_nrefinement_for_mesh(nrefinement_for_mesh);
        // Synchronise over all processors with max operation
        MPI_Allreduce(&local_nrefinement_for_mesh[0],
                      &nrefinement_for_mesh[0],
                      n_mesh,
                      MPI_UNSIGNED,
                      MPI_MAX,
                      Communicator_pt->mpi_comm());

#ifdef PARANOID
        // Check it: Reconciliation should only be required for
        // for processors on which no restart file was opened and
        // for which the meshes are therefore empty
        bool fail = false;
        std::ostringstream error_message;
        error_message << "Number of uniform refinements was not consistent \n"
                      << "for following meshes during restart on processor \n"
                      << "on which restart file could be opened:\n";
        for (unsigned i = 0; i < n_mesh; i++)
        {
          if ((local_nrefinement_for_mesh[i] != nrefinement_for_mesh[i]) &&
              restart_file_is_open)
          {
            fail = true;
            error_message << "Sub-mesh: " << i << "; local nrefinement: "
                          << local_nrefinement_for_mesh[i] << " "
                          << "; global/synced nrefinement: "
                          << nrefinement_for_mesh[i] << "\n";
          }
        }
        if (fail)
        {
          OomphLibWarning(
            error_message.str(), "Problem::read()", OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }
    }
#endif

    // Read line up to termination sign
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Check flag that indicates that we've read the final data
    unsigned tmp;
    tmp = std::atoi(input_string.c_str());

#ifdef PARANOID
    if (restart_file_is_open)
    {
      if (tmp != 9999)
      {
        std::ostringstream error_message;
        error_message
          << "Error in reading restart data: Uniform refinement when pruned \n"
          << "flags should be followed by 9999.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }

#else
    // Suppress comiler warnings about non-used variable
    tmp++;
    tmp--;
#endif


#ifdef OOMPH_HAS_MPI

    // Refine and prune if required
    if (refine_and_prune_required)
    {
      refine_uniformly(nrefinement_for_mesh);
      prune_halo_elements_and_nodes();
    }

    // target_domain_for_local_non_halo_element[e] contains the number
    // of the domain [0,1,...,nproc-1] to which non-halo element e on THE
    // CURRENT PROCESSOR ONLY has been assigned. The order of the non-halo
    // elements is the same as in the Problem's mesh, with the halo
    // elements being skipped.
    Vector<unsigned> target_domain_for_local_non_halo_element;

    // If a restart file has been generated using code compiled without MPI
    // then it will not have any of the base element data.
    // If we try to read in that file with code that has been compied using
    // MPI, even if running only one processor, then it will fail here.
    // The ideal fix is to edit the restart file so that it contains the two
    // lines
    //
    // 0 # Number of base elements; partitioning follows.
    // 8888 # Test flag for end of base element distribution
    //
    // after the end of the sub-meshes, but before the number of elements
    // However, we can determine that this is the problem if n_base = 0,
    // so there is a little bit of logic below to catch this case

    // Store current location in the file (before we are about to read
    // in either the base mesh or number of elements of the first mesh)
    std::streampos position_before_base_element = restart_file.tellg();
    // Boolean flag used to set whether to read in base element info
    bool read_in_base_element_info = true;

    // Read line up to termination sign
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Get number of base elements as recorded
    unsigned n_base_element_read_in = atoi(input_string.c_str());
    unsigned nbase = Base_mesh_element_pt.size();
    if (restart_file_is_open)
    {
      if (n_base_element_read_in != nbase)
      {
        // If we have zero base elements the problem could be that the
        // restart file was generated without MPI. Issue a warning
        // and continue anyway
        if (nbase == 0)
        {
          std::ostringstream warn_message;
          warn_message
            << "The number of base elements in the mesh is 0,\n"
            << " but the restart file indicates that there are "
            << n_base_element_read_in << ".\n"
            << "This could be because the restart file was \n"
            << "generated by using code without MPI.\n"
            << "\n"
            << "The best fix is to include two additional lines\n"
            << "in the restart file: \n\n"
            << "0 # Number of base elements; partitioning follows.\n"
            << "8888 # Test flag for end of base element distribution\n"
            << "\n"
            << "These lines go after the flag 9999 that indicates\n"
            << "the end of the submesh information.\n"
            << "\n"
            << "The file will now continue to be read assuming that\n"
            << "the base element information is not present.\n"
            << "If you get strange results then please look carefully\n"
            << "at the restart file. The safest thing to do is to \n"
            << "ensure that the restart file was generated by code\n"
            << "compiled and run with the same parallel options.\n";
          OomphLibWarning(warn_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
          // Set the skip flag to true
          // and rewind the file pointer
          read_in_base_element_info = false;
          restart_file.seekg(position_before_base_element);
        }
        // Otherwise throw a hard error
        else
        {
          std::ostringstream error_message;
          error_message << "About to read " << n_base_element_read_in
                        << " base elements \n"
                        << "though we only have " << nbase
                        << " base elements in mesh.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
    }

    // Read in the remaning base element information, if necessary
    if (read_in_base_element_info == true)
    {
      // Read in target_domain_for_base_element[e] for all base elements
      Vector<unsigned> target_domain_for_base_element(nbase);
      for (unsigned e = 0; e < nbase; e++)
      {
        // Read line
        getline(restart_file, input_string);

        // Get target domain
        target_domain_for_base_element[e] = atoi(input_string.c_str());
      }

      // Read line up to termination sign
      getline(restart_file, input_string, '#');

      // Ignore rest of line
      restart_file.ignore(80, '\n');

      // Check flag that indicates that we've read the final data
      tmp = std::atoi(input_string.c_str());


#ifdef PARANOID
      if (restart_file_is_open)
      {
        if (tmp != 8888)
        {
          std::ostringstream error_message;
          error_message
            << "Error in reading restart data: Target proc for base elements \n"
            << "should be followed by 8888.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
#endif

      // Loop over all elements (incl. any FaceElements) and assign
      // target domain for all local non-halo elements and check if
      // load balancing is required -- no need to do this if problem is
      // not distributed.
      unsigned load_balance_required_flag = 0;
      if (Problem_has_been_distributed)
      {
        // Working with TreeBasedRefineableMeshBase mesh
        unsigned local_load_balance_required_flag = 0;
        if (dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(0)))
        {
          const int my_rank = this->communicator_pt()->my_rank();
          unsigned nel = mesh_pt()->nelement();
          for (unsigned e = 0; e < nel; e++)
          {
            GeneralisedElement* el_pt = mesh_pt()->element_pt(e);
            if (!el_pt->is_halo())
            {
              // Get element number (plus one) in base element enumeration
              unsigned el_number_in_base_mesh_plus_one =
                Base_mesh_element_number_plus_one[el_pt];

              // If it's zero then we haven't found it, it may be a FaceElement
              // (in which case we move it to the same processor as its bulk
              // element
              if (el_number_in_base_mesh_plus_one == 0)
              {
                FaceElement* face_el_pt = dynamic_cast<FaceElement*>(el_pt);
                if (face_el_pt != 0)
                {
                  // Get corresponding bulk element
                  FiniteElement* bulk_el_pt = face_el_pt->bulk_element_pt();

                  // Use its element number (plus one) in base element
                  // enumeration
                  el_number_in_base_mesh_plus_one =
                    Base_mesh_element_number_plus_one[bulk_el_pt];

                  // If this is zero too we have a problem
                  if (el_number_in_base_mesh_plus_one == 0)
                  {
                    throw OomphLibError(
                      "el_number_in_base_mesh_plus_one=0 for bulk",
                      "Problem::read()",
                      OOMPH_EXCEPTION_LOCATION);
                  }
                }
              }

              // If we've made it here then we're not dealing with a
              // FaceElement but with an element that doesn't exist locally
              // --> WTF?
              if (el_number_in_base_mesh_plus_one == 0)
              {
                throw OomphLibError("el_number_in_base_mesh_plus_one=0",
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }

              // Assign target domain for next local non-halo element in
              // the order in which it's encountered in the global mesh
              target_domain_for_local_non_halo_element.push_back(
                target_domain_for_base_element[el_number_in_base_mesh_plus_one -
                                               1]);

              // Do elements on this processor to be moved elsewhere?
              if (int(target_domain_for_base_element
                        [el_number_in_base_mesh_plus_one - 1]) != my_rank)
              {
                local_load_balance_required_flag = 1;
              }
            }
          }

        } // if (working with TreeBasedRefineableMeshBase mesh)

        // Get overall need to load balance by max
        MPI_Allreduce(&local_load_balance_required_flag,
                      &load_balance_required_flag,
                      1,
                      MPI_UNSIGNED,
                      MPI_MAX,
                      this->communicator_pt()->mpi_comm());
      }

      // Do we need to load balance?
      if (load_balance_required_flag == 1)
      {
        oomph_info << "Doing load balancing after pruning\n";
        DocInfo doc_info;
        doc_info.disable_doc();
        bool report_stats = false;
        load_balance(
          doc_info, report_stats, target_domain_for_local_non_halo_element);
        oomph_info << "Done load balancing after pruning\n";
      }
      else
      {
        oomph_info << "No need for load balancing after pruning\n";
      }
    } // End of read in base element information
#endif


    // Boolean to record if any unstructured bulk meshes have
    // been read in (and therefore completely re-generated, with new
    // elements and nodes) from disk
    bool have_read_unstructured_mesh = false;

    // Call the actions before adaptation
    actions_before_adapt();

    // If there are unstructured meshes in the problem we need
    // to strip out any face elements that are attached to them
    // because restart of unstructured meshes re-creates their elements
    // and nodes from scratch, leading to dangling pointers from the
    // face elements to the old elements and nodes. This function is
    // virtual and (practically) empty in the Problem base class
    // but toggles a flag to indicate that it has been called. We can then
    // issue a warning below, prompting the user to consider overloading it
    // if the problem is found to contain unstructured bulk meshes.
    // Warning can be ignored if the bulk mesh is not associated with any
    // face elements.
    Empty_actions_before_read_unstructured_meshes_has_been_called = false;
    actions_before_read_unstructured_meshes();

    // Update number of submeshes
    n_mesh = nsub_mesh();

    // Single mesh:
    //------------
    if (n_mesh == 0)
    {
      // Refine single mesh (if it's refineable)
      if (TreeBasedRefineableMeshBase* mmesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(0)))
      {
        // When we get in here the problem has been constructed
        // by the constructor and the mesh is its original unrefined
        // form.
        // RefineableMeshBase::refine(...) reads the refinement pattern from the
        // specified file and performs refinements until the mesh has
        // reached the same level of refinement as the mesh that existed
        // when the problem was dumped to disk.
        mmesh_pt->refine(restart_file);
      }
#ifdef OOMPH_HAS_TRIANGLE_LIB
      // Regenerate mesh from triangulate IO if it's a triangular mesh
      TriangleMeshBase* mmesh_pt = dynamic_cast<TriangleMeshBase*>(mesh_pt(0));
      if (mmesh_pt != 0 && mmesh_pt->use_triangulateio_restart())
      {
#ifdef OOMPH_HAS_MPI
        // Check if the mesh is distributed, if that is the case then
        // additional info. needs to be read
        if (mmesh_pt->is_mesh_distributed())
        {
          // Dump the info. related with the distribution of the mesh
          mmesh_pt->read_distributed_info_for_restart(restart_file);
        }
#endif
        // The function reads the TriangulateIO data structure from the dump
        // file and then completely regenerates the mesh using the
        // data structure
        mmesh_pt->remesh_from_triangulateio(restart_file);
        have_read_unstructured_mesh = true;
#ifdef OOMPH_HAS_MPI
        // Check if the mesh is distributed, if that is the case then we
        // need to re-establish the halo/haloed scheme (similar as in the
        // RefineableTriangleMesh::adapt() method)
        if (mmesh_pt->is_mesh_distributed())
        {
          mmesh_pt->reestablish_distribution_info_for_restart(
            this->communicator_pt(), restart_file);
        }
#endif
        // Still left to update the polylines representation, that is performed
        // later since the nodes positions may still change when reading info.
        // for the mesh, see below
      }
#endif
    }

    // Multiple submeshes
    //------------------
    else
    {
      // Loop over submeshes
      for (unsigned imesh = 0; imesh < n_mesh; imesh++)
      {
        // Refine single mesh (if its refineable)
        if (TreeBasedRefineableMeshBase* mmesh_pt =
              dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(imesh)))
        {
          // When we get in here the problem has been constructed
          // by the constructor and the mesh is its original unrefined
          // form.
          // RefineableMeshBase::refine(...) reads the refinement pattern from
          // the specified file and performs refinements until the mesh has
          // reached the same level of refinement as the mesh that existed
          // when the problem was dumped to disk.
          mmesh_pt->refine(restart_file);
        }
#ifdef OOMPH_HAS_TRIANGLE_LIB
        // Regenerate mesh from triangulate IO if it's a triangular mesh
        TriangleMeshBase* mmesh_pt =
          dynamic_cast<TriangleMeshBase*>(mesh_pt(imesh));
        if (mmesh_pt != 0 && mmesh_pt->use_triangulateio_restart())
        {
#ifdef OOMPH_HAS_MPI
          // Check if the mesh is distributed, if that is the case then
          // additional info. needs to be read
          if (mmesh_pt->is_mesh_distributed())
          {
            // Dump the info. related with the distribution of the mesh
            mmesh_pt->read_distributed_info_for_restart(restart_file);
          }
#endif
          // The function reads the TriangulateIO data structure from the dump
          // file and then completely regenerates the mesh using the
          // data structure
          mmesh_pt->remesh_from_triangulateio(restart_file);
          have_read_unstructured_mesh = true;

#ifdef OOMPH_HAS_MPI
          // Check if the mesh is distributed, if that is the case then we
          // need to re-establish the halo/haloed scheme (similar as in the
          // RefineableTriangleMesh::adapt() method)
          if (mmesh_pt->is_mesh_distributed())
          {
            mmesh_pt->reestablish_distribution_info_for_restart(
              this->communicator_pt(), restart_file);
          }
#endif
          // Still left to update the polylines representation, that is
          // performed later since the nodes positions may still change when
          // reading info. for the mesh, see below
        }
#endif
      } // End of loop over submeshes


      // Rebuild the global mesh
      rebuild_global_mesh();
    }

    // Any actions after adapt
    actions_after_adapt();

    // Re-attach face elements (or whatever else needs to be done
    // following the total re-generation of the unstructured meshes
    Empty_actions_after_read_unstructured_meshes_has_been_called = false;
    actions_after_read_unstructured_meshes();


    // Issue warning:
    if (!Suppress_warning_about_actions_before_read_unstructured_meshes)
    {
      if (have_read_unstructured_mesh)
      {
        if (Empty_actions_before_read_unstructured_meshes_has_been_called ||
            Empty_actions_after_read_unstructured_meshes_has_been_called)
        {
          std::ostringstream warn_message;
          warn_message
            << "I've just read in some unstructured meshes and have, in\n"
            << "the process, totally re-generated their nodes and elements.\n"
            << "This may create dangling pointers that still point to the\n"
            << "old nodes and elements, e.g. because FaceElements were\n"
            << "attached to these meshes or pointers to nodes and elements\n"
            << "were stored somewhere. FaceElements should therefore be\n"
            << "removed before reading in these meshes, using an overloaded\n"
            << "version of the function\n\n"
            << "   Problem::actions_before_read_unstructured_meshes()\n\n"
            << "and then re-attached using an overloaded version of\n\n"
            << "   Problem::actions_after_read_unstructured_meshes().\n\n"
            << "The required content of these functions is likely to be "
               "similar\n"
            << "to the Problem::actions_before_adapt() and \n"
            << "Problem::actions_after_adapt() that would be required in\n"
            << "a spatially adaptive computation. If these functions already\n"
            << "exist and perform the required actions, the \n"
            << "actions_before/after_read_unstructured_meshes() functions\n"
            << "can remain empty because the former are called automatically.\n"
            << "In this case, this warning my be suppressed by setting the\n"
            << "public boolean\n\n"
            << "   "
               "Problem::Suppress_warning_about_actions_before_read_"
               "unstructured_meshes\n\n"
            << "to true." << std::endl;
          OomphLibWarning(warn_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
        }
      }
    }

    // Setup equation numbering scheme
    oomph_info << "\nNumber of equations in Problem::read(): "
               << assign_eqn_numbers() << std::endl
               << std::endl;
    // Read time info
    //---------------
    unsigned local_unsteady_restart_flag = 0;
    double local_time = -DBL_MAX;
    unsigned local_n_dt = 0;
#ifdef OOMPH_HAS_MPI
    unsigned local_sync_needed_flag = 0;
#endif
    Vector<double> local_dt;

    if (restart_file.is_open())
    {
      oomph_info << "Restart file exists" << std::endl;
#ifdef OOMPH_HAS_MPI
      local_sync_needed_flag = 0;
#endif
      // Read line up to termination sign
      getline(restart_file, input_string, '#');

      // Ignore rest of line
      restart_file.ignore(80, '\n');

      // Is the restart data from an unsteady run?
      local_unsteady_restart_flag = atoi(input_string.c_str());

      // Read line up to termination sign
      getline(restart_file, input_string, '#');

      // Ignore rest of line
      restart_file.ignore(80, '\n');

      // Read in initial time and set
      local_time = atof(input_string.c_str());

      // Read line up to termination sign
      getline(restart_file, input_string, '#');

      // Ignore rest of line
      restart_file.ignore(80, '\n');

      // Read & set number of timesteps
      local_n_dt = atoi(input_string.c_str());
      local_dt.resize(local_n_dt);

      // Read in timesteps:
      for (unsigned i = 0; i < local_n_dt; i++)
      {
        // Read line up to termination sign
        getline(restart_file, input_string, '#');

        // Ignore rest of line
        restart_file.ignore(80, '\n');

        // Read in initial time and set
        double prev_dt = atof(input_string.c_str());
        local_dt[i] = prev_dt;
      }
    }
    else
    {
      oomph_info << "Restart file does not exist" << std::endl;
#ifdef OOMPH_HAS_MPI
      local_sync_needed_flag = 1;
#endif
    }


    // No prepare global values, possibly via sync
    Vector<double> dt;

    // Do we need to sync?
    unsigned sync_needed_flag = 0;

#ifdef OOMPH_HAS_MPI
    if (Problem_has_been_distributed)
    {
      // Get need to sync by max
      MPI_Allreduce(&local_sync_needed_flag,
                    &sync_needed_flag,
                    1,
                    MPI_UNSIGNED,
                    MPI_MAX,
                    this->communicator_pt()->mpi_comm());
    }
#endif

    // Synchronise
    if (sync_needed_flag == 1)
    {
#ifdef OOMPH_HAS_MPI


#ifdef PARANOID
      if (!Problem_has_been_distributed)
      {
        std::ostringstream error_message;
        error_message << "Synchronisation of temporal restart data \n"
                      << "required even though Problem hasn't been distributed "
                         "-- very odd!\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Get unsteady restart flag by max-based reduction
      unsigned unsteady_restart_flag = 0;
      MPI_Allreduce(&local_unsteady_restart_flag,
                    &unsteady_restart_flag,
                    1,
                    MPI_UNSIGNED,
                    MPI_MAX,
                    this->communicator_pt()->mpi_comm());

      // So, is it an unsteady restart?
      unsteady_restart = false;
      if (unsteady_restart_flag == 1)
      {
        unsteady_restart = true;

        // Get time by max
        double time = -DBL_MAX;
        MPI_Allreduce(&local_time,
                      &time,
                      1,
                      MPI_DOUBLE,
                      MPI_MAX,
                      this->communicator_pt()->mpi_comm());
        time_pt()->time() = time;

        // Get number of timesteps by max-based reduction
        unsigned n_dt = 0;
        MPI_Allreduce(&local_n_dt,
                      &n_dt,
                      1,
                      MPI_UNSIGNED,
                      MPI_MAX,
                      this->communicator_pt()->mpi_comm());

        // Resize whatever needs resizing
        time_pt()->resize(n_dt);
        dt.resize(n_dt);
        if (local_dt.size() == 0)
        {
          local_dt.resize(n_dt, -DBL_MAX);
        }

        // Get timesteps increments by max-based reduction
        MPI_Allreduce(&local_dt[0],
                      &dt[0],
                      n_dt,
                      MPI_DOUBLE,
                      MPI_MAX,
                      this->communicator_pt()->mpi_comm());
      }

#else

      std::ostringstream error_message;
      error_message
        << "Synchronisation of temporal restart data \n"
        << "required even though we don't have mpi support -- very odd!\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);

#endif
    }
    // No sync needed -- just copy across
    else
    {
      unsteady_restart = false;
      if (local_unsteady_restart_flag == 1)
      {
        unsteady_restart = true;
        time_pt()->time() = local_time;
        time_pt()->resize(local_n_dt);
        dt.resize(local_n_dt);
        for (unsigned i = 0; i < local_n_dt; i++)
        {
          dt[i] = local_dt[i];
        }
      }
    }

    // Initialise timestep -- also sets the weights for all timesteppers
    // in the problem.
    if (unsteady_restart) initialise_dt(dt);

    // Loop over submeshes:
    unsigned nmesh = nsub_mesh();
    if (nmesh == 0) nmesh = 1;
    for (unsigned m = 0; m < nmesh; m++)
    {
      //    //---------------------------------------------------------
      //    // Keep this commented out code around to debug restarts
      //    //---------------------------------------------------------
      //    std::ofstream some_file;
      //    char filename[100];
      //    sprintf(filename,"read_mesh%i_on_proc%i.dat",m,
      //            this->communicator_pt()->my_rank());
      //    some_file.open(filename);
      //    mesh_pt(m)->output(some_file);
      //    some_file.close();

      //    sprintf(filename,"read_mesh%i_with_haloes_on_proc%i.dat",m,
      //            this->communicator_pt()->my_rank());
      //    mesh_pt(m)->enable_output_of_halo_elements();
      //    some_file.open(filename);
      //    mesh_pt(m)->output(some_file);
      //    mesh_pt(m)->disable_output_of_halo_elements();
      //    some_file.close();
      //    oomph_info << "Doced mesh " << m << " before reading\n";

      //    sprintf(filename,"read_nodes_mesh%i_on_proc%i.dat",m,
      //            this->communicator_pt()->my_rank());
      //    some_file.open(filename);
      //    unsigned nnod=mesh_pt(m)->nnode();
      //    for (unsigned j=0;j<nnod;j++)
      //     {
      //      Node* nod_pt=mesh_pt(m)->node_pt(j);
      //      unsigned n=nod_pt->ndim();
      //      for (unsigned i=0;i<n;i++)
      //       {
      //        some_file << nod_pt->x(i) << " ";
      //       }
      //      some_file << nod_pt->is_halo() << " "
      //                << nod_pt->nvalue() << " "
      //                << nod_pt->hang_code() << "\n";
      //     }
      //    some_file.close();
      //    oomph_info << "Doced mesh " << m << " before reading\n";
      //    //---------------------------------------------------------
      //    // End keep this commented out code around to debug restarts
      //    //---------------------------------------------------------

      mesh_pt(m)->read(restart_file);

#ifdef OOMPH_HAS_TRIANGLE_LIB
      // Here update the polyline representation if working with
      // triangle base meshes
      if (TriangleMeshBase* mmesh_pt =
            dynamic_cast<TriangleMeshBase*>(mesh_pt(m)))
      {
        // In charge of updating the polylines representation to the
        // current refinement/unrefinement level after restart, it
        // also update the shared boundaries in case of working with a
        // distributed mesh
        mmesh_pt->update_polyline_representation_from_restart();
      }
#endif // #ifdef OOMPH_HAS_TRIANGLE_LIB
    }

    // Read global data:
    //------------------

    // Number of global data
    unsigned Nglobal = Global_data_pt.size();

    // Read line up to termination sign
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Check # of nodes:
    unsigned long check_nglobal = atoi(input_string.c_str());


    if (restart_file_is_open)
    {
      if (check_nglobal != Nglobal)
      {
        std::ostringstream error_message;
        error_message << "The number of global data " << Nglobal
                      << " is not equal to that specified in the input file "
                      << check_nglobal << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }

    for (unsigned iglobal = 0; iglobal < Nglobal; iglobal++)
    {
      Global_data_pt[iglobal]->read(restart_file);
    }
  }

  //===================================================================
  /// Set all timesteps to the same value, dt, and assign
  /// weights for all timesteppers in the problem.
  //===================================================================
  void Problem::initialise_dt(const double& dt)
  {
    // Initialise the timesteps in the Problem's time object
    Time_pt->initialise_dt(dt);

    // Find out how many timesteppers there are
    unsigned n_time_steppers = ntime_stepper();

    // Loop over them all and set the weights
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      time_stepper_pt(i)->set_weights();
      if (time_stepper_pt(i)->adaptive_flag())
      {
        time_stepper_pt(i)->set_error_weights();
      }
    }
  }

  //=========================================================================
  /// Set the value of the timesteps to be equal to the values passed in
  /// a vector and assign weights for all timesteppers in the problem
  //========================================================================
  void Problem::initialise_dt(const Vector<double>& dt)
  {
    // Initialise the timesteps in the Problem's time object
    Time_pt->initialise_dt(dt);

    // Find out how many timesteppers there are
    unsigned n_time_steppers = ntime_stepper();

    // Loop over them all and set the weights
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      time_stepper_pt(i)->set_weights();
      if (time_stepper_pt(i)->adaptive_flag())
      {
        time_stepper_pt(i)->set_error_weights();
      }
    }
  }

  //========================================================
  /// Self-test: Check meshes and global data. Return 0 for OK
  //========================================================
  unsigned Problem::self_test()
  {
    // Initialise
    bool passed = true;

    // Are there any submeshes?
    unsigned Nmesh = nsub_mesh();

    // Just one mesh: Check it
    if (Nmesh == 0)
    {
      if (mesh_pt()->self_test() != 0)
      {
        passed = false;
        oomph_info
          << "\n ERROR: Failed Mesh::self_test() for single mesh in problem"
          << std::endl;
      }
    }
    // Loop over all submeshes and check them
    else
    {
      for (unsigned imesh = 0; imesh < Nmesh; imesh++)
      {
        if (mesh_pt(imesh)->self_test() != 0)
        {
          passed = false;
          oomph_info << "\n ERROR: Failed Mesh::self_test() for mesh imesh"
                     << imesh << std::endl;
        }
      }
    }


    // Check global data
    unsigned Nglobal = Global_data_pt.size();
    for (unsigned iglobal = 0; iglobal < Nglobal; iglobal++)
    {
      if (Global_data_pt[iglobal]->self_test() != 0)
      {
        passed = false;
        oomph_info
          << "\n ERROR: Failed Data::self_test() for global data iglobal"
          << iglobal << std::endl;
      }
    }


#ifdef OOMPH_HAS_MPI

    if (Problem_has_been_distributed)
    {
      // Note: This throws an error if it fails so no return is required.
      DocInfo tmp_doc_info;
      tmp_doc_info.disable_doc();
      check_halo_schemes(tmp_doc_info);
    }

#endif

    // Return verdict
    if (passed)
    {
      return 0;
    }
    else
    {
      return 1;
    }
  }

  //====================================================================
  /// A function that is used to adapt a bifurcation-tracking
  /// problem, which requires separate interpolation of the
  /// associated eigenfunction. The error measure is chosen to be
  /// a suitable combination of the errors in the base flow and the
  /// eigenfunction. The bifurcation type is passed as an argument
  //=====================================================================
  void Problem::bifurcation_adapt_helper(unsigned& n_refined,
                                         unsigned& n_unrefined,
                                         const unsigned& bifurcation_type,
                                         const bool& actually_adapt)
  {
    // Storage for eigenfunction from the problem
    Vector<DoubleVector> eigenfunction;
    // Get the eigenfunction from the problem
    this->get_bifurcation_eigenfunction(eigenfunction);

    // Get the bifurcation parameter
    double* parameter_pt = this->bifurcation_parameter_pt();

    // Get the frequency parameter if tracking a Hopf bifurcation
    double omega = 0.0;
    // If we're tracking a Hopf then also get the frequency
    if (bifurcation_type == 3)
    {
      omega = dynamic_cast<HopfHandler*>(assembly_handler_pt())->omega();
    }

    // If we're tracking a Pitchfork get the slack parameter (Hack)
    double sigma = 0.0;
    if (bifurcation_type == 2)
    {
      sigma = this->dof(this->ndof() - 1);
    }

    // We can now deactivate the bifurcation tracking in the problem
    // to restore the degrees of freedom to the unaugmented value
    this->deactivate_bifurcation_tracking();

    // Next, we create copies of the present problem
    // The number of copies depends on the number of eigenfunctions
    // One copy for each eigenfunction
    const unsigned n_copies = eigenfunction.size();
    Copy_of_problem_pt.resize(n_copies);

    // Loop over the number of copies
    for (unsigned c = 0; c < n_copies; c++)
    {
      // If we don't already have a copy
      if (Copy_of_problem_pt[c] == 0)
      {
        // Create the copy
        Copy_of_problem_pt[c] = this->make_copy();

        // Refine the copy to the same level as the current problem

        // Find number of submeshes
        const unsigned N_mesh = Copy_of_problem_pt[c]->nsub_mesh();
        // If there is only one mesh
        if (N_mesh == 0)
        {
          // Can we refine the mesh
          if (TreeBasedRefineableMeshBase* mmesh_pt =
                dynamic_cast<TreeBasedRefineableMeshBase*>(
                  Copy_of_problem_pt[c]->mesh_pt(0)))
          {
            // Is the adapt flag set
            if (mmesh_pt->is_adaptation_enabled())
            {
              // Now get the original problem's mesh if it's refineable
              if (TreeBasedRefineableMeshBase* original_mesh_pt =
                    dynamic_cast<TreeBasedRefineableMeshBase*>(
                      this->mesh_pt(0)))
              {
                mmesh_pt->refine_base_mesh_as_in_reference_mesh(
                  original_mesh_pt);
              }
              else
              {
                oomph_info
                  << "Info/Warning: Mesh in orginal problem is not refineable."
                  << std::endl;
              }
            }
            else
            {
              oomph_info << "Info/Warning: Mesh adaptation is disabled in copy."
                         << std::endl;
            }
          }
          else
          {
            oomph_info << "Info/Warning: Mesh cannot be adapted in copy."
                       << std::endl;
          }
        } // End of single mesh case
        // Otherwise loop over the submeshes
        else
        {
          for (unsigned m = 0; m < N_mesh; m++)
          {
            // Can we refine the submesh
            if (TreeBasedRefineableMeshBase* mmesh_pt =
                  dynamic_cast<TreeBasedRefineableMeshBase*>(
                    Copy_of_problem_pt[c]->mesh_pt(m)))
            {
              // Is the adapt flag set
              if (mmesh_pt->is_adaptation_enabled())
              {
                // Now get the original problem's mesh
                if (TreeBasedRefineableMeshBase* original_mesh_pt =
                      dynamic_cast<TreeBasedRefineableMeshBase*>(
                        this->mesh_pt(m)))
                {
                  mmesh_pt->refine_base_mesh_as_in_reference_mesh(
                    original_mesh_pt);
                }
                else
                {
                  oomph_info << "Info/Warning: Mesh in orginal problem is not "
                                "refineable."
                             << std::endl;
                }
              }
              else
              {
                oomph_info
                  << "Info/Warning: Mesh adaptation is disabled in copy."
                  << std::endl;
              }
            }
            else
            {
              oomph_info << "Info/Warning: Mesh cannot be adapted in copy."
                         << std::endl;
            }
          }
          // rebuild the global mesh in the copy
          Copy_of_problem_pt[c]->rebuild_global_mesh();

        } // End of multiple mesh case

        // Must call actions after adapt
        Copy_of_problem_pt[c]->actions_after_adapt();

        // Assign the equation numbers to the copy (quietly)
        (void)Copy_of_problem_pt[c]->assign_eqn_numbers();
      }
    } // End of creation of copies


    // Now check some numbers
    for (unsigned c = 0; c < n_copies; c++)
    {
      // Check that the dofs match for each copy
#ifdef PARANOID
      // If the problems don't match then complain
      if (Copy_of_problem_pt[c]->ndof() != this->ndof())
      {
        std::ostringstream error_stream;
        error_stream << "Number of unknowns in the problem copy " << c << " "
                     << "not equal to number in the original:\n"
                     << this->ndof() << " (original) "
                     << Copy_of_problem_pt[c]->ndof() << " (copy)\n";

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Assign the eigenfunction(s) to the copied problems
      Copy_of_problem_pt[c]->assign_eigenvector_to_dofs(eigenfunction[c]);
      // Set all pinned values to zero
      Copy_of_problem_pt[c]->set_pinned_values_to_zero();
    }

    // Symmetrise the problem if we are solving a pitchfork
    if (bifurcation_type == 2)
    {
      Copy_of_problem_pt[0]
        ->symmetrise_eigenfunction_for_adaptive_pitchfork_tracking();
    }

    // Find error estimates based on current problem and eigenproblem
    // Now we need to get the error estimates for both problems.
    Vector<Vector<double>> base_error, eigenfunction_error;
    this->get_all_error_estimates(base_error);
    // Loop over the copies
    for (unsigned c = 0; c < n_copies; c++)
    {
      // Get the error estimates for the copy
      Copy_of_problem_pt[c]->get_all_error_estimates(eigenfunction_error);

      // Find the number of meshes
      unsigned n_mesh = base_error.size();

#ifdef PARANOID
      if (n_mesh != eigenfunction_error.size())
      {
        std::ostringstream error_stream;
        error_stream << "Problems do not have the same number of meshes\n"
                     << "Base : " << n_mesh
                     << " : Eigenproblem : " << eigenfunction_error.size()
                     << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      for (unsigned m = 0; m < n_mesh; m++)
      {
        // Check the number of elements is the same
        unsigned n_element = base_error[m].size();
#ifdef PARANOID
        if (n_element != eigenfunction_error[m].size())
        {
          std::ostringstream error_stream;
          error_stream << "Mesh " << m
                       << " does not have the same number of elements in the "
                          "two problems:\n"
                       << "Base: " << n_element
                       << " :  Eigenproblem: " << eigenfunction_error[m].size()
                       << "\n";
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        // Now add all the error esimates together
        for (unsigned e = 0; e < n_element; e++)
        {
          // Add the error estimates (lazy)
          base_error[m][e] += eigenfunction_error[m][e];
        }
      }
    } // End of loop over copies

    // Then refine all problems based on the combined measure
    // if we are actually adapting (not just estimating the errors)
    if (actually_adapt)
    {
      this->adapt_based_on_error_estimates(n_refined, n_unrefined, base_error);
      for (unsigned c = 0; c < n_copies; c++)
      {
        Copy_of_problem_pt[c]->adapt_based_on_error_estimates(
          n_refined, n_unrefined, base_error);
      }
      // Symmetrise the problem (again) if we are solving for a pitchfork
      if (bifurcation_type == 2)
      {
        Copy_of_problem_pt[0]
          ->symmetrise_eigenfunction_for_adaptive_pitchfork_tracking();
      }

      // Now get the refined guess for the eigenvector
      for (unsigned c = 0; c < n_copies; c++)
      {
        Copy_of_problem_pt[c]->get_dofs(eigenfunction[c]);
      }
    }

    // Reactivate the tracking
    switch (bifurcation_type)
    {
        // Fold tracking
      case 1:
        this->activate_fold_tracking(parameter_pt);
        break;

        // Pitchfork
      case 2:
        this->activate_pitchfork_tracking(parameter_pt, eigenfunction[0]);
        // reset the slack parameter
        this->dof(this->ndof() - 1) = sigma;
        break;

        // Hopf
      case 3:
        this->activate_hopf_tracking(
          parameter_pt, omega, eigenfunction[0], eigenfunction[1]);
        break;

      default:
        std::ostringstream error_stream;
        error_stream << "Bifurcation type " << bifurcation_type
                     << " not known\n"
                     << "1: Fold, 2: Pitchfork, 3: Hopf\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //====================================================================
  /// A function that is used to document the errors when
  /// adapting a bifurcation-tracking
  /// problem, which requires separate interpolation of the
  /// associated eigenfunction. The error measure is chosen to be
  /// a suitable combination of the errors in the base flow and the
  /// eigenfunction. The bifurcation type is passed as an argument
  //=====================================================================
  void Problem::bifurcation_adapt_doc_errors(const unsigned& bifurcation_type)
  {
    // Dummy arguments
    unsigned n_refined, n_unrefined;
    // Just call the bifurcation helper without actually adapting
    bifurcation_adapt_helper(n_refined, n_unrefined, bifurcation_type, false);
  }


  //========================================================================
  /// Adapt problem:
  /// Perform mesh adaptation for (all) refineable (sub)mesh(es),
  /// based on their own error estimates and the target errors specified
  /// in the mesh(es). Following mesh adaptation,
  /// update global mesh, and re-assign equation numbers.
  /// Return # of refined/unrefined elements. On return from this
  /// function, Problem can immediately be solved again.
  //======================================================================
  void Problem::adapt(unsigned& n_refined, unsigned& n_unrefined)
  {
    double t_start_total = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start_total = TimingHelpers::timer();
    }

    // Get the bifurcation type
    int bifurcation_type = this->Assembly_handler_pt->bifurcation_type();

    bool continuation_problem = false;

    // If we have continuation data then we need to project that across to the
    // new mesh
    if (!Use_continuation_timestepper)
    {
      if (Dof_derivative.size() != 0)
      {
        continuation_problem = true;
      }
    }

    // If we are tracking a bifurcation then call the bifurcation adapt function
    if (bifurcation_type != 0)
    {
      this->bifurcation_adapt_helper(n_refined, n_unrefined, bifurcation_type);
      // Return immediately
      return;
    }

    if (continuation_problem)
    {
      // Create a copy of the problem
      Copy_of_problem_pt.resize(2);
      // If we don't already have a copy
      for (unsigned c = 0; c < 2; c++)
      {
        if (Copy_of_problem_pt[c] == 0)
        {
          // Create the copy
          Copy_of_problem_pt[c] = this->make_copy();

          // Refine the copy to the same level as the current problem
          // Must call actions before adapt
          Copy_of_problem_pt[c]->actions_before_adapt();

          // Find number of submeshes
          const unsigned N_mesh = Copy_of_problem_pt[c]->nsub_mesh();

          // If there is only one mesh
          if (N_mesh == 0)
          {
            // Can we refine the mesh
            if (TreeBasedRefineableMeshBase* mmesh_pt =
                  dynamic_cast<TreeBasedRefineableMeshBase*>(
                    Copy_of_problem_pt[c]->mesh_pt(0)))
            {
              // Is the adapt flag set
              if (mmesh_pt->is_adaptation_enabled())
              {
                // Now get the original problem's mesh if it's refineable
                if (TreeBasedRefineableMeshBase* original_mesh_pt =
                      dynamic_cast<TreeBasedRefineableMeshBase*>(
                        this->mesh_pt(0)))
                {
                  if (dynamic_cast<SolidMesh*>(original_mesh_pt) != 0)
                  {
                    oomph_info
                      << "Info/Warning: Adaptive Continuation is broken in "
                      << "SolidElement" << std::endl;
                  }
                  mmesh_pt->refine_base_mesh_as_in_reference_mesh(
                    original_mesh_pt);
                }
                else
                {
                  oomph_info << "Info/Warning: Mesh in orginal problem is not "
                                "refineable."
                             << std::endl;
                }
              }
              else
              {
                oomph_info
                  << "Info/Warning: Mesh adaptation is disabled in copy."
                  << std::endl;
              }
            }
            else if (TriangleMeshBase* tmesh_pt =
                       dynamic_cast<TriangleMeshBase*>(
                         Copy_of_problem_pt[c]->mesh_pt(0)))
            {
              if (TriangleMeshBase* original_mesh_pt =
                    dynamic_cast<TriangleMeshBase*>(this->mesh_pt(0)))
              {
                if (dynamic_cast<SolidMesh*>(original_mesh_pt) != 0)
                {
                  oomph_info
                    << "Info/Warning: Adaptive Continuation is broken in "
                    << "SolidElement" << std::endl;
                }

                // Remesh using the triangulateIO of the base mesh
                // Done via a file, so a bit hacky but this will be
                // superseded very soon
                std::ofstream tri_dump("triangle_mesh.dmp");
                original_mesh_pt->dump_triangulateio(tri_dump);
                tri_dump.close();
                std::ifstream tri_read("triangle_mesh.dmp");
                tmesh_pt->remesh_from_triangulateio(tri_read);
                tri_read.close();


                // Set the nodes to be at the same positions
                // as the original just in case the
                // triangulatio is out of sync with the real data
                const unsigned n_node = original_mesh_pt->nnode();
                for (unsigned n = 0; n < n_node; ++n)
                {
                  Node* const nod_pt = original_mesh_pt->node_pt(n);
                  Node* const new_node_pt = tmesh_pt->node_pt(n);
                  unsigned n_dim = nod_pt->ndim();
                  for (unsigned i = 0; i < n_dim; ++i)
                  {
                    new_node_pt->x(i) = nod_pt->x(i);
                  }
                }
              }
              else
              {
                oomph_info
                  << "Info/warning: Original Mesh is not TriangleBased\n"
                  << "... but the copy is!" << std::endl;
              }
            }
            else
            {
              oomph_info << "Info/Warning: Mesh cannot be adapted in copy."
                         << std::endl;
            }
          } // End of single mesh case
          // Otherwise loop over the submeshes
          else
          {
            for (unsigned m = 0; m < N_mesh; m++)
            {
              // Can we refine the submesh
              if (TreeBasedRefineableMeshBase* mmesh_pt =
                    dynamic_cast<TreeBasedRefineableMeshBase*>(
                      Copy_of_problem_pt[c]->mesh_pt(m)))
              {
                // Is the adapt flag set
                if (mmesh_pt->is_adaptation_enabled())
                {
                  // Now get the original problem's mesh
                  if (TreeBasedRefineableMeshBase* original_mesh_pt =
                        dynamic_cast<TreeBasedRefineableMeshBase*>(
                          this->mesh_pt(m)))
                  {
                    if (dynamic_cast<SolidMesh*>(original_mesh_pt) != 0)
                    {
                      oomph_info
                        << "Info/Warning: Adaptive Continuation is broken in "
                        << "SolidElement" << std::endl;
                    }

                    mmesh_pt->refine_base_mesh_as_in_reference_mesh(
                      original_mesh_pt);
                  }
                  else
                  {
                    oomph_info << "Info/Warning: Mesh in orginal problem is "
                                  "not refineable."
                               << std::endl;
                  }
                }
                else
                {
                  oomph_info
                    << "Info/Warning: Mesh adaptation is disabled in copy."
                    << std::endl;
                }
              }
              else if (TriangleMeshBase* tmesh_pt =
                         dynamic_cast<TriangleMeshBase*>(
                           Copy_of_problem_pt[c]->mesh_pt(m)))
              {
                if (TriangleMeshBase* original_mesh_pt =
                      dynamic_cast<TriangleMeshBase*>(this->mesh_pt(m)))
                {
                  if (dynamic_cast<SolidMesh*>(original_mesh_pt) != 0)
                  {
                    oomph_info
                      << "Info/Warning: Adaptive Continuation is broken in "
                      << "SolidElement" << std::endl;
                  }

                  // Remesh using the triangulateIO of the base mesh
                  // Done via a file, so a bit hacky but this will be
                  // superseded very soon
                  std::ofstream tri_dump("triangle_mesh.dmp");
                  original_mesh_pt->dump_triangulateio(tri_dump);
                  tri_dump.close();
                  std::ifstream tri_read("triangle_mesh.dmp");
                  tmesh_pt->remesh_from_triangulateio(tri_read);
                  tri_read.close();

                  // Set the nodes to be at the same positions
                  // as the original just in case the
                  // triangulatio is out of sync with the real data
                  const unsigned n_node = original_mesh_pt->nnode();
                  for (unsigned n = 0; n < n_node; ++n)
                  {
                    Node* const nod_pt = original_mesh_pt->node_pt(n);
                    Node* const new_node_pt = tmesh_pt->node_pt(n);
                    unsigned n_dim = nod_pt->ndim();
                    for (unsigned i = 0; i < n_dim; ++i)
                    {
                      new_node_pt->x(i) = nod_pt->x(i);
                    }
                  }
                }
                else
                {
                  oomph_info
                    << "Info/warning: Original Mesh is not TriangleBased\n"
                    << "... but the copy is!" << std::endl;
                }
              }
              else
              {
                oomph_info << "Info/Warning: Mesh cannot be adapted in copy."
                           << std::endl;
              }
            }


            // Must call actions after adapt
            Copy_of_problem_pt[c]->actions_after_adapt();

            // rebuild the global mesh in the copy
            Copy_of_problem_pt[c]->rebuild_global_mesh();

          } // End of multiple mesh case

          // Must call actions after adapt
          Copy_of_problem_pt[c]->actions_after_adapt();

          // Assign the equation numbers to the copy (quietly)
          (void)Copy_of_problem_pt[c]->assign_eqn_numbers();
        }

        // Check that the dofs match for each copy
#ifdef PARANOID
        // If the problems don't match then complain
        if (Copy_of_problem_pt[c]->ndof() != this->ndof())
        {
          std::ostringstream error_stream;
          error_stream << "Number of unknowns in the problem copy " << c << " "
                       << "not equal to number in the original:\n"
                       << this->ndof() << " (original) "
                       << Copy_of_problem_pt[c]->ndof() << " (copy)\n";

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }

      // Need to set the Dof derivatives to the copied problem
      // Assign the eigenfunction(s) to the copied problems
      unsigned ndof_local = Dof_distribution_pt->nrow_local();
      for (unsigned i = 0; i < ndof_local; i++)
      {
        Copy_of_problem_pt[0]->dof(i) = this->dof_derivative(i);
        Copy_of_problem_pt[1]->dof(i) = this->dof_current(i);
      }
      // Set all pinned values to zero
      Copy_of_problem_pt[0]->set_pinned_values_to_zero();
      // Don't need to for the current dofs that are actuall the dofs

      // Now adapt
      Vector<Vector<double>> base_error;
      this->get_all_error_estimates(base_error);
      this->adapt_based_on_error_estimates(n_refined, n_unrefined, base_error);
      Copy_of_problem_pt[0]->adapt_based_on_error_estimates(
        n_refined, n_unrefined, base_error);
      Copy_of_problem_pt[1]->adapt_based_on_error_estimates(
        n_refined, n_unrefined, base_error);

      // Now sort out the Dof pointer
      ndof_local = Dof_distribution_pt->nrow_local();
      if (Dof_derivative.size() != ndof_local)
      {
        Dof_derivative.resize(ndof_local, 0.0);
      }
      if (Dof_current.size() != ndof_local)
      {
        Dof_current.resize(ndof_local, 0.0);
      }
      for (unsigned i = 0; i < ndof_local; i++)
      {
        Dof_derivative[i] = Copy_of_problem_pt[0]->dof(i);
        Dof_current[i] = Copy_of_problem_pt[1]->dof(i);
      }
      // Return immediately
      return;
    }

    oomph_info << std::endl << std::endl;
    oomph_info << "Adapting problem:" << std::endl;
    oomph_info << "=================" << std::endl;

    double t_start = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start = TimingHelpers::timer();
    }

    // Call the actions before adaptation
    actions_before_adapt();

    double t_end = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for actions before adapt: " << t_end - t_start
                 << std::endl;
      t_start = TimingHelpers::timer();
    }

    // Initialise counters
    n_refined = 0;
    n_unrefined = 0;

    // Number of submeshes?
    unsigned Nmesh = nsub_mesh();

    // Single mesh:
    //------------
    if (Nmesh == 0)
    {
      // Refine single mesh if possible
      if (RefineableMeshBase* mmesh_pt =
            dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
      {
        if (mmesh_pt->is_adaptation_enabled())
        {
          double t_start = TimingHelpers::timer();

          // Get pointer to error estimator
          ErrorEstimator* error_estimator_pt =
            mmesh_pt->spatial_error_estimator_pt();

#ifdef PARANOID
          if (error_estimator_pt == 0)
          {
            throw OomphLibError("Error estimator hasn't been set yet",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif

          // Get error for all elements
          Vector<double> elemental_error(mmesh_pt->nelement());

          if (mmesh_pt->doc_info_pt() == 0)
          {
            error_estimator_pt->get_element_errors(mesh_pt(0), elemental_error);
          }
          else
          {
            error_estimator_pt->get_element_errors(
              mesh_pt(0), elemental_error, *mmesh_pt->doc_info_pt());
          }

          // Store max./min actual error
          mmesh_pt->max_error() = std::fabs(*std::max_element(
            elemental_error.begin(), elemental_error.end(), AbsCmp<double>()));

          mmesh_pt->min_error() = std::fabs(*std::min_element(
            elemental_error.begin(), elemental_error.end(), AbsCmp<double>()));

          oomph_info << "\n Max/min error: " << mmesh_pt->max_error() << " "
                     << mmesh_pt->min_error() << std::endl
                     << std::endl;


          if (Global_timings::Doc_comprehensive_timings)
          {
            t_end = TimingHelpers::timer();
            oomph_info << "Time for error estimation: " << t_end - t_start
                       << std::endl;
            t_start = TimingHelpers::timer();
          }

          // Adapt mesh
          mmesh_pt->adapt(elemental_error);

          // Add to counters
          n_refined += mmesh_pt->nrefined();
          n_unrefined += mmesh_pt->nunrefined();

          if (Global_timings::Doc_comprehensive_timings)
          {
            t_end = TimingHelpers::timer();
            oomph_info << "Time for complete mesh adaptation "
                       << "(but excluding comp of error estimate): "
                       << t_end - t_start << std::endl;
            t_start = TimingHelpers::timer();
          }
        }
        else
        {
          oomph_info << "Info/Warning: Mesh adaptation is disabled."
                     << std::endl;
        }
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be adapted" << std::endl;
      }
    }
    // Multiple submeshes
    //------------------
    else
    {
      // Loop over submeshes
      for (unsigned imesh = 0; imesh < Nmesh; imesh++)
      {
        // Refine single mesh uniformly if possible
        if (RefineableMeshBase* mmesh_pt =
              dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
        {
          double t_start = TimingHelpers::timer();

          // Get pointer to error estimator
          ErrorEstimator* error_estimator_pt =
            mmesh_pt->spatial_error_estimator_pt();

#ifdef PARANOID
          if (error_estimator_pt == 0)
          {
            throw OomphLibError("Error estimator hasn't been set yet",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif

          if (mmesh_pt->is_adaptation_enabled())
          {
            // Get error for all elements
            Vector<double> elemental_error(mmesh_pt->nelement());
            if (mmesh_pt->doc_info_pt() == 0)
            {
              error_estimator_pt->get_element_errors(mesh_pt(imesh),
                                                     elemental_error);
            }
            else
            {
              error_estimator_pt->get_element_errors(
                mesh_pt(imesh), elemental_error, *mmesh_pt->doc_info_pt());
            }

            // Store max./min error if the mesh has any elements
            if (mesh_pt(imesh)->nelement() > 0)
            {
              mmesh_pt->max_error() =
                std::fabs(*std::max_element(elemental_error.begin(),
                                            elemental_error.end(),
                                            AbsCmp<double>()));

              mmesh_pt->min_error() =
                std::fabs(*std::min_element(elemental_error.begin(),
                                            elemental_error.end(),
                                            AbsCmp<double>()));
            }

            oomph_info << "\n Max/min error: " << mmesh_pt->max_error() << " "
                       << mmesh_pt->min_error() << std::endl;


            if (Global_timings::Doc_comprehensive_timings)
            {
              t_end = TimingHelpers::timer();
              oomph_info << "Time for error estimation: " << t_end - t_start
                         << std::endl;
              t_start = TimingHelpers::timer();
            }

            // Adapt mesh
            mmesh_pt->adapt(elemental_error);

            // Add to counters
            n_refined += mmesh_pt->nrefined();
            n_unrefined += mmesh_pt->nunrefined();


            if (Global_timings::Doc_comprehensive_timings)
            {
              t_end = TimingHelpers::timer();
              oomph_info << "Time for complete mesh adaptation "
                         << "(but excluding comp of error estimate): "
                         << t_end - t_start << std::endl;
              t_start = TimingHelpers::timer();
            }
          }
          else
          {
            oomph_info << "Info/Warning: Mesh adaptation is disabled."
                       << std::endl;
          }
        }
        else
        {
          oomph_info << "Info/Warning: Mesh cannot be adapted." << std::endl;
        }

      } // End of loop over submeshes

      // Rebuild the global mesh
      rebuild_global_mesh();
    }


    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Total time for actual adaptation "
                 << "(all meshes; incl error estimates): " << t_end - t_start
                 << std::endl;
      t_start = TimingHelpers::timer();
    }

    // Any actions after adapt
    actions_after_adapt();


    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for actions after adapt: " << t_end - t_start
                 << std::endl;
      t_start = TimingHelpers::timer();

      oomph_info << "About to start re-assigning eqn numbers "
                 << "with Problem::assign_eqn_numbers() at end of "
                 << "Problem::adapt().\n";
    }

    // Attach the boundary conditions to the mesh
    oomph_info << "\nNumber of equations: " << assign_eqn_numbers() << std::endl
               << std::endl;


    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for re-assigning eqn numbers with "
                 << "Problem::assign_eqn_numbers() at end of Problem::adapt(): "
                 << t_end - t_start << std::endl;
      oomph_info << "Total time for adapt: " << t_end - t_start_total
                 << std::endl;
    }
  }

  //========================================================================
  /// p-adapt problem:
  /// Perform mesh adaptation for (all) refineable (sub)mesh(es),
  /// based on their own error estimates and the target errors specified
  /// in the mesh(es). Following mesh adaptation,
  /// update global mesh, and re-assign equation numbers.
  /// Return # of refined/unrefined elements. On return from this
  /// function, Problem can immediately be solved again.
  //======================================================================
  void Problem::p_adapt(unsigned& n_refined, unsigned& n_unrefined)
  {
    double t_start_total = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start_total = TimingHelpers::timer();
    }

    // Get the bifurcation type
    int bifurcation_type = this->Assembly_handler_pt->bifurcation_type();

    // If we are tracking a bifurcation then call the bifurcation adapt function
    if (bifurcation_type != 0)
    {
      this->bifurcation_adapt_helper(n_refined, n_unrefined, bifurcation_type);
      // Return immediately
      return;
    }

    oomph_info << std::endl << std::endl;
    oomph_info << "p-adapting problem:" << std::endl;
    oomph_info << "===================" << std::endl;

    double t_start = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start = TimingHelpers::timer();
    }

    // Call the actions before adaptation
    actions_before_adapt();

    double t_end = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for actions before adapt: " << t_end - t_start
                 << std::endl;
      t_start = TimingHelpers::timer();
    }

    // Initialise counters
    n_refined = 0;
    n_unrefined = 0;

    // Number of submeshes?
    unsigned Nmesh = nsub_mesh();

    // Single mesh:
    //------------
    if (Nmesh == 0)
    {
      // Refine single mesh if possible
      if (RefineableMeshBase* mmesh_pt =
            dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
      {
        if (mmesh_pt->is_p_adaptation_enabled())
        {
          double t_start = TimingHelpers::timer();

          // Get pointer to error estimator
          ErrorEstimator* error_estimator_pt =
            mmesh_pt->spatial_error_estimator_pt();

#ifdef PARANOID
          if (error_estimator_pt == 0)
          {
            throw OomphLibError("Error estimator hasn't been set yet",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif

          // Get error for all elements
          Vector<double> elemental_error(mmesh_pt->nelement());

          if (mmesh_pt->doc_info_pt() == 0)
          {
            error_estimator_pt->get_element_errors(mesh_pt(0), elemental_error);
          }
          else
          {
            error_estimator_pt->get_element_errors(
              mesh_pt(0), elemental_error, *mmesh_pt->doc_info_pt());
          }

          // Store max./min actual error
          mmesh_pt->max_error() = std::fabs(*std::max_element(
            elemental_error.begin(), elemental_error.end(), AbsCmp<double>()));

          mmesh_pt->min_error() = std::fabs(*std::min_element(
            elemental_error.begin(), elemental_error.end(), AbsCmp<double>()));

          oomph_info << "\n Max/min error: " << mmesh_pt->max_error() << " "
                     << mmesh_pt->min_error() << std::endl
                     << std::endl;


          if (Global_timings::Doc_comprehensive_timings)
          {
            t_end = TimingHelpers::timer();
            oomph_info << "Time for error estimation: " << t_end - t_start
                       << std::endl;
            t_start = TimingHelpers::timer();
          }

          // Adapt mesh
          mmesh_pt->p_adapt(elemental_error);

          // Add to counters
          n_refined += mmesh_pt->nrefined();
          n_unrefined += mmesh_pt->nunrefined();

          if (Global_timings::Doc_comprehensive_timings)
          {
            t_end = TimingHelpers::timer();
            oomph_info << "Time for complete mesh adaptation "
                       << "(but excluding comp of error estimate): "
                       << t_end - t_start << std::endl;
            t_start = TimingHelpers::timer();
          }
        }
        else
        {
          oomph_info << "Info/Warning: Mesh adaptation is disabled."
                     << std::endl;
        }
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be adapted" << std::endl;
      }
    }
    // Multiple submeshes
    //------------------
    else
    {
      // Loop over submeshes
      for (unsigned imesh = 0; imesh < Nmesh; imesh++)
      {
        // Refine single mesh uniformly if possible
        if (RefineableMeshBase* mmesh_pt =
              dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
        {
          double t_start = TimingHelpers::timer();

          // Get pointer to error estimator
          ErrorEstimator* error_estimator_pt =
            mmesh_pt->spatial_error_estimator_pt();

#ifdef PARANOID
          if (error_estimator_pt == 0)
          {
            throw OomphLibError("Error estimator hasn't been set yet",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif

          if (mmesh_pt->is_p_adaptation_enabled())
          {
            // Get error for all elements
            Vector<double> elemental_error(mmesh_pt->nelement());
            if (mmesh_pt->doc_info_pt() == 0)
            {
              error_estimator_pt->get_element_errors(mesh_pt(imesh),
                                                     elemental_error);
            }
            else
            {
              error_estimator_pt->get_element_errors(
                mesh_pt(imesh), elemental_error, *mmesh_pt->doc_info_pt());
            }

            // Store max./min error if the mesh has any elements
            if (mesh_pt(imesh)->nelement() > 0)
            {
              mmesh_pt->max_error() =
                std::fabs(*std::max_element(elemental_error.begin(),
                                            elemental_error.end(),
                                            AbsCmp<double>()));

              mmesh_pt->min_error() =
                std::fabs(*std::min_element(elemental_error.begin(),
                                            elemental_error.end(),
                                            AbsCmp<double>()));
            }

            oomph_info << "\n Max/min error: " << mmesh_pt->max_error() << " "
                       << mmesh_pt->min_error() << std::endl;


            if (Global_timings::Doc_comprehensive_timings)
            {
              t_end = TimingHelpers::timer();
              oomph_info << "Time for error estimation: " << t_end - t_start
                         << std::endl;
              t_start = TimingHelpers::timer();
            }

            // Adapt mesh
            mmesh_pt->p_adapt(elemental_error);

            // Add to counters
            n_refined += mmesh_pt->nrefined();
            n_unrefined += mmesh_pt->nunrefined();


            if (Global_timings::Doc_comprehensive_timings)
            {
              t_end = TimingHelpers::timer();
              oomph_info << "Time for complete mesh adaptation "
                         << "(but excluding comp of error estimate): "
                         << t_end - t_start << std::endl;
              t_start = TimingHelpers::timer();
            }
          }
          else
          {
            oomph_info << "Info/Warning: Mesh adaptation is disabled."
                       << std::endl;
          }
        }
        else
        {
          oomph_info << "Info/Warning: Mesh cannot be adapted." << std::endl;
        }

      } // End of loop over submeshes

      // Rebuild the global mesh
      rebuild_global_mesh();
    }


    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Total time for actual adaptation "
                 << "(all meshes; incl error estimates): " << t_end - t_start
                 << std::endl;
      t_start = TimingHelpers::timer();
    }

    // Any actions after adapt
    actions_after_adapt();


    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for actions after adapt: " << t_end - t_start
                 << std::endl;
      t_start = TimingHelpers::timer();

      oomph_info << "About to start re-assigning eqn numbers "
                 << "with Problem::assign_eqn_numbers() at end of "
                 << "Problem::adapt().\n";
    }

    // Attach the boundary conditions to the mesh
    oomph_info << "\nNumber of equations: " << assign_eqn_numbers() << std::endl
               << std::endl;


    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for re-assigning eqn numbers with "
                 << "Problem::assign_eqn_numbers() at end of Problem::adapt(): "
                 << t_end - t_start << std::endl;
      oomph_info << "Total time for adapt: " << t_end - t_start_total
                 << std::endl;
    }
  }

  //========================================================================
  /// Perform mesh adaptation for (all) refineable (sub)mesh(es),
  /// based on the error estimates in elemental_error
  /// and the target errors specified
  /// in the mesh(es). Following mesh adaptation,
  /// update global mesh, and re-assign equation numbers.
  /// Return # of refined/unrefined elements. On return from this
  /// function, Problem can immediately be solved again.
  //========================================================================
  void Problem::adapt_based_on_error_estimates(
    unsigned& n_refined,
    unsigned& n_unrefined,
    Vector<Vector<double>>& elemental_error)
  {
    oomph_info << std::endl << std::endl;
    oomph_info << "Adapting problem:" << std::endl;
    oomph_info << "=================" << std::endl;

    // Call the actions before adaptation
    actions_before_adapt();

    // Initialise counters
    n_refined = 0;
    n_unrefined = 0;

    // Number of submeshes?
    unsigned Nmesh = nsub_mesh();

    // Single mesh:
    //------------
    if (Nmesh == 0)
    {
      // Refine single mesh uniformly if possible
      if (RefineableMeshBase* mmesh_pt =
            dynamic_cast<RefineableMeshBase*>(Problem::mesh_pt(0)))
      {
        if (mmesh_pt->is_adaptation_enabled())
        {
          // Adapt mesh
          mmesh_pt->adapt(elemental_error[0]);

          // Add to counters
          n_refined += mmesh_pt->nrefined();
          n_unrefined += mmesh_pt->nunrefined();
        }
        else
        {
          oomph_info << "Info/Warning: Mesh adaptation is disabled."
                     << std::endl;
        }
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be adapted" << std::endl;
      }
    }

    // Multiple submeshes
    //------------------
    else
    {
      // Loop over submeshes
      for (unsigned imesh = 0; imesh < Nmesh; imesh++)
      {
        // Refine single mesh uniformly if possible
        if (RefineableMeshBase* mmesh_pt =
              dynamic_cast<RefineableMeshBase*>(Problem::mesh_pt(imesh)))
        {
          if (mmesh_pt->is_adaptation_enabled())
          {
            // Adapt mesh
            mmesh_pt->adapt(elemental_error[imesh]);

            // Add to counters
            n_refined += mmesh_pt->nrefined();
            n_unrefined += mmesh_pt->nunrefined();
          }
          else
          {
            oomph_info << "Info/Warning: Mesh adaptation is disabled."
                       << std::endl;
          }
        }
        else
        {
          oomph_info << "Info/Warning: Mesh cannot be adapted." << std::endl;
        }

      } // End of loop over submeshes

      // Rebuild the global mesh
      rebuild_global_mesh();
    }

    // Any actions after adapt
    actions_after_adapt();

    // Attach the boundary conditions to the mesh
    oomph_info << "\nNumber of equations: " << assign_eqn_numbers() << std::endl
               << std::endl;
  }


  //========================================================================
  /// Return the error estimates computed by (all) refineable
  /// (sub)mesh(es) in the elemental_error structure, which consists of
  /// a vector of elemental errors for each (sub)mesh.
  //========================================================================
  void Problem::get_all_error_estimates(Vector<Vector<double>>& elemental_error)
  {
    // Number of submeshes?
    const unsigned Nmesh = nsub_mesh();

    // Single mesh:
    //------------
    if (Nmesh == 0)
    {
      // There is only one mesh
      elemental_error.resize(1);
      // Refine single mesh uniformly if possible
      if (RefineableMeshBase* mmesh_pt =
            dynamic_cast<RefineableMeshBase*>(Problem::mesh_pt(0)))
      {
        // If we can adapt the mesh
        if (mmesh_pt->is_adaptation_enabled())
        {
          // Get pointer to error estimator
          ErrorEstimator* error_estimator_pt =
            mmesh_pt->spatial_error_estimator_pt();

#ifdef PARANOID
          if (error_estimator_pt == 0)
          {
            throw OomphLibError("Error estimator hasn't been set yet",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif

          // Get error for all elements
          elemental_error[0].resize(mmesh_pt->nelement());
          // Are we documenting the errors or not
          if (mmesh_pt->doc_info_pt() == 0)
          {
            error_estimator_pt->get_element_errors(Problem::mesh_pt(0),
                                                   elemental_error[0]);
          }
          else
          {
            error_estimator_pt->get_element_errors(Problem::mesh_pt(0),
                                                   elemental_error[0],
                                                   *mmesh_pt->doc_info_pt());
          }

          // Store max./min actual error
          mmesh_pt->max_error() =
            std::fabs(*std::max_element(elemental_error[0].begin(),
                                        elemental_error[0].end(),
                                        AbsCmp<double>()));

          mmesh_pt->min_error() =
            std::fabs(*std::min_element(elemental_error[0].begin(),
                                        elemental_error[0].end(),
                                        AbsCmp<double>()));

          oomph_info << "\n Max/min error: " << mmesh_pt->max_error() << " "
                     << mmesh_pt->min_error() << std::endl;
        }
        else
        {
          oomph_info << "Info/Warning: Mesh adaptation is disabled."
                     << std::endl;
        }
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be adapted" << std::endl;
      }
    }

    // Multiple submeshes
    //------------------
    else
    {
      // Resize to the number of submeshes
      elemental_error.resize(Nmesh);

      // Loop over submeshes
      for (unsigned imesh = 0; imesh < Nmesh; imesh++)
      {
        // Refine single mesh uniformly if possible
        if (RefineableMeshBase* mmesh_pt =
              dynamic_cast<RefineableMeshBase*>(Problem::mesh_pt(imesh)))
        {
          // Get pointer to error estimator
          ErrorEstimator* error_estimator_pt =
            mmesh_pt->spatial_error_estimator_pt();

#ifdef PARANOID
          if (error_estimator_pt == 0)
          {
            throw OomphLibError("Error estimator hasn't been set yet",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif
          // If we can adapt the mesh
          if (mmesh_pt->is_adaptation_enabled())
          {
            // Get error for all elements
            elemental_error[imesh].resize(mmesh_pt->nelement());
            if (mmesh_pt->doc_info_pt() == 0)
            {
              error_estimator_pt->get_element_errors(Problem::mesh_pt(imesh),
                                                     elemental_error[imesh]);
            }
            else
            {
              error_estimator_pt->get_element_errors(Problem::mesh_pt(imesh),
                                                     elemental_error[imesh],
                                                     *mmesh_pt->doc_info_pt());
            }

            // Store max./min error
            mmesh_pt->max_error() =
              std::fabs(*std::max_element(elemental_error[imesh].begin(),
                                          elemental_error[imesh].end(),
                                          AbsCmp<double>()));

            mmesh_pt->min_error() =
              std::fabs(*std::min_element(elemental_error[imesh].begin(),
                                          elemental_error[imesh].end(),
                                          AbsCmp<double>()));

            oomph_info << "\n Max/min error: " << mmesh_pt->max_error() << " "
                       << mmesh_pt->min_error() << std::endl;
          }
          else
          {
            oomph_info << "Info/Warning: Mesh adaptation is disabled."
                       << std::endl;
          }
        }
        else
        {
          oomph_info << "Info/Warning: Mesh cannot be adapted." << std::endl;
        }

      } // End of loop over submeshes
    }
  }

  //========================================================================
  /// \short Get max and min error for all elements in submeshes
  //========================================================================
  void Problem::doc_errors(DocInfo& doc_info)
  {
    // Get the bifurcation type
    int bifurcation_type = this->Assembly_handler_pt->bifurcation_type();
    // If we are tracking a bifurcation then call the bifurcation adapt function
    if (bifurcation_type != 0)
    {
      this->bifurcation_adapt_doc_errors(bifurcation_type);
      // Return immediately
      return;
    }

    // Number of submeshes?
    unsigned Nmesh = nsub_mesh();

    // Single mesh:
    //------------
    if (Nmesh == 0)
    {
      // Is the single mesh refineable?
      if (RefineableMeshBase* mmesh_pt =
            dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
      {
        // Get pointer to error estimator
        ErrorEstimator* error_estimator_pt =
          mmesh_pt->spatial_error_estimator_pt();

#ifdef PARANOID
        if (error_estimator_pt == 0)
        {
          throw OomphLibError("Error estimator hasn't been set yet",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Get error for all elements
        Vector<double> elemental_error(mmesh_pt->nelement());
        if (!doc_info.is_doc_enabled())
        {
          error_estimator_pt->get_element_errors(mesh_pt(0), elemental_error);
        }
        else
        {
          error_estimator_pt->get_element_errors(
            mesh_pt(0), elemental_error, doc_info);
        }

        // Store max./min actual error
        mmesh_pt->max_error() = std::fabs(*std::max_element(
          elemental_error.begin(), elemental_error.end(), AbsCmp<double>()));

        mmesh_pt->min_error() = std::fabs(*std::min_element(
          elemental_error.begin(), elemental_error.end(), AbsCmp<double>()));

        oomph_info << "\n Max/min error: " << mmesh_pt->max_error() << " "
                   << mmesh_pt->min_error() << std::endl;
      }
    }

    // Multiple submeshes
    //------------------
    else
    {
      // Loop over submeshes
      for (unsigned imesh = 0; imesh < Nmesh; imesh++)
      {
        // Is the single mesh refineable?
        if (RefineableMeshBase* mmesh_pt =
              dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
        {
          // Get pointer to error estimator
          ErrorEstimator* error_estimator_pt =
            mmesh_pt->spatial_error_estimator_pt();

#ifdef PARANOID
          if (error_estimator_pt == 0)
          {
            throw OomphLibError("Error estimator hasn't been set yet",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif

          // Get error for all elements
          Vector<double> elemental_error(mmesh_pt->nelement());
          if (mmesh_pt->doc_info_pt() == 0)
          {
            error_estimator_pt->get_element_errors(mesh_pt(imesh),
                                                   elemental_error);
          }
          else
          {
            error_estimator_pt->get_element_errors(
              mesh_pt(imesh), elemental_error, *mmesh_pt->doc_info_pt());
          }

          // Store max./min error if the mesh has any elements
          if (mesh_pt(imesh)->nelement() > 0)
          {
            mmesh_pt->max_error() =
              std::fabs(*std::max_element(elemental_error.begin(),
                                          elemental_error.end(),
                                          AbsCmp<double>()));

            mmesh_pt->min_error() =
              std::fabs(*std::min_element(elemental_error.begin(),
                                          elemental_error.end(),
                                          AbsCmp<double>()));
          }

          oomph_info << "\n Max/min error: " << mmesh_pt->max_error() << " "
                     << mmesh_pt->min_error() << std::endl;
        }

      } // End of loop over submeshes
    }
  }

  //========================================================================
  /// Refine (one and only!) mesh by splitting the elements identified
  /// by their numbers relative to the problems' only mesh, then rebuild
  /// the problem.
  //========================================================================
  void Problem::refine_selected_elements(
    const Vector<unsigned>& elements_to_be_refined)
  {
    actions_before_adapt();

    // Number of submeshes?
    unsigned Nmesh = nsub_mesh();

    // Single mesh:
    if (Nmesh == 0)
    {
      // Refine single mesh if possible
      if (TreeBasedRefineableMeshBase* mmesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(0)))
      {
        mmesh_pt->refine_selected_elements(elements_to_be_refined);
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be refined " << std::endl;
      }
    }
    // Multiple submeshes
    else
    {
      std::ostringstream error_message;
      error_message << "Problem::refine_selected_elements(...) only works for\n"
                    << "multiple-mesh problems if you specify the mesh\n"
                    << "number in the function argument before the Vector,\n"
                    << "or a Vector of Vectors for each submesh.\n"
                    << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Any actions after the adapatation phase
    actions_after_adapt();

    // Attach the boundary conditions to the mesh
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }

  //========================================================================
  /// Refine (one and only!) mesh by splitting the elements identified
  /// by their pointers, then rebuild the problem.
  //========================================================================
  void Problem::refine_selected_elements(
    const Vector<RefineableElement*>& elements_to_be_refined_pt)
  {
    actions_before_adapt();

    // Number of submeshes?
    unsigned Nmesh = nsub_mesh();

    // Single mesh:
    if (Nmesh == 0)
    {
      // Refine single mesh if possible
      if (TreeBasedRefineableMeshBase* mmesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(0)))
      {
        mmesh_pt->refine_selected_elements(elements_to_be_refined_pt);
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be refined " << std::endl;
      }
    }
    // Multiple submeshes
    else
    {
      std::ostringstream error_message;
      error_message << "Problem::refine_selected_elements(...) only works for\n"
                    << "multiple-mesh problems if you specify the mesh\n"
                    << "number in the function argument before the Vector,\n"
                    << "or a Vector of Vectors for each submesh.\n"
                    << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Any actions after the adapatation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }

  //========================================================================
  /// Refine specified submesh by splitting the elements identified
  /// by their numbers relative to the specified mesh, then rebuild the problem.
  //========================================================================
  void Problem::refine_selected_elements(
    const unsigned& i_mesh, const Vector<unsigned>& elements_to_be_refined)
  {
    actions_before_adapt();

    // Number of submeshes?
    unsigned n_mesh = nsub_mesh();

    if (i_mesh >= n_mesh)
    {
      std::ostringstream error_message;
      error_message << "Problem only has " << n_mesh
                    << " submeshes. Cannot refine submesh " << i_mesh
                    << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Refine single mesh if possible
    if (TreeBasedRefineableMeshBase* mmesh_pt =
          dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(i_mesh)))
    {
      mmesh_pt->refine_selected_elements(elements_to_be_refined);
    }
    else
    {
      oomph_info << "Info/Warning: Mesh cannot be refined " << std::endl;
    }

    if (n_mesh > 1)
    {
      // Rebuild the global mesh
      rebuild_global_mesh();
    }

    // Any actions after the adapatation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }


  //========================================================================
  /// Refine specified submesh by splitting the elements identified
  /// by their pointers, then rebuild the problem.
  //========================================================================
  void Problem::refine_selected_elements(
    const unsigned& i_mesh,
    const Vector<RefineableElement*>& elements_to_be_refined_pt)
  {
    actions_before_adapt();

    // Number of submeshes?
    unsigned n_mesh = nsub_mesh();

    if (i_mesh >= n_mesh)
    {
      std::ostringstream error_message;
      error_message << "Problem only has " << n_mesh
                    << " submeshes. Cannot refine submesh " << i_mesh
                    << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Refine single mesh if possible
    if (TreeBasedRefineableMeshBase* mmesh_pt =
          dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(i_mesh)))
    {
      mmesh_pt->refine_selected_elements(elements_to_be_refined_pt);
    }
    else
    {
      oomph_info << "Info/Warning: Mesh cannot be refined " << std::endl;
    }

    if (n_mesh > 1)
    {
      // Rebuild the global mesh
      rebuild_global_mesh();
    }

    // Any actions after the adapatation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }

  //========================================================================
  /// Refine all submeshes by splitting the elements identified by their
  /// numbers relative to each submesh in a Vector of Vectors, then
  /// rebuild the problem.
  //========================================================================
  void Problem::refine_selected_elements(
    const Vector<Vector<unsigned>>& elements_to_be_refined)
  {
    actions_before_adapt();

    // Number of submeshes?
    unsigned n_mesh = nsub_mesh();

    // Refine all submeshes if possible
    for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
    {
      if (TreeBasedRefineableMeshBase* mmesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(i_mesh)))
      {
        mmesh_pt->refine_selected_elements(elements_to_be_refined[i_mesh]);
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be refined " << std::endl;
      }
    }

    // Rebuild the global mesh
    rebuild_global_mesh();

    // Any actions after the adapatation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }

  //========================================================================
  /// Refine all submeshes by splitting the elements identified by their
  /// pointers within each submesh in a Vector of Vectors, then
  /// rebuild the problem.
  //========================================================================
  void Problem::refine_selected_elements(
    const Vector<Vector<RefineableElement*>>& elements_to_be_refined_pt)
  {
    actions_before_adapt();

    // Number of submeshes?
    unsigned n_mesh = nsub_mesh();

    // Refine all submeshes if possible
    for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
    {
      if (TreeBasedRefineableMeshBase* mmesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(i_mesh)))
      {
        mmesh_pt->refine_selected_elements(elements_to_be_refined_pt[i_mesh]);
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be refined " << std::endl;
      }
    }

    // Rebuild the global mesh
    rebuild_global_mesh();

    // Any actions after the adapatation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }

  //========================================================================
  /// p-refine (one and only!) mesh by refining the elements identified
  /// by their numbers relative to the problems' only mesh, then rebuild
  /// the problem.
  //========================================================================
  void Problem::p_refine_selected_elements(
    const Vector<unsigned>& elements_to_be_refined)
  {
    actions_before_adapt();

    // Number of submeshes?
    unsigned Nmesh = nsub_mesh();

    // Single mesh:
    if (Nmesh == 0)
    {
      // Refine single mesh if possible
      if (TreeBasedRefineableMeshBase* mmesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(0)))
      {
        mmesh_pt->p_refine_selected_elements(elements_to_be_refined);
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be refined " << std::endl;
      }
    }
    // Multiple submeshes
    else
    {
      std::ostringstream error_message;
      error_message
        << "Problem::p_refine_selected_elements(...) only works for\n"
        << "multiple-mesh problems if you specify the mesh\n"
        << "number in the function argument before the Vector,\n"
        << "or a Vector of Vectors for each submesh.\n"
        << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Any actions after the adapatation phase
    actions_after_adapt();

    // Attach the boundary conditions to the mesh
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }

  //========================================================================
  /// p-refine (one and only!) mesh by refining the elements identified
  /// by their pointers, then rebuild the problem.
  //========================================================================
  void Problem::p_refine_selected_elements(
    const Vector<PRefineableElement*>& elements_to_be_refined_pt)
  {
    actions_before_adapt();

    // Number of submeshes?
    unsigned Nmesh = nsub_mesh();

    // Single mesh:
    if (Nmesh == 0)
    {
      // Refine single mesh if possible
      if (TreeBasedRefineableMeshBase* mmesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(0)))
      {
        mmesh_pt->p_refine_selected_elements(elements_to_be_refined_pt);
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be refined " << std::endl;
      }
    }
    // Multiple submeshes
    else
    {
      std::ostringstream error_message;
      error_message
        << "Problem::p_refine_selected_elements(...) only works for\n"
        << "multiple-mesh problems if you specify the mesh\n"
        << "number in the function argument before the Vector,\n"
        << "or a Vector of Vectors for each submesh.\n"
        << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Any actions after the adapatation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }

  //========================================================================
  /// p-refine specified submesh by refining the elements identified
  /// by their numbers relative to the specified mesh, then rebuild the problem.
  //========================================================================
  void Problem::p_refine_selected_elements(
    const unsigned& i_mesh, const Vector<unsigned>& elements_to_be_refined)
  {
    OomphLibWarning(
      "p-refinement for multiple submeshes has not yet been tested.",
      "Problem::p_refine_selected_elements()",
      OOMPH_EXCEPTION_LOCATION);

    actions_before_adapt();

    // Number of submeshes?
    unsigned n_mesh = nsub_mesh();

    if (i_mesh >= n_mesh)
    {
      std::ostringstream error_message;
      error_message << "Problem only has " << n_mesh
                    << " submeshes. Cannot p-refine submesh " << i_mesh
                    << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Refine single mesh if possible
    if (TreeBasedRefineableMeshBase* mmesh_pt =
          dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(i_mesh)))
    {
      mmesh_pt->p_refine_selected_elements(elements_to_be_refined);
    }
    else
    {
      oomph_info << "Info/Warning: Mesh cannot be refined " << std::endl;
    }

    if (n_mesh > 1)
    {
      // Rebuild the global mesh
      rebuild_global_mesh();
    }

    // Any actions after the adapatation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }


  //========================================================================
  /// p-refine specified submesh by refining the elements identified
  /// by their pointers, then rebuild the problem.
  //========================================================================
  void Problem::p_refine_selected_elements(
    const unsigned& i_mesh,
    const Vector<PRefineableElement*>& elements_to_be_refined_pt)
  {
    OomphLibWarning(
      "p-refinement for multiple submeshes has not yet been tested.",
      "Problem::p_refine_selected_elements()",
      OOMPH_EXCEPTION_LOCATION);

    actions_before_adapt();

    // Number of submeshes?
    unsigned n_mesh = nsub_mesh();

    if (i_mesh >= n_mesh)
    {
      std::ostringstream error_message;
      error_message << "Problem only has " << n_mesh
                    << " submeshes. Cannot p-refine submesh " << i_mesh
                    << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Refine single mesh if possible
    if (TreeBasedRefineableMeshBase* mmesh_pt =
          dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(i_mesh)))
    {
      mmesh_pt->p_refine_selected_elements(elements_to_be_refined_pt);
    }
    else
    {
      oomph_info << "Info/Warning: Mesh cannot be refined " << std::endl;
    }

    if (n_mesh > 1)
    {
      // Rebuild the global mesh
      rebuild_global_mesh();
    }

    // Any actions after the adapatation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }

  //========================================================================
  /// p-refine all submeshes by refining the elements identified by their
  /// numbers relative to each submesh in a Vector of Vectors, then
  /// rebuild the problem.
  //========================================================================
  void Problem::p_refine_selected_elements(
    const Vector<Vector<unsigned>>& elements_to_be_refined)
  {
    OomphLibWarning(
      "p-refinement for multiple submeshes has not yet been tested.",
      "Problem::p_refine_selected_elements()",
      OOMPH_EXCEPTION_LOCATION);

    actions_before_adapt();

    // Number of submeshes?
    unsigned n_mesh = nsub_mesh();

    // Refine all submeshes if possible
    for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
    {
      if (TreeBasedRefineableMeshBase* mmesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(i_mesh)))
      {
        mmesh_pt->p_refine_selected_elements(elements_to_be_refined[i_mesh]);
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be refined " << std::endl;
      }
    }

    // Rebuild the global mesh
    rebuild_global_mesh();

    // Any actions after the adapatation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }

  //========================================================================
  /// p-refine all submeshes by refining the elements identified by their
  /// pointers within each submesh in a Vector of Vectors, then
  /// rebuild the problem.
  //========================================================================
  void Problem::p_refine_selected_elements(
    const Vector<Vector<PRefineableElement*>>& elements_to_be_refined_pt)
  {
    OomphLibWarning(
      "p-refinement for multiple submeshes has not yet been tested.",
      "Problem::p_refine_selected_elements()",
      OOMPH_EXCEPTION_LOCATION);

    actions_before_adapt();

    // Number of submeshes?
    unsigned n_mesh = nsub_mesh();

    // Refine all submeshes if possible
    for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
    {
      if (TreeBasedRefineableMeshBase* mmesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(i_mesh)))
      {
        mmesh_pt->p_refine_selected_elements(elements_to_be_refined_pt[i_mesh]);
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be refined " << std::endl;
      }
    }

    // Rebuild the global mesh
    rebuild_global_mesh();

    // Any actions after the adapatation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }


  //========================================================================
  /// Helper function to do compund refinement of (all) refineable
  /// (sub)mesh(es) uniformly as many times as specified in vector and
  /// rebuild problem; doc refinement process. Set boolean argument
  /// to true if you want to prune immediately after refining the meshes
  /// individually.
  //========================================================================
  void Problem::refine_uniformly_aux(const Vector<unsigned>& nrefine_for_mesh,
                                     DocInfo& doc_info,
                                     const bool& prune)
  {
    double t_start = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start = TimingHelpers::timer();
    }

    actions_before_adapt();

    double t_end = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info
        << "Time for actions before adapt in Problem::refine_uniformly_aux(): "
        << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }

    // Number of submeshes?
    unsigned n_mesh = nsub_mesh();

    // Single mesh:
    if (n_mesh == 0)
    {
      // Refine single mesh uniformly if possible
      if (RefineableMeshBase* mmesh_pt =
            dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
      {
        unsigned nref = nrefine_for_mesh[0];
        for (unsigned i = 0; i < nref; i++)
        {
          mmesh_pt->refine_uniformly(doc_info);
        }
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be refined uniformly "
                   << std::endl;
      }
    }
    // Multiple submeshes
    else
    {
      // Loop over submeshes
      for (unsigned imesh = 0; imesh < n_mesh; imesh++)
      {
        // Refine i-th submesh uniformly if possible
        if (RefineableMeshBase* mmesh_pt =
              dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
        {
          unsigned nref = nrefine_for_mesh[imesh];
          for (unsigned i = 0; i < nref; i++)
          {
            mmesh_pt->refine_uniformly(doc_info);
          }
        }
        else
        {
          oomph_info << "Info/Warning: Cannot refine mesh " << imesh
                     << std::endl;
        }
      }
      // Rebuild the global mesh
      rebuild_global_mesh();
    }

    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for mesh-level mesh refinement in "
                 << "Problem::refine_uniformly_aux(): " << t_end - t_start
                 << std::endl;
      t_start = TimingHelpers::timer();
    }

    // Any actions after the adaptation phase
    actions_after_adapt();


    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info
        << "Time for actions after adapt  Problem::refine_uniformly_aux(): "
        << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }


#ifdef OOMPH_HAS_MPI

    // Prune it?
    if (prune)
    {
      // Note: This calls assign eqn numbers already...
      Bypass_increase_in_dof_check_during_pruning = true;
      prune_halo_elements_and_nodes();
      Bypass_increase_in_dof_check_during_pruning = false;

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for Problem::prune_halo_elements_and_nodes() in "
                   << "Problem::refine_uniformly_aux(): " << t_end - t_start
                   << std::endl;
      }
    }
    else
#else
    if (prune)
    {
      std::ostringstream error_message;
      error_message
        << "Requested pruning in serial build. Ignoring the request.\n";
      OomphLibWarning(error_message.str(),
                      "Problem::refine_uniformly_aux()",
                      OOMPH_EXCEPTION_LOCATION);
    }
#endif
    {
      // Do equation numbering
      oomph_info
        << "Number of equations after Problem::refine_uniformly_aux(): "
        << assign_eqn_numbers() << std::endl;

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for Problem::assign_eqn_numbers() in "
                   << "Problem::refine_uniformly_aux(): " << t_end - t_start
                   << std::endl;
      }
    }
  }


  //========================================================================
  /// Helper function to do compund p-refinement of (all) p-refineable
  /// (sub)mesh(es) uniformly as many times as specified in vector and
  /// rebuild problem; doc refinement process. Set boolean argument
  /// to true if you want to prune immediately after refining the meshes
  /// individually.
  //========================================================================
  void Problem::p_refine_uniformly_aux(const Vector<unsigned>& nrefine_for_mesh,
                                       DocInfo& doc_info,
                                       const bool& prune)
  {
    double t_start = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start = TimingHelpers::timer();
    }

    actions_before_adapt();

    double t_end = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for actions before adapt in "
                    "Problem::p_refine_uniformly_aux(): "
                 << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }

    // Number of submeshes?
    unsigned n_mesh = nsub_mesh();

    // Single mesh:
    if (n_mesh == 0)
    {
      // Refine single mesh uniformly if possible
      if (RefineableMeshBase* mmesh_pt =
            dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
      {
        unsigned nref = nrefine_for_mesh[0];
        for (unsigned i = 0; i < nref; i++)
        {
          mmesh_pt->p_refine_uniformly(doc_info);
        }
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be p-refined uniformly "
                   << std::endl;
      }
    }
    // Multiple submeshes
    else
    {
      OomphLibWarning(
        "p-refinement for multiple submeshes has not yet been tested.",
        "Problem::p_refine_uniformly_aux()",
        OOMPH_EXCEPTION_LOCATION);

      // Loop over submeshes
      for (unsigned imesh = 0; imesh < n_mesh; imesh++)
      {
        // Refine i-th submesh uniformly if possible
        if (RefineableMeshBase* mmesh_pt =
              dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
        {
          unsigned nref = nrefine_for_mesh[imesh];
          for (unsigned i = 0; i < nref; i++)
          {
            mmesh_pt->p_refine_uniformly(doc_info);
          }
        }
        else
        {
          oomph_info << "Info/Warning: Cannot p-refine mesh " << imesh
                     << std::endl;
        }
      }
      // Rebuild the global mesh
      rebuild_global_mesh();
    }

    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for mesh-level mesh refinement in "
                 << "Problem::p_refine_uniformly_aux(): " << t_end - t_start
                 << std::endl;
      t_start = TimingHelpers::timer();
    }

    // Any actions after the adaptation phase
    actions_after_adapt();


    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info
        << "Time for actions after adapt  Problem::p_refine_uniformly_aux(): "
        << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }


#ifdef OOMPH_HAS_MPI

    // Prune it?
    if (prune)
    {
      // Note: This calls assign eqn numbers already...
      Bypass_increase_in_dof_check_during_pruning = true;
      prune_halo_elements_and_nodes();
      Bypass_increase_in_dof_check_during_pruning = false;

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for Problem::prune_halo_elements_and_nodes() in "
                   << "Problem::p_refine_uniformly_aux(): " << t_end - t_start
                   << std::endl;
      }
    }
    else
#else
    if (prune)
    {
      std::ostringstream error_message;
      error_message
        << "Requested pruning in serial build. Ignoring the request.\n";
      OomphLibWarning(error_message.str(),
                      "Problem::p_refine_uniformly_aux()",
                      OOMPH_EXCEPTION_LOCATION);
    }
#endif
    {
      // Do equation numbering
      oomph_info
        << "Number of equations after Problem::p_refine_uniformly_aux(): "
        << assign_eqn_numbers() << std::endl;

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for Problem::assign_eqn_numbers() in "
                   << "Problem::p_refine_uniformly_aux(): " << t_end - t_start
                   << std::endl;
      }
    }
  }

  //========================================================================
  /// Refine submesh i_mesh uniformly and rebuild problem;
  /// doc refinement process.
  //========================================================================
  void Problem::refine_uniformly(const unsigned& i_mesh, DocInfo& doc_info)
  {
    actions_before_adapt();

#ifdef PARANOID
    // Number of submeshes?
    if (i_mesh >= nsub_mesh())
    {
      std::ostringstream error_message;
      error_message << "imesh " << i_mesh
                    << " is greater than the number of sub meshes "
                    << nsub_mesh() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Refine single mesh uniformly if possible
    if (RefineableMeshBase* mmesh_pt =
          dynamic_cast<RefineableMeshBase*>(mesh_pt(i_mesh)))
    {
      mmesh_pt->refine_uniformly(doc_info);
    }
    else
    {
      oomph_info << "Info/Warning: Mesh cannot be refined uniformly "
                 << std::endl;
    }

    // Rebuild the global mesh
    rebuild_global_mesh();

    // Any actions after the adaptation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }

  //========================================================================
  /// p-refine submesh i_mesh uniformly and rebuild problem;
  /// doc refinement process.
  //========================================================================
  void Problem::p_refine_uniformly(const unsigned& i_mesh, DocInfo& doc_info)
  {
    actions_before_adapt();

#ifdef PARANOID
    // Number of submeshes?
    if (i_mesh >= nsub_mesh())
    {
      std::ostringstream error_message;
      error_message << "imesh " << i_mesh
                    << " is greater than the number of sub meshes "
                    << nsub_mesh() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Refine single mesh uniformly if possible
    if (RefineableMeshBase* mmesh_pt =
          dynamic_cast<RefineableMeshBase*>(mesh_pt(i_mesh)))
    {
      mmesh_pt->p_refine_uniformly(doc_info);
    }
    else
    {
      oomph_info << "Info/Warning: Mesh cannot be refined uniformly "
                 << std::endl;
    }

    // Rebuild the global mesh
    rebuild_global_mesh();

    // Any actions after the adaptation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }


  //========================================================================
  /// Unrefine (all) refineable (sub)mesh(es) uniformly and rebuild problem.
  /// Return 0 for success,
  /// 1 for failure (if unrefinement has reached the coarsest permitted
  /// level)
  //========================================================================
  unsigned Problem::unrefine_uniformly()
  {
    // Call actions_before_adapt()
    actions_before_adapt();

    // Has unrefinement been successful?
    unsigned success_flag = 0;

    // Number of submeshes?
    unsigned n_mesh = nsub_mesh();

    // Single mesh:
    if (n_mesh == 0)
    {
      // Unrefine single mesh uniformly if possible
      if (RefineableMeshBase* mmesh_pt =
            dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
      {
        success_flag += mmesh_pt->unrefine_uniformly();
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be unrefined uniformly "
                   << std::endl;
      }
    }
    // Multiple submeshes
    else
    {
      // Loop over submeshes
      for (unsigned imesh = 0; imesh < n_mesh; imesh++)
      {
        // Unrefine i-th submesh uniformly if possible
        if (RefineableMeshBase* mmesh_pt =
              dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
        {
          success_flag += mmesh_pt->unrefine_uniformly();
        }
        else
        {
          oomph_info << "Info/Warning: Cannot unrefine mesh " << imesh
                     << std::endl;
        }
      }
      // Rebuild the global mesh
      rebuild_global_mesh();
    }

    // Any actions after the adaptation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << " Number of equations: " << assign_eqn_numbers() << std::endl;

    // Judge success
    if (success_flag > 0)
    {
      return 1;
    }
    else
    {
      return 0;
    }
  }

  //========================================================================
  /// Unrefine submesh i_mesh uniformly and rebuild problem.
  /// Return 0 for success,
  /// 1 for failure (if unrefinement has reached the coarsest permitted
  /// level)
  //========================================================================
  unsigned Problem::unrefine_uniformly(const unsigned& i_mesh)
  {
    actions_before_adapt();

    // Has unrefinement been successful?
    unsigned success_flag = 0;

#ifdef PARANOID
    // Number of submeshes?
    if (i_mesh >= nsub_mesh())
    {
      std::ostringstream error_message;
      error_message << "imesh " << i_mesh
                    << " is greater than the number of sub meshes "
                    << nsub_mesh() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Unrefine single mesh uniformly if possible
    if (RefineableMeshBase* mmesh_pt =
          dynamic_cast<RefineableMeshBase*>(mesh_pt(i_mesh)))
    {
      success_flag += mmesh_pt->unrefine_uniformly();
    }
    else
    {
      oomph_info << "Info/Warning: Mesh cannot be unrefined uniformly "
                 << std::endl;
    }

    // Rebuild the global mesh
    rebuild_global_mesh();

    // Any actions after the adaptation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;

    // Judge success
    if (success_flag > 0)
    {
      return 1;
    }
    else
    {
      return 0;
    }
  }


  //========================================================================
  /// p-unrefine (all) p-refineable (sub)mesh(es) uniformly and rebuild problem;
  /// doc refinement process.
  //========================================================================
  void Problem::p_unrefine_uniformly(DocInfo& doc_info)
  {
    actions_before_adapt();

    // Number of submeshes?
    unsigned n_mesh = nsub_mesh();

    // Single mesh:
    if (n_mesh == 0)
    {
      // Unrefine single mesh uniformly if possible
      if (RefineableMeshBase* mmesh_pt =
            dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
      {
        mmesh_pt->p_unrefine_uniformly(doc_info);
      }
      else
      {
        oomph_info << "Info/Warning: Mesh cannot be p-unrefined uniformly "
                   << std::endl;
      }
    }
    // Multiple submeshes
    else
    {
      // Not tested:
      throw OomphLibError("This functionality has not yet been tested.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
      // Loop over submeshes
      for (unsigned imesh = 0; imesh < n_mesh; imesh++)
      {
        // Unrefine i-th submesh uniformly if possible
        if (RefineableMeshBase* mmesh_pt =
              dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
        {
          mmesh_pt->p_unrefine_uniformly(doc_info);
        }
        else
        {
          oomph_info << "Info/Warning: Cannot p-unrefine mesh " << imesh
                     << std::endl;
        }
      }
      // Rebuild the global mesh
      rebuild_global_mesh();
    }

    // Any actions after the adaptation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }

  //========================================================================
  /// p-unrefine submesh i_mesh uniformly and rebuild problem;
  /// doc refinement process.
  //========================================================================
  void Problem::p_unrefine_uniformly(const unsigned& i_mesh, DocInfo& doc_info)
  {
    actions_before_adapt();

#ifdef PARANOID
    // Number of submeshes?
    if (i_mesh >= nsub_mesh())
    {
      std::ostringstream error_message;
      error_message << "imesh " << i_mesh
                    << " is greater than the number of sub meshes "
                    << nsub_mesh() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Refine single mesh uniformly if possible
    if (RefineableMeshBase* mmesh_pt =
          dynamic_cast<RefineableMeshBase*>(mesh_pt(i_mesh)))
    {
      mmesh_pt->p_unrefine_uniformly(doc_info);
    }
    else
    {
      oomph_info << "Info/Warning: Mesh cannot be p-unrefined uniformly "
                 << std::endl;
    }

    // Rebuild the global mesh
    rebuild_global_mesh();

    // Any actions after the adaptation phase
    actions_after_adapt();

    // Do equation numbering
    oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
  }


  //========================================================================
  /// Do one timestep, dt, forward  using Newton's method with specified
  /// tolerance and linear solver specified via member data.
  /// Keep adapting on all meshes to criteria specified in
  /// these meshes (up to max_adapt adaptations are performed).
  /// If first_timestep==true, re-set initial conditions after mesh adaptation.
  /// Shifting of time can be suppressed by overwriting the
  /// default value of shift (true). [Shifting must be done
  /// if first_timestep==true because we're constantly re-assigning
  /// the initial conditions; if first_timestep==true and shift==false
  /// shifting is performed anyway and a warning is issued.
  //========================================================================
  void Problem::unsteady_newton_solve(const double& dt,
                                      const unsigned& max_adapt,
                                      const bool& first_timestep,
                                      const bool& shift)
  {
    // Do shifting or not?
    bool shift_it = shift;

    // Warning:
    if (first_timestep && (!shift) && (!Default_set_initial_condition_called))
    {
      shift_it = true;
      oomph_info
        << "\n\n===========================================================\n";
      oomph_info << "                  ********  WARNING *********** \n";
      oomph_info
        << "===========================================================\n";
      oomph_info << "Problem::unsteady_newton_solve() called with "
                 << std::endl;
      oomph_info << "first_timestep: " << first_timestep << std::endl;
      oomph_info << "shift: " << shift << std::endl;
      oomph_info << "This doesn't make sense (shifting does have to be done"
                 << std::endl;
      oomph_info << "since we're constantly re-assigning the initial conditions"
                 << std::endl;
      oomph_info
        << "\n===========================================================\n\n";
    }


    // Find the initial time
    double initial_time = time_pt()->time();

    // Max number of solves
    unsigned max_solve = max_adapt + 1;

    // Adaptation loop
    //----------------
    for (unsigned isolve = 0; isolve < max_solve; isolve++)
    {
      // Only adapt after the first solve has been done!
      if (isolve > 0)
      {
        unsigned n_refined;
        unsigned n_unrefined;

        // Adapt problem
        adapt(n_refined, n_unrefined);

#ifdef OOMPH_HAS_MPI
        // Adaptation only converges if ALL the processes have no
        // refinement or unrefinement to perform
        unsigned total_refined = 0;
        unsigned total_unrefined = 0;
        if (Problem_has_been_distributed)
        {
          MPI_Allreduce(&n_refined,
                        &total_refined,
                        1,
                        MPI_UNSIGNED,
                        MPI_SUM,
                        this->communicator_pt()->mpi_comm());
          n_refined = total_refined;
          MPI_Allreduce(&n_unrefined,
                        &total_unrefined,
                        1,
                        MPI_UNSIGNED,
                        MPI_SUM,
                        this->communicator_pt()->mpi_comm());
          n_unrefined = total_unrefined;
        }
#endif

        oomph_info << "---> " << n_refined << " elements were refined, and "
                   << n_unrefined << " were unrefined, in total." << std::endl;

        // Check convergence of adaptation cycle
        if ((n_refined == 0) && (n_unrefined == 0))
        {
          oomph_info << "\n \n Solution is fully converged in "
                     << "Problem::unsteady_newton_solver() \n \n ";
          break;
        }

        // Reset the time
        time_pt()->time() = initial_time;

        // Reset the inital condition on refined meshes. Note that because we
        // have reset the global time to the initial time, the initial
        // conditions are reset at time t=0 rather than at time t=dt
        if (first_timestep)
        {
          // Reset default set_initial_condition has been called flag to false
          Default_set_initial_condition_called = false;

          oomph_info << "Re-setting initial condition " << std::endl;
          set_initial_condition();

          // If the default set_initial_condition function has been called,
          // we must not shift the timevalues on the first timestep, as we
          // will NOT be constantly re-assigning the initial condition
          if (Default_set_initial_condition_called)
          {
            shift_it = false;
          }
        }
      }

      // Now do the actual unsteady timestep
      // If it's the first time around the loop, or the first timestep
      // shift the timevalues, otherwise don't
      // Note: we need to shift if it's the first timestep because
      // we're constantly re-assigning the initial condition above!
      // The only exception to this is if the default set_initial_condition
      // function has been called, in which case we must NOT shift!
      if ((isolve == 0) || (first_timestep))
      {
        Problem::unsteady_newton_solve(dt, shift_it);
      }
      // Subsequent solve: Have shifted already -- don't do it again.
      else
      {
        shift_it = false;
        Problem::unsteady_newton_solve(dt, shift_it);
      }

      if (isolve == max_solve - 1)
      {
        oomph_info
          << std::endl
          << "----------------------------------------------------------"
          << std::endl
          << "Reached max. number of adaptations in \n"
          << "Problem::unsteady_newton_solver().\n"
          << "----------------------------------------------------------"
          << std::endl
          << std::endl;
      }

    } // End of adaptation loop
  }


  //========================================================================
  /// \short Adaptive Newton solver.
  /// The linear solver takes a pointer to the problem (which defines
  /// the Jacobian \b J and the residual Vector \b r) and returns
  /// the solution \b x of the system
  /// \f[ {\bf J} {\bf x} = - \bf{r} \f].
  /// Performs at most max_adapt adaptations on all meshes.
  //========================================================================
  void Problem::newton_solve(const unsigned& max_adapt)
  {
    // Max number of solves
    unsigned max_solve = max_adapt + 1;

    // Adaptation loop
    //----------------
    for (unsigned isolve = 0; isolve < max_solve; isolve++)
    {
      // Only adapt after the first solve has been done!
      if (isolve > 0)
      {
        unsigned n_refined;
        unsigned n_unrefined;

        // Adapt problem
        adapt(n_refined, n_unrefined);

#ifdef OOMPH_HAS_MPI
        // Adaptation only converges if ALL the processes have no
        // refinement or unrefinement to perform
        unsigned total_refined = 0;
        unsigned total_unrefined = 0;
        if (Problem_has_been_distributed)
        {
          MPI_Allreduce(&n_refined,
                        &total_refined,
                        1,
                        MPI_UNSIGNED,
                        MPI_SUM,
                        this->communicator_pt()->mpi_comm());
          n_refined = total_refined;
          MPI_Allreduce(&n_unrefined,
                        &total_unrefined,
                        1,
                        MPI_UNSIGNED,
                        MPI_SUM,
                        this->communicator_pt()->mpi_comm());
          n_unrefined = total_unrefined;
        }
#endif

        oomph_info << "---> " << n_refined << " elements were refined, and "
                   << n_unrefined << " were unrefined"
#ifdef OOMPH_HAS_MPI
                   << ", in total (over all processors).\n";
#else
                   << ".\n";
#endif


        // Check convergence of adaptation cycle
        if ((n_refined == 0) && (n_unrefined == 0))
        {
          oomph_info << "\n \n Solution is fully converged in "
                     << "Problem::newton_solver(). \n \n ";
          break;
        }
      }


      // Do actual solve
      //----------------
      {
        // Now update anything that needs updating
        // NOT NEEDED -- IS CALLED IN newton_solve BELOW! #
        // actions_before_newton_solve();

        try
        {
          // Solve the non-linear problem for this timestep with Newton's method
          newton_solve();
        }
        // Catch any exceptions thrown in the Newton solver
        catch (NewtonSolverError& error)
        {
          oomph_info << std::endl
                     << "USER-DEFINED ERROR IN NEWTON SOLVER " << std::endl;
          // Check to see whether we have reached Max_iterations
          if (error.iterations == Max_newton_iterations)
          {
            oomph_info << "MAXIMUM NUMBER OF ITERATIONS (" << error.iterations
                       << ") REACHED WITHOUT CONVERGENCE " << std::endl;
          }
          // If not, it must be that we have exceeded the maximum residuals
          else
          {
            oomph_info << "MAXIMUM RESIDUALS: " << error.maxres
                       << "EXCEEDS PREDEFINED MAXIMUM " << Max_residuals
                       << std::endl;
          }

          // Die horribly!!
          std::ostringstream error_stream;
          error_stream << "Error occured in adaptive Newton solver. "
                       << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        // Now update anything that needs updating
        // NOT NEEDED -- WAS CALLED IN newton_solve ABOVE
        // !actions_after_newton_solve();

      } // End of solve block


      if (isolve == max_solve - 1)
      {
        oomph_info
          << std::endl
          << "----------------------------------------------------------"
          << std::endl
          << "Reached max. number of adaptations in \n"
          << "Problem::newton_solver().\n"
          << "----------------------------------------------------------"
          << std::endl
          << std::endl;
      }

    } // End of adaptation loop
  }

  //========================================================================
  /// Delete any external storage for any submeshes
  /// NB this would ordinarily take place within the adaptation procedure
  /// for each submesh (See RefineableMesh::adapt_mesh(...)), but there
  /// are instances where the actions_before/after_adapt routines are used
  /// and no adaptive routines are called in between (e.g. when doc-ing
  /// errors at the end of an adaptive newton solver)
  //========================================================================
  void Problem::delete_all_external_storage()
  {
    // Number of submeshes
    unsigned n_mesh = nsub_mesh();

    // External storage will only exist if there is more than one (sub)mesh
    if (n_mesh > 1)
    {
      for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
      {
        mesh_pt(i_mesh)->delete_all_external_storage();
      }
    }
  }


#ifdef OOMPH_HAS_MPI

  //====================================================================
  /// Get all the halo data stored on this processor and store pointers
  /// to the data in a map, indexed by the gobal eqn number
  //====================================================================
  void Problem::get_all_halo_data(std::map<unsigned, double*>& map_of_halo_data)
  {
    // Halo data is stored in the meshes, so kick the problem down to that
    // level

    // Find the number of meshes
    unsigned n_mesh = this->nsub_mesh();
    // If there are no submeshes it's only the main mesh
    if (n_mesh == 0)
    {
      mesh_pt()->get_all_halo_data(map_of_halo_data);
    }
    // Otherwise loop over all the submeshes
    else
    {
      for (unsigned imesh = 0; imesh < n_mesh; ++imesh)
      {
        mesh_pt(imesh)->get_all_halo_data(map_of_halo_data);
      }
    }
  }


  //========================================================================
  /// Check the halo/haloed/shared node/element schemes.
  //========================================================================
  void Problem::check_halo_schemes(DocInfo& doc_info)
  {
    // The bulk of the stuff that was in this routine is mesh-based, and
    // should therefore drop into the Mesh base class.  All that needs to remain
    // here is a "wrapper" which calls the function dependent upon the number
    // of (sub)meshes that may have been distributed.

    unsigned n_mesh = nsub_mesh();

    if (n_mesh == 0)
    {
      oomph_info << "Checking halo schemes on single mesh" << std::endl;
      doc_info.label() = "_one_and_only_mesh_";
      mesh_pt()->check_halo_schemes(doc_info,
                                    Max_permitted_error_for_halo_check);
    }
    else // there are submeshes
    {
      for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
      {
        oomph_info << "Checking halo schemes on submesh " << i_mesh
                   << std::endl;
        std::stringstream tmp;
        tmp << "_mesh" << i_mesh << "_";
        doc_info.label() = tmp.str();
        mesh_pt(i_mesh)->check_halo_schemes(doc_info,
                                            Max_permitted_error_for_halo_check);
      }
    }
  }


  //========================================================================
  /// Synchronise all dofs by calling the appropriate synchronisation
  /// routines for all meshes and the assembly handler
  //========================================================================
  void Problem::synchronise_all_dofs()
  {
    // Synchronise dofs themselves
    bool do_halos = true;
    bool do_external_halos = false;
    this->synchronise_dofs(do_halos, do_external_halos);


    do_halos = false;
    do_external_halos = true;
    this->synchronise_dofs(do_halos, do_external_halos);

    // Now perform any synchronisation required by the assembly handler
    this->assembly_handler_pt()->synchronise();
  }


  //========================================================================
  /// Synchronise the degrees of freedom by overwriting
  /// the haloed values with their non-halo counterparts held
  /// on other processors. Bools control if we deal with data associated with
  /// external halo/ed elements/nodes or the "normal" halo/ed ones.
  //========================================================================
  void Problem::synchronise_dofs(const bool& do_halos,
                                 const bool& do_external_halos)
  {
    // Do we have submeshes?
    unsigned n_mesh_loop = 1;
    unsigned nmesh = nsub_mesh();
    if (nmesh > 0)
    {
      n_mesh_loop = nmesh;
    }

    // Local storage for number of processors and current processor
    const int n_proc = this->communicator_pt()->nproc();

    // If only one processor then return
    if (n_proc == 1)
    {
      return;
    }

    const int my_rank = this->communicator_pt()->my_rank();

    // Storage for number of data to be sent to each processor
    Vector<int> send_n(n_proc, 0);

    // Storage for all values to be sent to all processors
    Vector<double> send_data;

    // Start location within send_data for data to be sent to each processor
    Vector<int> send_displacement(n_proc, 0);

    // Loop over all processors
    for (int rank = 0; rank < n_proc; rank++)
    {
      // Set the offset for the current processor
      send_displacement[rank] = send_data.size();

      // Don't bother to do anything if the processor in the loop is the
      // current processor
      if (rank != my_rank)
      {
        // Deal with sub-meshes one-by-one if required
        Mesh* my_mesh_pt = 0;

        // Loop over submeshes
        for (unsigned imesh = 0; imesh < n_mesh_loop; imesh++)
        {
          if (nmesh == 0)
          {
            my_mesh_pt = mesh_pt();
          }
          else
          {
            my_mesh_pt = mesh_pt(imesh);
          }

          if (do_halos)
          {
            // How many of my nodes are haloed by the processor whose values
            // are updated?
            unsigned n_nod = my_mesh_pt->nhaloed_node(rank);
            for (unsigned n = 0; n < n_nod; n++)
            {
              // Add the data for each haloed node to the vector
              my_mesh_pt->haloed_node_pt(rank, n)->add_values_to_vector(
                send_data);
            }

            // Now loop over haloed elements and prepare to add their
            // internal data to the big vector to be sent
            Vector<GeneralisedElement*> haloed_elem_pt =
              my_mesh_pt->haloed_element_pt(rank);
            unsigned nelem_haloed = haloed_elem_pt.size();
            for (unsigned e = 0; e < nelem_haloed; e++)
            {
              haloed_elem_pt[e]->add_internal_data_values_to_vector(send_data);
            }
          }

          if (do_external_halos)
          {
            // How many of my nodes are externally haloed by the processor whose
            // values are updated?  NB these nodes are on the external mesh.
            unsigned n_ext_nod = my_mesh_pt->nexternal_haloed_node(rank);
            for (unsigned n = 0; n < n_ext_nod; n++)
            {
              // Add data from each external haloed node to the vector
              my_mesh_pt->external_haloed_node_pt(rank, n)
                ->add_values_to_vector(send_data);
            }

            // Now loop over haloed elements and prepare to send internal data
            unsigned next_elem_haloed =
              my_mesh_pt->nexternal_haloed_element(rank);
            for (unsigned e = 0; e < next_elem_haloed; e++)
            {
              my_mesh_pt->external_haloed_element_pt(rank, e)
                ->add_internal_data_values_to_vector(send_data);
            }
          }
        } // end of loop over meshes
      }

      // Find the number of data added to the vector
      send_n[rank] = send_data.size() - send_displacement[rank];
    }


    // Storage for the number of data to be received from each processor
    Vector<int> receive_n(n_proc, 0);

    // Now send numbers of data to be sent between all processors
    MPI_Alltoall(&send_n[0],
                 1,
                 MPI_INT,
                 &receive_n[0],
                 1,
                 MPI_INT,
                 this->communicator_pt()->mpi_comm());

    // We now prepare the data to be received
    // by working out the displacements from the received data
    Vector<int> receive_displacement(n_proc, 0);
    int receive_data_count = 0;
    for (int rank = 0; rank < n_proc; ++rank)
    {
      // Displacement is number of data received so far
      receive_displacement[rank] = receive_data_count;
      receive_data_count += receive_n[rank];
    }

    // Now resize the receive buffer for all data from all processors
    // Make sure that it has a size of at least one
    if (receive_data_count == 0)
    {
      ++receive_data_count;
    }
    Vector<double> receive_data(receive_data_count);

    // Make sure that the send buffer has size at least one
    // so that we don't get a segmentation fault
    if (send_data.size() == 0)
    {
      send_data.resize(1);
    }

    // Now send the data between all the processors
    MPI_Alltoallv(&send_data[0],
                  &send_n[0],
                  &send_displacement[0],
                  MPI_DOUBLE,
                  &receive_data[0],
                  &receive_n[0],
                  &receive_displacement[0],
                  MPI_DOUBLE,
                  this->communicator_pt()->mpi_comm());

    // Now use the received data to update the halo nodes
    for (int send_rank = 0; send_rank < n_proc; send_rank++)
    {
      // Don't bother to do anything for the processor corresponding to the
      // current processor or if no data were received from this processor
      if ((send_rank != my_rank) && (receive_n[send_rank] != 0))
      {
        // Counter for the data within the large array
        unsigned count = receive_displacement[send_rank];

        // Deal with sub-meshes one-by-one if required
        Mesh* my_mesh_pt = 0;

        // Loop over submeshes
        for (unsigned imesh = 0; imesh < n_mesh_loop; imesh++)
        {
          if (nmesh == 0)
          {
            my_mesh_pt = mesh_pt();
          }
          else
          {
            my_mesh_pt = mesh_pt(imesh);
          }

          if (do_halos)
          {
            // How many of my nodes are halos whose non-halo counter
            // parts live on processor send_rank?
            unsigned n_nod = my_mesh_pt->nhalo_node(send_rank);
            for (unsigned n = 0; n < n_nod; n++)
            {
              // Read in values for each halo node
              my_mesh_pt->halo_node_pt(send_rank, n)
                ->read_values_from_vector(receive_data, count);
            }

            // Get number of halo elements whose non-halo is
            // on process send_rank
            Vector<GeneralisedElement*> halo_elem_pt =
              my_mesh_pt->halo_element_pt(send_rank);

            unsigned nelem_halo = halo_elem_pt.size();
            for (unsigned e = 0; e < nelem_halo; e++)
            {
              halo_elem_pt[e]->read_internal_data_values_from_vector(
                receive_data, count);
            }
          }

          if (do_external_halos)
          {
            // How many of my nodes are external halos whose external non-halo
            // counterparts live on processor send_rank?
            unsigned n_ext_nod = my_mesh_pt->nexternal_halo_node(send_rank);

            // Copy into the values of the external halo nodes
            // on the present processors
            for (unsigned n = 0; n < n_ext_nod; n++)
            {
              // Read the data from the array into each halo node
              my_mesh_pt->external_halo_node_pt(send_rank, n)
                ->read_values_from_vector(receive_data, count);
            }

            // Get number of halo elements whose non-halo is
            // on process send_rank
            unsigned next_elem_halo =
              my_mesh_pt->nexternal_halo_element(send_rank);
            for (unsigned e = 0; e < next_elem_halo; e++)
            {
              my_mesh_pt->external_halo_element_pt(send_rank, e)
                ->read_internal_data_values_from_vector(receive_data, count);
            }
          }

        } // end of loop over meshes
      }
    } // End of data is received
  } // End of synchronise


  //========================================================================
  ///  Synchronise equation numbers and return the total
  /// number of degrees of freedom in the overall problem
  //========================================================================
  long Problem::synchronise_eqn_numbers(const bool& assign_local_eqn_numbers)
  {
    // number of equations on this processor, which at this stage is only known
    // by counting the number of dofs that have been added to the problem
    unsigned my_n_eqn = Dof_pt.size();

    // my rank
    unsigned my_rank = Communicator_pt->my_rank();

    // number of processors
    unsigned nproc = Communicator_pt->nproc();

    //  // Time alternative communication
    //  Vector<unsigned> n_eqn(nproc);
    //  {
    //   double t_start = TimingHelpers::timer();

    //   // Gather numbers of equations (enumerated independently on all procs)
    //   MPI_Allgather(&my_n_eqn,1,MPI_UNSIGNED,&n_eqn[0],
    //                 1,MPI_INT,Communicator_pt->mpi_comm());

    //   double t_end = TimingHelpers::timer();
    //   oomph_info << "Time for allgather-based exchange of eqn numbers: "
    //              << t_end-t_start << std::endl;
    //  }

    double t_start = TimingHelpers::timer();

    // send my_n_eqn to with rank greater than my_rank
    unsigned n_send = nproc - my_rank - 1;
    Vector<MPI_Request> send_req(n_send);
    for (unsigned p = my_rank + 1; p < nproc; p++)
    {
      MPI_Isend(&my_n_eqn,
                1,
                MPI_UNSIGNED,
                p,
                0,
                Communicator_pt->mpi_comm(),
                &send_req[p - my_rank - 1]);
    }

    // recv n_eqn from processors with rank less than my_rank
    Vector<unsigned> n_eqn_on_proc(my_rank);
    for (unsigned p = 0; p < my_rank; p++)
    {
      MPI_Recv(&n_eqn_on_proc[p],
               1,
               MPI_UNSIGNED,
               p,
               0,
               Communicator_pt->mpi_comm(),
               MPI_STATUS_IGNORE);
    }

    double t_end = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for send and receive stuff: " << t_end - t_start
                 << std::endl;
      t_start = TimingHelpers::timer();
    }

    // determine the number of equation on processors with rank
    // less than my_rank
    unsigned my_eqn_num_base = 0;
    for (unsigned p = 0; p < my_rank; p++)
    {
      my_eqn_num_base += n_eqn_on_proc[p];
      //   if (n_eqn_on_proc[p]!=n_eqn[p])
      //     {
      //      std::cout << "proc " << my_rank << "clash in eqn numbers: "
      //                << p << " " << n_eqn_on_proc[p] << " " << n_eqn[p]
      //                << std::endl;
      //     }
    }

    // Loop over all internal data (on elements) and bump up their
    // equation numbers if they exist
    unsigned nelem = mesh_pt()->nelement();
    for (unsigned e = 0; e < nelem; e++)
    {
      GeneralisedElement* el_pt = mesh_pt()->element_pt(e);

      unsigned nintern_data = el_pt->ninternal_data();
      for (unsigned iintern = 0; iintern < nintern_data; iintern++)
      {
        Data* int_data_pt = el_pt->internal_data_pt(iintern);
        unsigned nval = int_data_pt->nvalue();
        for (unsigned ival = 0; ival < nval; ival++)
        {
          int old_eqn_number = int_data_pt->eqn_number(ival);
          if (old_eqn_number >= 0) // i.e. it's being used
          {
            // Bump up eqn number
            int new_eqn_number = old_eqn_number + my_eqn_num_base;
            int_data_pt->eqn_number(ival) = new_eqn_number;
          }
        }
      }
    }

    // Loop over all nodes on current processor and bump up their
    // equation numbers if they're not pinned!
    unsigned nnod = mesh_pt()->nnode();
    for (unsigned j = 0; j < nnod; j++)
    {
      Node* nod_pt = mesh_pt()->node_pt(j);

      // loop over ALL eqn numbers - variable number of values
      unsigned nval = nod_pt->nvalue();

      for (unsigned ival = 0; ival < nval; ival++)
      {
        int old_eqn_number = nod_pt->eqn_number(ival);
        // Include all eqn numbers
        if (old_eqn_number >= 0)
        {
          // Bump up eqn number
          int new_eqn_number = old_eqn_number + my_eqn_num_base;
          nod_pt->eqn_number(ival) = new_eqn_number;
        }
      }

      // Is this a solid node? If so, need to bump up its equation number(s)
      SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);

      if (solid_nod_pt != 0)
      {
        // Find equation numbers
        unsigned nval = solid_nod_pt->variable_position_pt()->nvalue();
        for (unsigned ival = 0; ival < nval; ival++)
        {
          int old_eqn_number =
            solid_nod_pt->variable_position_pt()->eqn_number(ival);
          // include all eqn numbers

          if (old_eqn_number >= 0)
          {
            // Bump up eqn number
            int new_eqn_number = old_eqn_number + my_eqn_num_base;
            solid_nod_pt->variable_position_pt()->eqn_number(ival) =
              new_eqn_number;
          }
        }
      }
    }

    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for bumping: " << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }


    // Now copy the haloed eqn numbers across
    // This has to include the internal data equation numbers as well
    // as the solid node equation numbers
    bool do_halos = true;
    bool do_external_halos = false;
    copy_haloed_eqn_numbers_helper(do_halos, do_external_halos);

    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for copy_haloed_eqn_numbers_helper for halos: "
                 << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }

    // Now do external halo stuff
    do_halos = false;
    do_external_halos = true;
    copy_haloed_eqn_numbers_helper(do_halos, do_external_halos);

    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info
        << "Time for copy_haloed_eqn_numbers_helper for external halos: "
        << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }

    // Now the global equation numbers have been updated.
    //---------------------------------------------------
    // Setup the local equation numbers again.
    //----------------------------------------
    if (assign_local_eqn_numbers)
    {
      // Loop over the submeshes: Note we need to call the submeshes' own
      // assign_*_eqn_number() otherwise we miss additional functionality
      // that is implemented (e.g.) in SolidMeshes!
      unsigned n_sub_mesh = nsub_mesh();
      if (n_sub_mesh == 0)
      {
        mesh_pt()->assign_local_eqn_numbers(Store_local_dof_pt_in_elements);
      }
      else
      {
        for (unsigned i = 0; i < n_sub_mesh; i++)
        {
          mesh_pt(i)->assign_local_eqn_numbers(Store_local_dof_pt_in_elements);
        }
      }
    }

    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for assign_local_eqn_numbers in sync: "
                 << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }

    // wait for the sends to complete
    if (n_send > 0)
    {
      Vector<MPI_Status> send_status(n_send);
      MPI_Waitall(n_send, &send_req[0], &send_status[0]);
    }

    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for waitall: " << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }

    // build the Dof distribution pt
    Dof_distribution_pt->build(Communicator_pt, my_eqn_num_base, my_n_eqn);

    // and return the total number of equations in the problem
    return (long)Dof_distribution_pt->nrow();
  }


  //=======================================================================
  /// A private helper function to
  /// copy the haloed equation numbers into the halo equation numbers,
  /// either for the problem's one and only mesh or for all of its
  /// submeshes. Bools control if we deal with data associated with
  /// external halo/ed elements/nodes or the "normal" halo/ed ones.
  //===================================================================
  void Problem::copy_haloed_eqn_numbers_helper(const bool& do_halos,
                                               const bool& do_external_halos)
  {
    // Do we have submeshes?
    unsigned n_mesh_loop = 1;
    unsigned nmesh = nsub_mesh();
    if (nmesh > 0)
    {
      n_mesh_loop = nmesh;
    }

    // Storage for number of processors and current processor
    int n_proc = this->communicator_pt()->nproc();

    // If only one processor then return
    if (n_proc == 1)
    {
      return;
    }
    int my_rank = this->communicator_pt()->my_rank();

    // Storage for number of data to be sent to each processor
    Vector<int> send_n(n_proc, 0);
    // Storage for all equation numbers to be sent to all processors
    Vector<long> send_data;
    // Start location within send_data for data to be sent to each processor
    Vector<int> send_displacement(n_proc, 0);


    // Loop over all processors whose eqn numbers are to be updated
    for (int rank = 0; rank < n_proc; rank++)
    {
      // Set the displacement of the current processor in the loop
      send_displacement[rank] = send_data.size();

      // If I'm not the processor whose halo eqn numbers are updated,
      // some of my nodes may be haloed: Stick their
      // eqn numbers into the vector
      if (rank != my_rank)
      {
        // Deal with sub-meshes one-by-one if required
        Mesh* my_mesh_pt = 0;

        // Loop over submeshes
        for (unsigned imesh = 0; imesh < n_mesh_loop; imesh++)
        {
          if (nmesh == 0)
          {
            my_mesh_pt = mesh_pt();
          }
          else
          {
            my_mesh_pt = mesh_pt(imesh);
          }

          if (do_halos)
          {
            // Add equation numbers for each haloed node
            unsigned n_nod = my_mesh_pt->nhaloed_node(rank);
            for (unsigned n = 0; n < n_nod; n++)
            {
              my_mesh_pt->haloed_node_pt(rank, n)->add_eqn_numbers_to_vector(
                send_data);
            }

            // Add the equation numbers associated with internal data
            // in the haloed elements
            Vector<GeneralisedElement*> haloed_elem_pt =
              my_mesh_pt->haloed_element_pt(rank);
            unsigned nelem_haloed = haloed_elem_pt.size();
            for (unsigned e = 0; e < nelem_haloed; e++)
            {
              haloed_elem_pt[e]->add_internal_eqn_numbers_to_vector(send_data);
            }
          }

          if (do_external_halos)
          {
            // Add equation numbers associated with external haloed nodes
            unsigned n_ext_nod = my_mesh_pt->nexternal_haloed_node(rank);
            for (unsigned n = 0; n < n_ext_nod; n++)
            {
              my_mesh_pt->external_haloed_node_pt(rank, n)
                ->add_eqn_numbers_to_vector(send_data);
            }

            // Add the equation numbers associated with internal data in
            // each external haloed element
            unsigned next_elem_haloed =
              my_mesh_pt->nexternal_haloed_element(rank);
            for (unsigned e = 0; e < next_elem_haloed; e++)
            {
              // how many internal data values for this element?
              my_mesh_pt->external_haloed_element_pt(rank, e)
                ->add_internal_eqn_numbers_to_vector(send_data);
            }
          }

        } // end of loop over meshes
      }

      // Find the number of data added to the vector by this processor
      send_n[rank] = send_data.size() - send_displacement[rank];
    }

    // Storage for the number of data to be received from each processor
    Vector<int> receive_n(n_proc, 0);

    // Communicate all numbers of data to be sent between all processors
    MPI_Alltoall(&send_n[0],
                 1,
                 MPI_INT,
                 &receive_n[0],
                 1,
                 MPI_INT,
                 this->communicator_pt()->mpi_comm());

    // We now prepare the data to be received
    // by working out the displacements from the received data
    Vector<int> receive_displacement(n_proc, 0);
    int receive_data_count = 0;
    for (int rank = 0; rank < n_proc; ++rank)
    {
      // Displacement is number of data received so far
      receive_displacement[rank] = receive_data_count;
      receive_data_count += receive_n[rank];
    }

    // Now resize the receive buffer
    // Make sure that it has a size of at least one
    if (receive_data_count == 0)
    {
      ++receive_data_count;
    }
    Vector<long> receive_data(receive_data_count);

    // Make sure that the send buffer has size at least one
    // so that we don't get a segmentation fault
    if (send_data.size() == 0)
    {
      send_data.resize(1);
    }

    // Now send the data between all the processors
    MPI_Alltoallv(&send_data[0],
                  &send_n[0],
                  &send_displacement[0],
                  MPI_LONG,
                  &receive_data[0],
                  &receive_n[0],
                  &receive_displacement[0],
                  MPI_LONG,
                  this->communicator_pt()->mpi_comm());


    // Loop over all other processors to receive their
    // eqn numbers
    for (int send_rank = 0; send_rank < n_proc; send_rank++)
    {
      // Don't do anything for the processor corresponding to the
      // current processor or if no data were received from this processor
      if ((send_rank != my_rank) && (receive_n[send_rank] != 0))
      {
        // Counter for the data within the large array
        unsigned count = receive_displacement[send_rank];

        // Deal with sub-meshes one-by-one if required
        Mesh* my_mesh_pt = 0;

        // Loop over submeshes
        for (unsigned imesh = 0; imesh < n_mesh_loop; imesh++)
        {
          if (nmesh == 0)
          {
            my_mesh_pt = mesh_pt();
          }
          else
          {
            my_mesh_pt = mesh_pt(imesh);
          }

          if (do_halos)
          {
            // How many of my nodes are halos whose non-halo counter
            // parts live on processor send_rank?
            unsigned n_nod = my_mesh_pt->nhalo_node(send_rank);
            for (unsigned n = 0; n < n_nod; n++)
            {
              // Generalise to variable number of values per node
              my_mesh_pt->halo_node_pt(send_rank, n)
                ->read_eqn_numbers_from_vector(receive_data, count);
            }

            // Get number of halo elements whose non-halo is on
            // process send_rank
            Vector<GeneralisedElement*> halo_elem_pt =
              my_mesh_pt->halo_element_pt(send_rank);
            unsigned nelem_halo = halo_elem_pt.size();
            for (unsigned e = 0; e < nelem_halo; e++)
            {
              halo_elem_pt[e]->read_internal_eqn_numbers_from_vector(
                receive_data, count);
            }
          }

          if (do_external_halos)
          {
            // How many of my nodes are external halos whose external non-halo
            // counterparts live on processor send_rank?
            unsigned n_ext_nod = my_mesh_pt->nexternal_halo_node(send_rank);
            for (unsigned n = 0; n < n_ext_nod; n++)
            {
              my_mesh_pt->external_halo_node_pt(send_rank, n)
                ->read_eqn_numbers_from_vector(receive_data, count);
            }

            // Get number of external halo elements whose external haloed
            // counterpart is on process send_rank
            unsigned next_elem_halo =
              my_mesh_pt->nexternal_halo_element(send_rank);
            for (unsigned e = 0; e < next_elem_halo; e++)
            {
              my_mesh_pt->external_halo_element_pt(send_rank, e)
                ->read_internal_eqn_numbers_from_vector(receive_data, count);
            }
          }

        } // end of loop over meshes
      }
    } // End of loop over processors
  }

  //==========================================================================
  /// Balance the load of a (possibly non-uniformly refined) problem that has
  /// already been distributed, by re-distributing elements over the processors.
  /// Produce explicit stats of load balancing process if boolean, report_stats,
  /// is set to true and doc various bits of data (mainly for debugging)
  /// in directory specified by DocInfo object.
  //==========================================================================
  void Problem::load_balance(
    DocInfo& doc_info,
    const bool& report_stats,
    const Vector<unsigned>& input_target_domain_for_local_non_halo_element)
  {
    double start_t = TimingHelpers::timer();

    // Number of processes
    const unsigned n_proc = this->communicator_pt()->nproc();

    // Don't do anything if this is a single-process job
    if (n_proc == 1)
    {
      if (report_stats)
      {
        std::ostringstream warn_message;
        warn_message << "WARNING: You've tried to load balance a problem over\n"
                     << "only one processor: ignoring your request.\n";
        OomphLibWarning(warn_message.str(),
                        "Problem::load_balance()",
                        OOMPH_EXCEPTION_LOCATION);
      }
    }
    // Multiple processors
    else
    {
      // This will only work if the problem has already been distributed
      if (!Problem_has_been_distributed)
      {
        // Throw an error
        std::ostringstream error_stream;
        error_stream << "You have called Problem::load_balance()\n"
                     << "on a non-distributed problem. This doesn't\n"
                     << "make sense -- go distribute your problem first."
                     << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Timings
      double t_start = 0.0;
      double t_metis = 0.0;
      double t_partition = 0.0;
      double t_distribute = 0.0;
      double t_refine = 0.0;
      double t_copy_solution = 0.0;

      if (report_stats)
      {
        t_start = TimingHelpers::timer();
      }


#ifdef PARANOID
      unsigned old_ndof = ndof();
#endif

      // Store pointers to the old mesh(es) so we retain a handle
      //---------------------------------------------------------
      // to them for deletion
      //---------------------
      Vector<Mesh*> old_mesh_pt;
      unsigned n_mesh = nsub_mesh();
      if (n_mesh == 0)
      {
        // Resize the container
        old_mesh_pt.resize(1);
        old_mesh_pt[0] = mesh_pt();
      }
      else
      {
        // Resize the container
        old_mesh_pt.resize(n_mesh);
        for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
        {
          old_mesh_pt[i_mesh] = mesh_pt(i_mesh);
        }
      }


      // Partition the global mesh in its current state
      //-----------------------------------------------

      // target_domain_for_local_non_halo_element[e] contains the number
      // of the domain [0,1,...,nproc-1] to which non-halo element e on THE
      // CURRENT PROCESSOR ONLY has been assigned. The order of the non-halo
      // elements is the same as in the Problem's mesh, with the halo
      // elements being skipped.
      Vector<unsigned> target_domain_for_local_non_halo_element;

      // Do any of the processors want to go through externally imposed
      // partitioning? If so, we'd better do it here too (even if the processor
      // is empty, e.g. following a restart on a larger number of procs) or
      // we hang.
      unsigned local_ntarget =
        input_target_domain_for_local_non_halo_element.size();
      unsigned global_ntarget = 0;
      MPI_Allreduce(&local_ntarget,
                    &global_ntarget,
                    1,
                    MPI_UNSIGNED,
                    MPI_MAX,
                    Communicator_pt->mpi_comm());

      // External prescribed partitioning
      if (global_ntarget > 0)
      {
        target_domain_for_local_non_halo_element =
          input_target_domain_for_local_non_halo_element;
      }
      else
      {
        // Metis does not always produce repeatable results which is
        // a disaster for validation runs -- this bypasses metis and
        // comes up with a stupid but repeatable partioning.
        if (Use_default_partition_in_load_balance)
        {
          // Bypass METIS to perform the partitioning
          unsigned objective = 0;
          bool bypass_metis = true;
          METIS::partition_distributed_mesh(
            this,
            objective,
            target_domain_for_local_non_halo_element,
            bypass_metis);
        }
        else
        {
          // Use METIS to perform the partitioning
          unsigned objective = 0;
          METIS::partition_distributed_mesh(
            this, objective, target_domain_for_local_non_halo_element);
        }
      }

      if (report_stats)
      {
        t_metis = TimingHelpers::timer();
      }

      // Setup map linking element with target domain
      std::map<GeneralisedElement*, unsigned>
        target_domain_for_local_non_halo_element_map;
      unsigned n_elem = mesh_pt()->nelement();
      unsigned count_non_halo_el = 0;
      for (unsigned e = 0; e < n_elem; e++)
      {
        GeneralisedElement* el_pt = mesh_pt()->element_pt(e);
        if (!el_pt->is_halo())
        {
          target_domain_for_local_non_halo_element_map[el_pt] =
            target_domain_for_local_non_halo_element[count_non_halo_el];
          count_non_halo_el++;
        }
      }

      // Load balancing is equivalent to distribution so call the
      // appropriate "actions before". NOTE: This acts on the
      // current, refined, distributed, etc. problem object
      // before it's being wiped. This step is therefore not
      // a duplicate of the call below, which acts on the
      // new, not-yet refined, distributed etc. problem!
      actions_before_distribute();

      // Re-setup target domains for remaining elements (FaceElements
      // are likely to have been stripped out in actions_before_distribute()
      n_elem = mesh_pt()->nelement();
      target_domain_for_local_non_halo_element.clear();
      target_domain_for_local_non_halo_element.reserve(n_elem);
      count_non_halo_el = 0;
      for (unsigned e = 0; e < n_elem; e++)
      {
        GeneralisedElement* el_pt = mesh_pt()->element_pt(e);
        if (!el_pt->is_halo())
        {
          target_domain_for_local_non_halo_element.push_back(
            target_domain_for_local_non_halo_element_map[el_pt]);
        }
      } // for (e < n_elem)

      // Re-setup the number of sub-meshes since some of them may have
      // been stripped out in actions_before_distribute(), but save the
      // number of old sub-meshes
      const unsigned n_old_sub_meshes = n_mesh;
      n_mesh = nsub_mesh();

      // Now get the target domains for each of the submeshes, we only
      // get the target domains for the nonhalo elements
      Vector<Vector<unsigned>> target_domain_for_local_non_halo_element_submesh(
        n_mesh);
      // If we have no sub-meshes then we do not need to copy the target areas
      // of the submeshes
      if (n_mesh != 0)
      {
        // Counter to copy the target domains from the global vector
        unsigned count_td = 0;
        for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
        {
          // Get the number of elements (considering halo elements)
          const unsigned nsub_ele = mesh_pt(i_mesh)->nelement();
          // Now copy that number of data from the global target domains
          for (unsigned i = 0; i < nsub_ele; i++)
          {
            // Get the element
            GeneralisedElement* ele_pt = mesh_pt(i_mesh)->element_pt(i);
            // ... and check if it is a nonhalo element
            if (!ele_pt->is_halo())
            {
              // Get the target domain for the current element
              const unsigned target_domain =
                target_domain_for_local_non_halo_element[count_td++];
              // Add the target domain for the nonhalo element in the
              // submesh
              target_domain_for_local_non_halo_element_submesh[i_mesh]
                .push_back(target_domain);
            } // if (!ele_pt->is_halo())
          } // for (i < nsub_ele)
        } // for (imesh < n_mesh)

#ifdef PARANOID
        // Check that the total number of copied data be the same as the
        // total number of nonhalo elements in the (sub)-mesh(es)
        const unsigned ntarget_domain =
          target_domain_for_local_non_halo_element.size();
        if (count_td != ntarget_domain)
        {
          std::ostringstream error_stream;
          error_stream
            << "The number of nonhalo elements (" << count_td
            << ") found in (all)\n"
            << "the (sub)-mesh(es) is different from the number of target "
               "domains\n"
            << "(" << ntarget_domain << ") for the nonhalo elements.\n"
            << "Please ensure that you called the rebuild_global_mesh() method "
            << "after the\npossible deletion of FaceElements in "
            << "actions_before_distribute()!!!\n\n";
          throw OomphLibError(error_stream.str(),
                              "Problem::load_balance()",
                              OOMPH_EXCEPTION_LOCATION);
        } // if (count_td != ntarget_domain)
#endif

      } // if (n_mesh != 0)

      // Check if we have different type of submeshes (unstructured
      // and/or structured). Identify to which type each submesh belongs.
      // If we have only one mesh then identify to which type that mesh
      // belongs.

      // The load balancing strategy acts in the structured meshes and
      // then acts in the unstructured meshes

      // Vector to temporaly store pointers to unstructured meshes
      // (TriangleMeshBase)
      Vector<TriangleMeshBase*> unstructured_mesh_pt;
      std::vector<bool> is_unstructured_mesh;

      // Flag to indicate that there are unstructured meshes as part of
      // the problem
      bool are_there_unstructured_meshes = false;

      // We have only one mesh
      if (n_mesh == 0)
      {
        // Check if it is a TriangleMeshBase mesh
        if (TriangleMeshBase* tri_mesh_pt =
              dynamic_cast<TriangleMeshBase*>(old_mesh_pt[0]))
        {
          // Add the pointer to the unstructured meshes container
          unstructured_mesh_pt.push_back(tri_mesh_pt);
          // Indicate that it is an unstructured mesh
          is_unstructured_mesh.push_back(true);
          // Indicate that there are unstructured meshes as part of the
          // problem
          are_there_unstructured_meshes = true;
        }
        else
        {
          // Add the pointer to the unstructured meshes container (null
          // pointer)
          unstructured_mesh_pt.push_back(tri_mesh_pt);
          // Indicate that it is not an unstructured mesh
          is_unstructured_mesh.push_back(false);
        }
      } // if (n_mesh == 0)
      else // We have sub-meshes
      {
        // Check which sub-meshes are unstructured meshes (work with the
        // old sub-meshes number)
        for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
        {
          // Is it a TriangleMeshBase mesh
          if (TriangleMeshBase* tri_mesh_pt =
                dynamic_cast<TriangleMeshBase*>(old_mesh_pt[i_mesh]))
          {
            // Add the pointer to the unstructured meshes container
            unstructured_mesh_pt.push_back(tri_mesh_pt);
            // Indicate that it is an unstructured mesh
            is_unstructured_mesh.push_back(true);
            // Indicate that there are unstructured meshes as part of the
            // problem
            are_there_unstructured_meshes = true;
          }
          else
          {
            // Add the pointer to the unstructured meshes container (null
            // pointer)
            unstructured_mesh_pt.push_back(tri_mesh_pt);
            // Indicate that it is not an unstructured mesh
            is_unstructured_mesh.push_back(false);
          }
        } // for (i_mesh < n_mesh)
      } // else if (n_mesh == 0) // We have sub-meshes

      // Extract data to be sent to various processors after the
      //--------------------------------------------------------
      // problem has been rebuilt/re-distributed
      //----------------------------------------

      // Storage for number of data to be sent to each processor
      Vector<int> send_n(n_proc, 0);

      // Storage for all values to be sent to all processors
      Vector<double> send_data;

      // Start location within send_data for data to be sent to each processor
      Vector<int> send_displacement(n_proc, 0);

      // Old and new domains for each base element (available for all, for
      // convenience)
      Vector<unsigned> old_domain_for_base_element;
      Vector<unsigned> new_domain_for_base_element;

      // Flat-packed refinement info, labeled by id of locally
      // available root elements
      std::map<unsigned, Vector<unsigned>> flat_packed_refinement_info_for_root;

      // Max. level of refinement
      unsigned max_refinement_level_overall = 0;

      // Prepare the input for the get_data...() method, only copy the
      // data from the structured meshes, TreeBaseMesh meshes
      Vector<unsigned>
        target_domain_for_local_non_halo_element_in_structured_mesh;
      if (n_mesh == 0)
      {
        // Check if the mesh is an structured mesh
        if (!is_unstructured_mesh[0])
        {
          const unsigned nele_mesh =
            target_domain_for_local_non_halo_element.size();
          for (unsigned e = 0; e < nele_mesh; e++)
          {
            const unsigned target_domain =
              target_domain_for_local_non_halo_element[e];
            target_domain_for_local_non_halo_element_in_structured_mesh
              .push_back(target_domain);
          } // for (e < nele_mesh)
        } // if (!is_unstructured_mesh[0])
      } // if (n_mesh == 0)
      else
      {
        // Copy the target domains from the structured meshes only
        for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
        {
          // Check if the mesh is an structured mesh
          if (!is_unstructured_mesh[i_mesh])
          {
            const unsigned nele_sub_mesh =
              target_domain_for_local_non_halo_element_submesh[i_mesh].size();
            for (unsigned e = 0; e < nele_sub_mesh; e++)
            {
              const unsigned target_domain =
                target_domain_for_local_non_halo_element_submesh[i_mesh][e];
              target_domain_for_local_non_halo_element_in_structured_mesh
                .push_back(target_domain);
            } // for (e < nele_sub_mesh)
          } // if (!is_triangle_mesh_base[i_mesh])
        } // for (i_mesh < n_mesh)
      } // else if (n_mesh == 0)

      // Extract data from current problem
      // sorted into data to be sent to various processors after
      // rebuilding the meshes in a load-balanced form
      get_data_to_be_sent_during_load_balancing(
        target_domain_for_local_non_halo_element_in_structured_mesh,
        send_n,
        send_data,
        send_displacement,
        old_domain_for_base_element,
        new_domain_for_base_element,
        max_refinement_level_overall);

      // Extract flat-packed refinement pattern
      get_flat_packed_refinement_pattern_for_load_balancing(
        old_domain_for_base_element,
        new_domain_for_base_element,
        max_refinement_level_overall,
        flat_packed_refinement_info_for_root);

      if (report_stats)
      {
        t_partition = TimingHelpers::timer();
        oomph_info << "CPU for partition calculation for roots: "
                   << t_partition - t_metis << std::endl;
      }


      // Flush and delete old submeshes and null the global mesh
      //--------------------------------------------------------
      // and rebuild the new (not yet distributed, refined etc.) mesh
      //-------------------------------------------------------------
      // that will be distributed in the new, improved way determined
      //-------------------------------------------------------------
      // by METIS
      //---------
      Vector<unsigned> pruned_refinement_level(
        std::max(int(n_old_sub_meshes), 1));
      if (n_mesh == 0)
      {
        pruned_refinement_level[0] = 0;
        TreeBasedRefineableMeshBase* ref_mesh_pt =
          dynamic_cast<TreeBasedRefineableMeshBase*>(old_mesh_pt[0]);
        if (ref_mesh_pt != 0)
        {
          pruned_refinement_level[0] =
            ref_mesh_pt->uniform_refinement_level_when_pruned();
        }

        // If the mesh is an unstructured mesh (TriangleMeshBase mesh)
        // then we should not delete it since the load balance strategy
        // requires the mesh

        // Delete the mesh if it is not an unstructured mesh
        if (!is_unstructured_mesh[0])
        {
          delete old_mesh_pt[0];
          old_mesh_pt[0] = 0;
        } // if (!is_unstructured_mesh[0])
      } // if (n_mesh==0)
      else
      {
        // Loop over the number of old meshes (required to delete the
        // pointers of structured meshes in the old_mesh_pt structure)
        for (unsigned i_mesh = 0; i_mesh < n_old_sub_meshes; i_mesh++)
        {
          pruned_refinement_level[i_mesh] = 0;
          TreeBasedRefineableMeshBase* ref_mesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(old_mesh_pt[i_mesh]);
          if (ref_mesh_pt != 0)
          {
            pruned_refinement_level[i_mesh] =
              ref_mesh_pt->uniform_refinement_level_when_pruned();
          }

          // If the mesh is an unstructured mesh (TriangleMeshBase mesh)
          // then we should NOT delete it since the load balance strategy
          // requires the mesh

          // Delete the mesh if it is not an unstructured mesh
          if (!is_unstructured_mesh[i_mesh])
          {
            delete old_mesh_pt[i_mesh];
            old_mesh_pt[i_mesh] = 0;
          } // if (!is_unstructured_mesh[i_mesh])

        } // for (i_mesh<n_mesh)

        // Empty storage for sub-meshes
        flush_sub_meshes();

        // Flush the storage for nodes and elements in compound mesh
        // (they've already been deleted in the sub-meshes)
        mesh_pt()->flush_element_and_node_storage();

        // Kill
        delete mesh_pt();
        mesh_pt() = 0;
      } // else if (n_mesh==0)

      bool some_mesh_has_been_pruned = false;
      unsigned n = pruned_refinement_level.size();
      for (unsigned i = 0; i < n; i++)
      {
        if (pruned_refinement_level[i] > 0) some_mesh_has_been_pruned = true;
      }

      // (Re-)build the new mesh(es) -- this must get the problem into the
      // state it was in when it was first distributed!
      build_mesh();

      // Has one of the meshes been pruned; if so refine to the
      // common refinement level
      if (some_mesh_has_been_pruned)
      {
        // Do actions before adapt
        actions_before_adapt();

        // Re-assign number of submeshes -- when this was first
        // set, the problem may have had face meshes that have now
        // disappeared.
        n_mesh = nsub_mesh();

        // Now adapt meshes manually
        if (n_mesh == 0)
        {
          TreeBasedRefineableMeshBase* ref_mesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt());
          if (ref_mesh_pt != 0)
          {
            // Get min and max refinement level
            unsigned local_min_ref = 0;
            unsigned local_max_ref = 0;
            ref_mesh_pt->get_refinement_levels(local_min_ref, local_max_ref);

            // Reconcile between processors: If (e.g. following
            // distribution/pruning) the mesh has no elements on this
            // processor) then ignore its contribution to the poll of
            // max/min refinement levels
            int int_local_min_ref = local_min_ref;
            if (ref_mesh_pt->nelement() == 0)
            {
              int_local_min_ref = INT_MAX;
            }
            int int_min_ref = 0;
            MPI_Allreduce(&int_local_min_ref,
                          &int_min_ref,
                          1,
                          MPI_INT,
                          MPI_MIN,
                          Communicator_pt->mpi_comm());

            // Overall min refinement level over all meshes
            unsigned min_ref = unsigned(int_min_ref);

            // Refine as many times as required to get refinement up to
            // uniform refinement level after last prune
            unsigned nref = pruned_refinement_level[0] - min_ref;
            oomph_info << "Refining one-and-only mesh uniformly " << nref
                       << " times\n";
            for (unsigned i = 0; i < nref; i++)
            {
              ref_mesh_pt->refine_uniformly();
            }
          }
        }
        else
        {
          for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
          {
            TreeBasedRefineableMeshBase* ref_mesh_pt =
              dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(i_mesh));
            if (ref_mesh_pt != 0)
            {
              // Get min and max refinement level
              unsigned local_min_ref = 0;
              unsigned local_max_ref = 0;
              ref_mesh_pt->get_refinement_levels(local_min_ref, local_max_ref);

              // Reconcile between processors: If (e.g. following
              // distribution/pruning) the mesh has no elements on this
              // processor) then ignore its contribution to the poll of
              // max/min refinement levels
              int int_local_min_ref = local_min_ref;
              if (ref_mesh_pt->nelement() == 0)
              {
                int_local_min_ref = INT_MAX;
              }
              int int_min_ref = 0;
              MPI_Allreduce(&int_local_min_ref,
                            &int_min_ref,
                            1,
                            MPI_INT,
                            MPI_MIN,
                            Communicator_pt->mpi_comm());

              // Overall min refinement level over all meshes
              unsigned min_ref = unsigned(int_min_ref);

              // Refine as many times as required to get refinement up to
              // uniform refinement level after last prune
              unsigned nref = pruned_refinement_level[i_mesh] - min_ref;
              oomph_info << "Refining sub-mesh " << i_mesh << " uniformly "
                         << nref << " times\n";
              for (unsigned i = 0; i < nref; i++)
              {
                ref_mesh_pt->refine_uniformly();
              }
            }
          }
          // Rebuild the global mesh
          rebuild_global_mesh();
        }

        // Do actions after adapt
        actions_after_adapt();

        // Re-assign number of submeshes -- when this was first
        // set, the problem may have had face meshes that have now
        // disappeared.
        n_mesh = nsub_mesh();
      } // if (some_mesh_has_been_pruned)


      // Perform any actions before distribution but now for the new mesh
      // NOTE: This does NOT replicate the actions_before_distribute()
      // call made above for the previous mesh!
      actions_before_distribute();

      // Do some book-keeping
      //---------------------

      // Re-assign number of submeshes -- when this was first
      // set, the problem may have had face meshes that have now
      // disappeared.
      n_mesh = nsub_mesh();

      // The submeshes, if they exist, need to know their own element
      // domains.
      // NOTE: This vector only stores the target domains or the
      // element partition for structured meshes
      Vector<Vector<unsigned>> submesh_element_partition(n_mesh);
      if (n_mesh != 0)
      {
        unsigned count = 0;
        for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
        {
          // Only work with structured meshes
          if (!is_unstructured_mesh[i_mesh])
          {
            // Get the number of element in the mesh
            const unsigned nsub_ele = mesh_pt(i_mesh)->nelement();
            submesh_element_partition[i_mesh].resize(nsub_ele);
            for (unsigned e = 0; e < nsub_ele; e++)
            {
              submesh_element_partition[i_mesh][e] =
                new_domain_for_base_element[count++];
            } // for (e<nsub_elem)
          } // if (sub_mesh_pt!=0)
        } // for (i_mesh<n_mesh)

#ifdef PARANOID
        const unsigned nnew_domain_for_base_element =
          new_domain_for_base_element.size();
        if (count != nnew_domain_for_base_element)
        {
          std::ostringstream error_stream;
          error_stream
            << "The number of READ target domains for nonhalo elements\n"
            << " is (" << count << "), but the number of target domains for\n"
            << "nonhalo elements is (" << nnew_domain_for_base_element
            << ")!\n";
          throw OomphLibError(error_stream.str(),
                              "Problem::load_balance()",
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      } // if (n_mesh!=0)

      // Setup the map between "root" element and number in global mesh
      // again

      // This map is only established for structured meshes, then we
      // need to check here the type of mesh
      if (n_mesh == 0)
      {
        // Check if the only one mesh is an stuctured mesh
        if (!is_unstructured_mesh[0])
        {
          const unsigned n_ele = mesh_pt()->nelement();
          Base_mesh_element_pt.resize(n_ele);
          Base_mesh_element_number_plus_one.clear();
          for (unsigned e = 0; e < n_ele; e++)
          {
            GeneralisedElement* el_pt = mesh_pt()->element_pt(e);
            Base_mesh_element_number_plus_one[el_pt] = e + 1;
            Base_mesh_element_pt[e] = el_pt;
          } // for (e<n_ele)
        } // if (!is_triangle_mesh_base[0])
      } // if (n_mesh==0)
      else
      {
        // If we have submeshes then we only add those elements that
        // belong to structured meshes, but first compute the number of
        // total elements in the structured sub-meshes
        unsigned nglobal_element = 0;
        for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
        {
          // Check if mesh is an structured mesh
          if (!is_unstructured_mesh[i_mesh])
          {
            nglobal_element += mesh_pt(i_mesh)->nelement();
          } // if (!is_triangle_mesh_base[i_mesh])
        } // for (i_mesh<n_mesh)

        // Once computed the number of elements, then resize the
        // structure
        Base_mesh_element_pt.resize(nglobal_element);
        Base_mesh_element_number_plus_one.clear();
        unsigned counter_base = 0;
        for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
        {
          // Check if mesh is a structured mesh
          if (!is_unstructured_mesh[i_mesh])
          {
            const unsigned n_ele = mesh_pt(i_mesh)->nelement();
            for (unsigned e = 0; e < n_ele; e++)
            {
              GeneralisedElement* el_pt = mesh_pt(i_mesh)->element_pt(e);
              Base_mesh_element_number_plus_one[el_pt] = counter_base + 1;
              Base_mesh_element_pt[counter_base] = el_pt;
              // Inrease the global element number
              counter_base++;
            } // for (e<n_ele)
          } // if (!is_triangle_mesh_base[i_mesh])
        } // for (i_mesh<n_mesh)

#ifdef PARANOID
        if (counter_base != nglobal_element)
        {
          std::ostringstream error_stream;
          error_stream << "The number of global elements (" << nglobal_element
                       << ") is not the same as the number of\nadded elements ("
                       << counter_base << ") to the Base_mesh_element_pt data "
                       << "structure!!!\n\n";
          throw OomphLibError(error_stream.str(),
                              "Problem::load_balance()",
                              OOMPH_EXCEPTION_LOCATION);
        } // if (counter_base != nglobal_element)
#endif // #ifdef PARANOID

      } // else if (n_mesh==0)

      // Storage for the number of face elements in the base mesh --
      // element is identified by number of bulk element and face index
      // so we can reconstruct it if and when the FaceElements have been wiped
      // in actions_before_distribute().
      // NOTE: Not really clear (any more) why this is required. Typically
      //       FaceElements get wiped in actions_before_distribute() so
      //       at this point there shouldn't be any of them left.
      //       This is certainly the case in all our currently existing
      //       test codes. However, I'm too scared to take this out
      //       in case it does matter (we're not insisting that FaceElements
      //       are always removed in actions_before_distribute()...).
      std::map<unsigned, std::map<int, unsigned>> face_element_number;
      unsigned n_element = mesh_pt()->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        FaceElement* face_el_pt =
          dynamic_cast<FaceElement*>(mesh_pt()->finite_element_pt(e));
        if (face_el_pt != 0)
        {
#ifdef PARANOID
          std::stringstream info;
          info << "================================================\n";
          info << "INFO: I've come across a FaceElement while \n";
          info << "       load-balancing a problem. \n";
          info << "================================================\n";
          oomph_info << info.str() << std::endl;
#endif
          FiniteElement* bulk_elem_pt = face_el_pt->bulk_element_pt();
          unsigned e_bulk = Base_mesh_element_number_plus_one[bulk_elem_pt];
#ifdef PARANOID
          if (e_bulk == 0)
          {
            throw OomphLibError("Base_mesh_element_number_plus_one[...]=0",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif
          e_bulk -= 1;
          int face_index = face_el_pt->face_index();
          face_element_number[e_bulk][face_index] = e;
        }
      }

      // Distribute the (sub)meshes
      //---------------------------
      Vector<GeneralisedElement*> deleted_element_pt;
      if (n_mesh == 0)
      {
        // Only distribute (load balance strategy) if this is an
        // structured mesh
        if (!is_unstructured_mesh[0])
        {
#ifdef PARANOID
          if (mesh_pt()->nelement() != new_domain_for_base_element.size())
          {
            std::ostringstream error_stream;
            error_stream << "Distributing one-and-only mesh containing "
                         << mesh_pt()->nelement() << " elements with info for "
                         << new_domain_for_base_element.size() << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif

          if (report_stats)
          {
            oomph_info << "Distributing one and only mesh\n"
                       << "------------------------------" << std::endl;
          }

          // No pre-set distribution from restart that may leave some
          // processors empty so no need to overrule deletion of elements
          bool overrule_keep_as_halo_element_status = false;

          mesh_pt()->distribute(this->communicator_pt(),
                                new_domain_for_base_element,
                                deleted_element_pt,
                                doc_info,
                                report_stats,
                                overrule_keep_as_halo_element_status);

        } // if (!is_unstructured_mesh[0])

      } // if (n_mesh==0)
      else // There are submeshes, "distribute" each one separately
      {
        // Rebuild the mesh only if one of the meshes was modified
        bool need_to_rebuild_mesh = false;
        for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
        {
          // Perform the load balancing based on distribution in the
          // structured meshes only
          if (!is_unstructured_mesh[i_mesh])
          {
            if (report_stats)
            {
              oomph_info << "Distributing submesh " << i_mesh << " of "
                         << n_mesh << " in total\n"
                         << "---------------------------------------------"
                         << std::endl;
            }

            // Set the doc_info number to reflect the submesh
            doc_info.number() = i_mesh;

            // No pre-set distribution from restart that may leave some
            // processors empty so no need to overrule deletion of elements
            bool overrule_keep_as_halo_element_status = false;
            mesh_pt(i_mesh)->distribute(this->communicator_pt(),
                                        submesh_element_partition[i_mesh],
                                        deleted_element_pt,
                                        doc_info,
                                        report_stats,
                                        overrule_keep_as_halo_element_status);

            // Set the flag to rebuild the global mesh
            need_to_rebuild_mesh = true;

          } // if (!is_unstructured_mesh[i_mesh])

        } // for (i_mesh<n_mesh)

        if (need_to_rebuild_mesh)
        {
          // Rebuild the global mesh
          rebuild_global_mesh();
        } // if (need_to_rebuild_mesh)

      } // else if (n_mesh==0)

      // Null out information associated with deleted elements
      unsigned n_del = deleted_element_pt.size();
      for (unsigned e = 0; e < n_del; e++)
      {
        GeneralisedElement* el_pt = deleted_element_pt[e];
        unsigned old_el_number = Base_mesh_element_number_plus_one[el_pt] - 1;
        Base_mesh_element_number_plus_one[el_pt] = 0;
        Base_mesh_element_pt[old_el_number] = 0;
      }

      // Has one of the meshes been pruned before distribution? If so
      // then prune here now
      if (some_mesh_has_been_pruned)
      {
        Vector<GeneralisedElement*> deleted_element_pt;
        if (n_mesh == 0)
        {
          TreeBasedRefineableMeshBase* ref_mesh_pt =
            dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt());
          if (ref_mesh_pt != 0)
          {
            ref_mesh_pt->prune_halo_elements_and_nodes(
              deleted_element_pt, doc_info, report_stats);
          }
        }
        else
        {
          for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
          {
            TreeBasedRefineableMeshBase* ref_mesh_pt =
              dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt(i_mesh));
            if (ref_mesh_pt != 0)
            {
              ref_mesh_pt->prune_halo_elements_and_nodes(
                deleted_element_pt, doc_info, report_stats);
            }
          }
          // Rebuild the global mesh
          rebuild_global_mesh();
        }

        // Null out information associated with deleted elements
        unsigned n_del = deleted_element_pt.size();
        for (unsigned e = 0; e < n_del; e++)
        {
          GeneralisedElement* el_pt = deleted_element_pt[e];
          unsigned old_el_number = Base_mesh_element_number_plus_one[el_pt] - 1;
          Base_mesh_element_number_plus_one[el_pt] = 0;
          Base_mesh_element_pt[old_el_number] = 0;
        }

        // Setup the map between "root" element and number in global mesh again
        setup_base_mesh_info_after_pruning();
      }

      if (report_stats)
      {
        t_distribute = TimingHelpers::timer();
        oomph_info << "CPU for build and distribution of new mesh(es): "
                   << t_distribute - t_partition << std::endl;
      }


      // Send refinement info to other processors
      //-----------------------------------------

      // Storage for refinement pattern:  Given ID of root element,
      // root_element_id, and current refinement level, level, the e-th entry in
      // refinement_info_for_root_elements[root_element_id][level][e] is equal

      // to 2 if the e-th element (using the enumeration when the mesh has been
      // refined to the level-th level) is to be refined during the next
      // refinement; it's 1 if it's not to be refined.
      Vector<Vector<Vector<unsigned>>> refinement_info_for_root_elements;


      // Send refinement information between processors, using flat-packed
      // information accumulated earlier
      send_refinement_info_helper(old_domain_for_base_element,
                                  new_domain_for_base_element,
                                  max_refinement_level_overall,
                                  flat_packed_refinement_info_for_root,
                                  refinement_info_for_root_elements);

      // Refine each mesh based upon refinement information stored for each root
      //------------------------------------------------------------------------
      refine_distributed_base_mesh(refinement_info_for_root_elements,
                                   max_refinement_level_overall);

      if (report_stats)
      {
        t_refine = TimingHelpers::timer();
        oomph_info << "CPU for refinement of base mesh: "
                   << t_refine - t_distribute << std::endl;
      }

      // NOTE: The following two calls are important e.g. when
      //       FaceElements that resize nodes are attached/detached
      //       after/before adaptation. If we don't attach them
      //       on the newly built/refined mesh, there isn't enough
      //       storage for the nodal values that are sent around
      //       (in a flat-packed format) resulting in total disaster.
      //       So we attach them first, but then immediatly strip
      //       them out again because the FaceElements themselves
      //       will have been stripped out before distribution/adaptation.

      // Do actions after adapt because we have just adapted the mesh.
      actions_after_adapt();

      // Now strip it back out to get problem into the same state
      // it was in when data to be sent was recorded.
      actions_before_adapt();

      // Send the stored values in each root from the old mesh into the new mesh
      //------------------------------------------------------------------------
      send_data_to_be_sent_during_load_balancing(
        send_n, send_data, send_displacement);

      // If there are unstructured meshes here we perform the load
      // balancing of those meshes
      if (are_there_unstructured_meshes)
      {
        // Delete any storage of external elements and nodes
        this->delete_all_external_storage();

        if (n_mesh == 0)
        {
          // Before doing the load balancing delete the mesh created at
          // calling build_mesh(), and restore the pointer to the old
          // mesh

          // It MUST be an unstructured mesh, otherwise we should not be
          // here
          if (is_unstructured_mesh[0])
          {
            // Delete the new created mesh
            delete mesh_pt();
            // Re-assign the pointer to the old mesh
            this->mesh_pt() = old_mesh_pt[0];
          } // if (is_unstructured_mesh[0])
#ifdef PARANOID
          else
          {
            std::ostringstream error_stream;
            error_stream << "The only one mesh in the problem is not an "
                            "unstructured mesh,\n"
                         << "but the flag 'are_there_unstructures_meshes' ("
                         << are_there_unstructured_meshes
                         << ") was turned on,\n"
                         << "this is weird. Please check for any  condition "
                            "that may have\n"
                         << "turned on this flag!!!!\n\n";
            throw OomphLibError(error_stream.str(),
                                "Problem::load_balance()",
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif

          unstructured_mesh_pt[0]->load_balance(
            target_domain_for_local_non_halo_element);
        } // if (n_mesh == 0)
        else
        {
          // Before doing the load balancing delete the meshes created
          // at calling build_mesh(), and restore the pointer to the
          // old meshes
          for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
          {
            if (is_unstructured_mesh[i_mesh])
            {
              // Delete the new created mesh
              delete mesh_pt(i_mesh);
              // Now point it to nothing
              mesh_pt(i_mesh) = 0;
              // ... and re-assign the pointer to the old mesh
              this->mesh_pt(i_mesh) = old_mesh_pt[i_mesh];
            } // if (is_unstructured_mesh[i_mesh])

          } // for (i_mesh<n_mesh)

          // Empty storage for sub-meshes
          // flush_sub_meshes();

          // Flush the storage for nodes and elements in compound mesh
          // (they've already been deleted in the sub-meshes)
          mesh_pt()->flush_element_and_node_storage();

          // Now we can procede with the load balancing thing
          for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
          {
            if (is_unstructured_mesh[i_mesh])
            {
              // Get the number of elements in the "i_mesh" (the old one)
              const unsigned n_element = old_mesh_pt[i_mesh]->nelement();

              // Perform the load balancing if there are elements in the
              // mesh. We check for this case because the meshes created
              // from face elements have been cleaned in
              // "actions_before_distribute()"
              if (n_element > 0 && is_unstructured_mesh[i_mesh])
              {
                unstructured_mesh_pt[i_mesh]->load_balance(
                  target_domain_for_local_non_halo_element_submesh[i_mesh]);
              } // if (n_element > 0)
            } // if (is_unstructured_mesh[i_mesh)]
          } // for (i_mesh < n_mesh)

          // Rebuild the global mesh
          rebuild_global_mesh();

        } // else if (n_mesh == 0)

      } // if (are_there_unstructured_meshes)

      if (report_stats)
      {
        t_copy_solution = TimingHelpers::timer();
        oomph_info << "CPU for transferring solution to new mesh(es): "
                   << t_copy_solution - t_refine << std::endl;
        oomph_info << "CPU for load balancing: " << t_copy_solution - t_start
                   << std::endl;
      }

      // Do actions after distribution
      actions_after_distribute();

      // Re-assign equation numbers
#ifdef PARANOID
      unsigned n_dof = assign_eqn_numbers();
#else
      assign_eqn_numbers();
#endif

      if (report_stats)
      {
        oomph_info
          << "Total number of elements on this processor after load balance: "
          << mesh_pt()->nelement() << std::endl;

        oomph_info << "Number of non-halo elements on this processor after "
                      "load balance: "
                   << mesh_pt()->nnon_halo_element() << std::endl;
      }

#ifdef PARANOID
      if (n_dof != old_ndof)
      {
        std::ostringstream error_stream;
        error_stream
          << "Number of dofs in load_balance() has changed from " << old_ndof
          << " to " << n_dof << "\n"
          << "Check that you've implemented any necessary "
             "actions_before/after\n"
          << "adapt/distribute functions, e.g. to pin redundant pressure dofs"
          << " etc.\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    // Finally synchronise all dofs to allow halo check to pass
    synchronise_all_dofs();

    double end_t = TimingHelpers::timer();
    oomph_info << "Time for load_balance() [sec]    : " << end_t - start_t
               << std::endl;
  }


  //==========================================================================
  /// Send refinement information between processors
  //==========================================================================
  void Problem::send_refinement_info_helper(
    Vector<unsigned>& old_domain_for_base_element,
    Vector<unsigned>& new_domain_for_base_element,
    const unsigned& max_refinement_level_overall,
    std::map<unsigned, Vector<unsigned>>& flat_packed_refinement_info_for_root,
    Vector<Vector<Vector<unsigned>>>& refinement_info_for_root_elements)
  {
    // Number of processes etc.
    const int n_proc = this->communicator_pt()->nproc();
    const int my_rank = this->communicator_pt()->my_rank();

    // Make space
    unsigned n_base_element = old_domain_for_base_element.size();
    refinement_info_for_root_elements.resize(n_base_element);

    // Make space for list of domains that the refinement info
    // is to be forwarded to
    std::map<unsigned, Vector<unsigned>> halo_domain_of_haloed_base_element;

    // Find out haloed elements in new, redistributed problem
    //-------------------------------------------------------

    // halo_domains[e][j] = j-th halo domain associated with (haloed) element e
    std::map<unsigned, Vector<unsigned>> halo_domains;

    // Loop over sub meshes
    unsigned n_sub_mesh = nsub_mesh();
    unsigned max_mesh = std::max(n_sub_mesh, unsigned(1));
    for (unsigned i_mesh = 0; i_mesh < max_mesh; i_mesh++)
    {
      // Choose the right mesh
      Mesh* my_mesh_pt = 0;
      if (n_sub_mesh == 0)
      {
        my_mesh_pt = mesh_pt();
      }
      else
      {
        my_mesh_pt = mesh_pt(i_mesh);
      }

      // Only work with structured meshes
      TriangleMeshBase* sub_mesh_pt =
        dynamic_cast<TriangleMeshBase*>(mesh_pt(i_mesh));
      if (!(sub_mesh_pt != 0))
      {
        // Loop over processors to find haloed elements -- need to
        // send their refinement patterns processors that hold their
        // halo counterparts!
        for (int p = 0; p < n_proc; p++)
        {
          Vector<GeneralisedElement*> haloed_elem_pt =
            my_mesh_pt->haloed_element_pt(p);
          unsigned nhaloed = haloed_elem_pt.size();
          for (unsigned h = 0; h < nhaloed; h++)
          {
            // This element must send its refinement information to processor p
            unsigned e = Base_mesh_element_number_plus_one[haloed_elem_pt[h]];
#ifdef PARANOID
            if (e == 0)
            {
              throw OomphLibError("Base_mesh_element_number_plus_one[...]=0",
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
#endif
            e -= 1;
            halo_domains[e].push_back(p);
          }
        }
      } // if (!(sub_mesh_pt!=0))
    } // for (i_mesh<max_mesh)

    // Accumulate relevant flat-packed refinement data to be sent to
    //--------------------------------------------------------------
    // various processors
    //-------------------

    // Map to accumulate unsigned data to be sent to each processor
    // (map for sparsity)
    std::map<unsigned, Vector<unsigned>> data_for_proc;

    // Number of base elements to be sent to specified domain
    Vector<unsigned> nbase_elements_for_proc(n_proc, 0);

    // Total number of entries in send vector
    unsigned count = 0;

    // Loop over all base elements
    //----------------------------
    for (unsigned e = 0; e < n_base_element; e++)
    {
      // Is it one of mine (i.e. was it a non-halo element on this
      //----------------------------------------------------------
      // processor before re-distribution, and do I therefore hold
      //----------------------------------------------------------
      // refinement information for it)?
      //--------------------------------
      if (int(old_domain_for_base_element[e]) == my_rank)
      {
        // Where does it go?
        unsigned new_domain = new_domain_for_base_element[e];

        // Keep counting
        nbase_elements_for_proc[new_domain]++;

        // If it stays local, deal with it here
        if (int(new_domain) == my_rank)
        {
          // Record on which other procs/domains the refinement info for
          // this element is required because it's haloed.
          unsigned nhalo = halo_domains[e].size();
          halo_domain_of_haloed_base_element[e].resize(nhalo);
          for (unsigned j = 0; j < nhalo; j++)
          {
            halo_domain_of_haloed_base_element[e][j] = halo_domains[e][j];
          }

          // Provide storage for refinement pattern
          refinement_info_for_root_elements[e].resize(
            max_refinement_level_overall);

#ifdef PARANOID
          // Get number of additional data sent for check
          unsigned n_additional_data =
            flat_packed_refinement_info_for_root[e].size();
#endif

          // Get number of tree nodes
          unsigned n_tree_nodes = flat_packed_refinement_info_for_root[e][0];

          // Counter for entries to be processed locally
          unsigned local_count = 1; // (have already processed zero-th entry)

          // Loop over levels and number of nodes in tree
          for (unsigned level = 0; level < max_refinement_level_overall;
               level++)
          {
            for (unsigned ee = 0; ee < n_tree_nodes; ee++)
            {
              // Element exists at this level
              if (flat_packed_refinement_info_for_root[e][local_count] == 1)
              {
                local_count++;

                // Element should be refined
                if (flat_packed_refinement_info_for_root[e][local_count] == 1)
                {
                  refinement_info_for_root_elements[e][level].push_back(2);
                  local_count++;
                }
                // Element should not be refined
                else
                {
                  refinement_info_for_root_elements[e][level].push_back(1);
                  local_count++;
                }
              }
              // Element does not exist at this level
              else
              {
                refinement_info_for_root_elements[e][level].push_back(0);
                local_count++;
              }
            }
          }

#ifdef PARANOID
          if (n_additional_data != local_count)
          {
            std::stringstream error_message;
            error_message << "Number of additional data: " << n_additional_data
                          << " doesn't match that actually send: "
                          << local_count << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif
        }
        // Element in question is not one of mine so prepare for sending
        //--------------------------------------------------------------
        else
        {
          // Make space
          unsigned current_size = data_for_proc[new_domain].size();
          unsigned n_additional_data =
            flat_packed_refinement_info_for_root[e].size();
          data_for_proc[new_domain].reserve(current_size + n_additional_data +
                                            2);

          // Keep counting
          count += n_additional_data + 2;

          // Add base element number
          data_for_proc[new_domain].push_back(e);

#ifdef PARANOID
          // Add number of flat-packed instructions to follow
          data_for_proc[new_domain].push_back(n_additional_data);
#endif

          // Add flat packed refinement data
          for (unsigned j = 0; j < n_additional_data; j++)
          {
            data_for_proc[new_domain].push_back(
              flat_packed_refinement_info_for_root[e][j]);
          }
        }
      }
    }


    // Now do the actual send/receive
    //-------------------------------

    // Storage for number of data to be sent to each processor
    Vector<int> send_n(n_proc, 0);

    // Storage for all values to be sent to all processors
    Vector<unsigned> send_data;
    send_data.reserve(count);

    // Start location within send_data for data to be sent to each processor
    Vector<int> send_displacement(n_proc, 0);

    // Loop over all processors
    for (int rank = 0; rank < n_proc; rank++)
    {
      // Set the offset for the current processor
      send_displacement[rank] = send_data.size();

      // Don't bother to do anything if the processor in the loop is the
      // current processor
      if (rank != my_rank)
      {
        // Record how many base elements are to be sent
        send_data.push_back(nbase_elements_for_proc[rank]);

        // Add data
        unsigned n_data = data_for_proc[rank].size();
        for (unsigned j = 0; j < n_data; j++)
        {
          send_data.push_back(data_for_proc[rank][j]);
        }
      }

      // Find the number of data added to the vector
      send_n[rank] = send_data.size() - send_displacement[rank];
    }

    // Storage for the number of data to be received from each processor
    Vector<int> receive_n(n_proc, 0);

    // Now send numbers of data to be sent between all processors
    MPI_Alltoall(&send_n[0],
                 1,
                 MPI_INT,
                 &receive_n[0],
                 1,
                 MPI_INT,
                 this->communicator_pt()->mpi_comm());

    // We now prepare the data to be received
    // by working out the displacements from the received data
    Vector<int> receive_displacement(n_proc, 0);
    int receive_data_count = 0;
    for (int rank = 0; rank < n_proc; ++rank)
    {
      // Displacement is number of data received so far
      receive_displacement[rank] = receive_data_count;
      receive_data_count += receive_n[rank];
    }

    // Now resize the receive buffer for all data from all processors
    // Make sure that it has a size of at least one
    if (receive_data_count == 0)
    {
      ++receive_data_count;
    }
    Vector<unsigned> receive_data(receive_data_count);

    // Make sure that the send buffer has size at least one
    // so that we don't get a segmentation fault
    if (send_data.size() == 0)
    {
      send_data.resize(1);
    }

    // Now send the data between all the processors
    MPI_Alltoallv(&send_data[0],
                  &send_n[0],
                  &send_displacement[0],
                  MPI_UNSIGNED,
                  &receive_data[0],
                  &receive_n[0],
                  &receive_displacement[0],
                  MPI_UNSIGNED,
                  this->communicator_pt()->mpi_comm());


    // Now use the received data to update
    //-----------------------------------
    for (int send_rank = 0; send_rank < n_proc; send_rank++)
    {
      // Don't bother to do anything for the processor corresponding to the
      // current processor or if no data were received from this processor
      if ((send_rank != my_rank) && (receive_n[send_rank] != 0))
      {
        // Counter for the data within the large array
        unsigned count = receive_displacement[send_rank];

        // Loop over base elements
        unsigned nbase_element = receive_data[count];
        count++;
        for (unsigned b = 0; b < nbase_element; b++)
        {
          //  Get base element number
          unsigned base_element_number = receive_data[count];
          count++;

          // Record on which other procs/domains the refinement info for
          // this element is required because it's haloed.
          unsigned nhalo = halo_domains[base_element_number].size();
          halo_domain_of_haloed_base_element[base_element_number].resize(nhalo);
          for (unsigned j = 0; j < nhalo; j++)
          {
            halo_domain_of_haloed_base_element[base_element_number][j] =
              halo_domains[base_element_number][j];
          }

          // Provide storage for refinement pattern
          refinement_info_for_root_elements[base_element_number].resize(
            max_refinement_level_overall);

          // Get number of flat-packed instructions to follow
          // (only used for check)
#ifdef PARANOID
          unsigned n_additional_data = receive_data[count];
          count++;

          // Counter for number of additional data (validation only)
          unsigned check_count = 0;
#endif

          // Get number of tree nodes
          unsigned n_tree_nodes = receive_data[count];
          count++;

#ifdef PARANOID
          check_count++;
#endif

          // Loop over levels and number of nodes in tree
          for (unsigned level = 0; level < max_refinement_level_overall;
               level++)
          {
            for (unsigned e = 0; e < n_tree_nodes; e++)
            {
              // Element exists at this level
              if (receive_data[count] == 1)
              {
                count++;

#ifdef PARANOID
                check_count++;
#endif

                // Element should be refined
                if (receive_data[count] == 1)
                {
                  refinement_info_for_root_elements[base_element_number][level]
                    .push_back(2);
                  count++;

#ifdef PARANOID
                  check_count++;
#endif
                }
                // Element should not be refined
                else
                {
                  refinement_info_for_root_elements[base_element_number][level]
                    .push_back(1);
                  count++;

#ifdef PARANOID
                  check_count++;
#endif
                }
              }
              // Element does not exist at this level
              else
              {
                refinement_info_for_root_elements[base_element_number][level]
                  .push_back(0);
                count++;

#ifdef PARANOID
                check_count++;
#endif
              }
            }
          }

#ifdef PARANOID
          if (n_additional_data != check_count)
          {
            std::stringstream error_message;
            error_message << "Number of additional data: " << n_additional_data
                          << " doesn't match that actually send: "
                          << check_count << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif
        }
      }
    }


    // Now send the fully assembled refinement info to halo elements
    //---------------------------------------------------------------
    {
      // Accumulate data to be sent
      //---------------------------

      // Map to accumulate data to be sent to other procs
      // (map for sparsity)
      std::map<unsigned, Vector<unsigned>> data_for_proc;

      // Number of base elements to be sent to specified domain
      Vector<unsigned> nbase_elements_for_proc(n_proc, 0);

      // Loop over all haloed root elements and find out which
      // processors they have haloes on
      for (std::map<unsigned, Vector<unsigned>>::iterator it =
             halo_domain_of_haloed_base_element.begin();
           it != halo_domain_of_haloed_base_element.end();
           it++)
      {
        // Get base element number
        unsigned base_element_number = (*it).first;

        // Loop over target domains
        Vector<unsigned> domains = (*it).second;
        unsigned nd = domains.size();
        for (unsigned jd = 0; jd < nd; jd++)
        {
          // Actual number of domain
          unsigned d = domains[jd];

          // Keep counting number of base elemements for domain
          nbase_elements_for_proc[d]++;

          // Write base element number
          data_for_proc[d].push_back(base_element_number);

          // Write refinement info in flat-packed form
          for (unsigned level = 0; level < max_refinement_level_overall;
               level++)
          {
            // Number of entries at each level
            unsigned n =
              refinement_info_for_root_elements[base_element_number][level]
                .size();
            data_for_proc[d].push_back(n);
            for (unsigned j = 0; j < n; j++)
            {
              data_for_proc[d].push_back(
                refinement_info_for_root_elements[base_element_number][level]
                                                 [j]);
            }
          }
        }
      }


      // Do the actual send
      //-------------------

      // Storage for number of data to be sent to each processor
      Vector<int> send_n(n_proc, 0);

      // Storage for all values to be sent to all processors
      Vector<unsigned> send_data;
      send_data.reserve(count);

      // Start location within send_data for data to be sent to each processor
      Vector<int> send_displacement(n_proc, 0);

      // Loop over all processors
      for (int rank = 0; rank < n_proc; rank++)
      {
        // Set the offset for the current processor
        send_displacement[rank] = send_data.size();

        // Don't bother to do anything if the processor in the loop is the
        // current processor
        if (rank != my_rank)
        {
          // Record how many base elements are to be sent
          send_data.push_back(nbase_elements_for_proc[rank]);

          // Add data
          unsigned n_data = data_for_proc[rank].size();
          for (unsigned j = 0; j < n_data; j++)
          {
            send_data.push_back(data_for_proc[rank][j]);
          }
        }
        // Find the number of data added to the vector
        send_n[rank] = send_data.size() - send_displacement[rank];
      }

      // Storage for the number of data to be received from each processor
      Vector<int> receive_n(n_proc, 0);

      // Now send numbers of data to be sent between all processors
      MPI_Alltoall(&send_n[0],
                   1,
                   MPI_INT,
                   &receive_n[0],
                   1,
                   MPI_INT,
                   this->communicator_pt()->mpi_comm());

      // We now prepare the data to be received
      // by working out the displacements from the received data
      Vector<int> receive_displacement(n_proc, 0);
      int receive_data_count = 0;
      for (int rank = 0; rank < n_proc; ++rank)
      {
        // Displacement is number of data received so far
        receive_displacement[rank] = receive_data_count;
        receive_data_count += receive_n[rank];
      }

      // Now resize the receive buffer for all data from all processors
      // Make sure that it has a size of at least one
      if (receive_data_count == 0)
      {
        ++receive_data_count;
      }
      Vector<unsigned> receive_data(receive_data_count);

      // Make sure that the send buffer has size at least one
      // so that we don't get a segmentation fault
      if (send_data.size() == 0)
      {
        send_data.resize(1);
      }

      // Now send the data between all the processors
      MPI_Alltoallv(&send_data[0],
                    &send_n[0],
                    &send_displacement[0],
                    MPI_UNSIGNED,
                    &receive_data[0],
                    &receive_n[0],
                    &receive_displacement[0],
                    MPI_UNSIGNED,
                    this->communicator_pt()->mpi_comm());


      // Now use the received data
      //------------------------
      for (int send_rank = 0; send_rank < n_proc; send_rank++)
      {
        // Don't bother to do anything for the processor corresponding to the
        // current processor or if no data were received from this processor
        if ((send_rank != my_rank) && (receive_n[send_rank] != 0))
        {
          // Counter for the data within the large array
          unsigned count = receive_displacement[send_rank];

          // Read number of base elements
          unsigned nbase_element = receive_data[count];
          count++;

          for (unsigned e = 0; e < nbase_element; e++)
          {
            // Read base element number
            unsigned base_element_number = receive_data[count];
            count++;

            // Provide storage for refinement pattern
            refinement_info_for_root_elements[base_element_number].resize(
              max_refinement_level_overall);

            // Read refinement info in flat-packed form
            for (unsigned level = 0; level < max_refinement_level_overall;
                 level++)
            {
              // Read number of entries at each level
              unsigned n = receive_data[count];
              count++;

              // Read entries
              for (unsigned j = 0; j < n; j++)
              {
                refinement_info_for_root_elements[base_element_number][level]
                  .push_back(receive_data[count]);
                count++;
              }
            }
          }
        }
      }
    }
  }

  //==========================================================================
  /// Load balance helper routine: Send data to other
  /// processors during load balancing.
  /// - send_n: Input, number of data to be sent to each processor
  /// - send_data: Input, storage for all values to be sent to all processors
  /// - send_displacement: Input, start location within send_data for data to
  ///   be sent to each processor
  //==========================================================================
  void Problem::send_data_to_be_sent_during_load_balancing(
    Vector<int>& send_n,
    Vector<double>& send_data,
    Vector<int>& send_displacement)
  {
    // Communicator info
    OomphCommunicator* comm_pt = this->communicator_pt();
    const int n_proc = comm_pt->nproc();

    // Storage for the number of data to be received from each processor
    Vector<int> receive_n(n_proc, 0);

    // Now send numbers of data to be sent between all processors
    MPI_Alltoall(&send_n[0],
                 1,
                 MPI_INT,
                 &receive_n[0],
                 1,
                 MPI_INT,
                 this->communicator_pt()->mpi_comm());

    // We now prepare the data to be received
    // by working out the displacements from the received data
    Vector<int> receive_displacement(n_proc, 0);
    int receive_data_count = 0;
    for (int rank = 0; rank < n_proc; ++rank)
    {
      // Displacement is number of data received so far
      receive_displacement[rank] = receive_data_count;
      receive_data_count += receive_n[rank];
    }

    // Now resize the receive buffer for all data from all processors
    // Make sure that it has a size of at least one
    if (receive_data_count == 0)
    {
      ++receive_data_count;
    }
    Vector<double> receive_data(receive_data_count);

    // Make sure that the send buffer has size at least one
    // so that we don't get a segmentation fault
    if (send_data.size() == 0)
    {
      send_data.resize(1);
    }

    // Now send the data between all the processors
    MPI_Alltoallv(&send_data[0],
                  &send_n[0],
                  &send_displacement[0],
                  MPI_DOUBLE,
                  &receive_data[0],
                  &receive_n[0],
                  &receive_displacement[0],
                  MPI_DOUBLE,
                  this->communicator_pt()->mpi_comm());

    unsigned el_count = 0;

    // Only do each node once
    Vector<std::map<Node*, bool>> node_done(n_proc);

    // Now use the received data to update the halo nodes
    for (int send_rank = 0; send_rank < n_proc; send_rank++)
    {
      // Don't bother to do anything if no data were received from this
      // processor
      // NOTE: We do have to loop over our own processor number to process
      //       the data locally.
      if (receive_n[send_rank] != 0)
      {
        // Counter for the data within the large array
        unsigned count = receive_displacement[send_rank];

        // How many batches are there for current rank
        unsigned nbatch = unsigned(receive_data[count]);
        count++;

        // Loop over batches (containing leaves associated with root elements)
        for (unsigned b = 0; b < nbatch; b++)
        {
          // How many elements were received for this batch?
          unsigned nel = unsigned(receive_data[count]);
          count++;

          // Get the unique base/root element number of this batch
          // in unrefined mesh
          unsigned base_el_no = unsigned(receive_data[count]);
          count++;

          // Get pointer to base/root element from reverse lookup scheme
          GeneralisedElement* root_el_pt = Base_mesh_element_pt[base_el_no];

          // Vector for pointers to associated elements in batch
          Vector<GeneralisedElement*> batch_el_pt;

          // Is it a refineable element?
          RefineableElement* ref_root_el_pt =
            dynamic_cast<RefineableElement*>(root_el_pt);
          if (ref_root_el_pt != 0)
          {
            // Get all leaves associated with this base/root element
            Vector<Tree*> all_leaf_nodes_pt;
            ref_root_el_pt->tree_pt()->stick_leaves_into_vector(
              all_leaf_nodes_pt);

            // How many leaves are there?
            unsigned n_leaf = all_leaf_nodes_pt.size();

#ifdef PARANOID
            if (n_leaf != nel)
            {
              std::ostringstream error_message;
              error_message
                << "Number of leaves: " << n_leaf << " "
                << " doesn't match number of elements sent in batch: " << nel
                << "\n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
#endif

            // Loop over batch of elements associated with this base/root
            // element
            batch_el_pt.resize(n_leaf);
            for (unsigned e = 0; e < n_leaf; e++)
            {
              batch_el_pt[e] = all_leaf_nodes_pt[e]->object_pt();
            }
          }
          // Not refineable -- the batch contains just the root element itself
          else
          {
#ifdef PARANOID
            if (1 != nel)
            {
              std::ostringstream error_message;
              error_message
                << "Non-refineable root element should only be associated with"
                << " one element but nel=" << nel << "\n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
#endif
            batch_el_pt.push_back(root_el_pt);
          }

          // Now loop  over all elements in batch
          for (unsigned e = 0; e < nel; e++)
          {
            GeneralisedElement* el_pt = batch_el_pt[e];
            el_count++;

            // FE?
            FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(el_pt);
            if (fe_pt != 0)
            {
              // Loop over nodes
              unsigned nnod = fe_pt->nnode();
              for (unsigned j = 0; j < nnod; j++)
              {
                Node* nod_pt = fe_pt->node_pt(j);
                if (!node_done[send_rank][nod_pt])
                {
                  node_done[send_rank][nod_pt] = true;


                  // Read number of values (as double) to allow for resizing
                  // before read (req'd in case we store data that
                  // got introduced by attaching FaceElements to bulk)
                  unsigned nval = unsigned(receive_data[count]);
                  count++;

#ifdef PARANOID
                  // Does the size match?
                  if (nval < nod_pt->nvalue())
                  {
                    std::ostringstream error_message;
                    error_message
                      << "Node has more values, namely " << nod_pt->nvalue()
                      << ", than we're about to receive, namely " << nval
                      << ". Something's wrong!\n";
                    throw OomphLibError(error_message.str(),
                                        OOMPH_CURRENT_FUNCTION,
                                        OOMPH_EXCEPTION_LOCATION);
                  }
#endif


#ifdef PARANOID
                  // Check if it's been sent as a boundary node
                  unsigned is_boundary_node = unsigned(receive_data[count]);
                  count++;
#endif

                  // Check if it's actually a boundary node
                  BoundaryNodeBase* bnod_pt =
                    dynamic_cast<BoundaryNodeBase*>(nod_pt);
                  if (bnod_pt != 0)
                  {
#ifdef PARANOID
                    // Check if local and received status are consistent
                    if (is_boundary_node != 1)
                    {
                      std::ostringstream error_message;
                      error_message << "Local node is boundary node but "
                                       "information sent is\n"
                                    << "for non-boundary node\n";
                      throw OomphLibError(error_message.str(),
                                          OOMPH_CURRENT_FUNCTION,
                                          OOMPH_EXCEPTION_LOCATION);
                    }
#endif

                    // Do we have entries in the map?
                    unsigned n_entry = unsigned(receive_data[count]);
                    count++;
                    if (n_entry > 0)
                    {
                      // Create storage, if it doesn't already exist, for the
                      // map that will contain the position of the first entry
                      // of this face element's additional values,
                      if (
                        bnod_pt
                          ->index_of_first_value_assigned_by_face_element_pt() ==
                        0)
                      {
                        bnod_pt
                          ->index_of_first_value_assigned_by_face_element_pt() =
                          new std::map<unsigned, unsigned>;
                      }

                      // Get pointer to the map of indices associated with
                      // additional values created by face elements
                      std::map<unsigned, unsigned>* map_pt =
                        bnod_pt
                          ->index_of_first_value_assigned_by_face_element_pt();

                      // Loop over number of entries in map
                      for (unsigned i = 0; i < n_entry; i++)
                      {
                        // Read out pairs...
                        unsigned first = unsigned(receive_data[count]);
                        count++;
                        unsigned second = unsigned(receive_data[count]);
                        count++;

                        // ...and assign
                        (*map_pt)[first] = second;
                      }
                    }
                  }
#ifdef PARANOID
                  // Not a boundary node
                  else
                  {
                    // Check if local and received status are consistent
                    if (is_boundary_node != 0)
                    {
                      std::ostringstream error_message;
                      error_message << "Local node is not a boundary node but "
                                       "information \n"
                                    << "sent is for boundary node.\n";
                      throw OomphLibError(error_message.str(),
                                          OOMPH_CURRENT_FUNCTION,
                                          OOMPH_EXCEPTION_LOCATION);
                    }
                  }
#endif

                  // Do we have to resize? This can happen if node was
                  // resized (due to a FaceElement that hasn't been attached
                  // yet here) when the send data was written. If so make space
                  // for the data here
                  if (nval > nod_pt->nvalue())
                  {
                    nod_pt->resize(nval);
                  }

                  // Now read the actual values
                  nod_pt->read_values_from_vector(receive_data, count);
                }
              }
            }

            // Now add internal data
            el_pt->read_internal_data_values_from_vector(receive_data, count);
          }
        }
      }
    }

    // Now that this is done, we need to synchronise dofs to get
    // the halo element and node values correct
    bool do_halos = true;
    bool do_external_halos = false;
    this->synchronise_dofs(do_halos, do_external_halos);

    // Now rebuild global mesh if required
    unsigned n_mesh = nsub_mesh();
    if (n_mesh != 0)
    {
      bool do_halos = false;
      bool do_external_halos = true;
      this->synchronise_dofs(do_halos, do_external_halos);
      rebuild_global_mesh();
    }
  }


  //==========================================================================
  /// Load balance helper routine: Get data to be sent to other
  /// processors during load balancing and other information about
  /// re-distribution.
  /// - target_domain_for_local_non_halo_element: Input, generated by METIS.
  ///   target_domain_for_local_non_halo_element[e] contains the number
  ///   of the domain [0,1,...,nproc-1] to which non-halo element e on THE
  ///   CURRENT PROCESSOR ONLY has been assigned. The order of the non-halo
  ///   elements is the same as in the Problem's mesh, with the halo
  ///   elements being skipped.
  /// - send_n: Output, number of data to be sent to each processor
  /// - send_data: Output, storage for all values to be sent to all processors
  /// - send_displacement: Output, start location within send_data for data to
  ///   be sent to each processor
  /// - max_refinement_level_overall: Output, max. refinement level of any
  ///   element
  //==========================================================================
  void Problem::get_data_to_be_sent_during_load_balancing(
    const Vector<unsigned>& target_domain_for_local_non_halo_element,
    Vector<int>& send_n,
    Vector<double>& send_data,
    Vector<int>& send_displacement,
    Vector<unsigned>& old_domain_for_base_element,
    Vector<unsigned>& new_domain_for_base_element,
    unsigned& max_refinement_level_overall)
  {
    // Communicator info
    OomphCommunicator* comm_pt = this->communicator_pt();
    const int n_proc = comm_pt->nproc();
    const int my_rank = this->communicator_pt()->my_rank();

    //------------------------------------------------------------------------
    // Overall strategy: Loop over all elements (in structured meshes),
    // identify their corresponding root elements and move all associated
    // leaves together, collecting the leaves in batches.
    // ------------------------------------------------------------------------

    // Map to store whether the root element has been visited yet
    std::map<RefineableElement*, bool> root_el_done;

#ifdef PARANOID

    // Map for checking if all elements associated with same root
    // have the same target processor
    std::map<RefineableElement*, unsigned> target_plus_one_for_root;

#endif

    // Storage for maximum refinement level
    unsigned max_refinement_level = 0;

    // Storage for (vector of) elements associated with target domain
    // (stored in map for sparsity): element_for_processor[d][e] is pointer
    // to e-th element that's supposed to move onto processor (domain) d.
    std::map<unsigned, Vector<GeneralisedElement*>> element_for_processor;

    // Storage for the number of elements in a specified batch of leaf
    // elements, all of which are associated with the same root/base element:
    // nelement_batch_for_processor[d][j] is the number of (leaf)
    // elements (all associated with the same root) to be moved together to
    // domain/processor d, in the j-th batch of elements.
    std::map<unsigned, Vector<unsigned>> nelement_batch_for_processor;

    // Storage for the unique number of the root element (in the unrefined
    // base mesh) whose leaves are moved together in a batch:
    // base_element_for_element_batch_for_processo[d][j] is the number of
    // unique number of the root element (in the unrefined
    // base mesh) of all leaf elements (associated with that root),
    // to be moved together to domain/processor d, in the j-th batch of
    // elements.
    std::map<unsigned, Vector<unsigned>>
      base_element_for_element_batch_for_processor;

    // Record old and new domains for non-halo root elements (will be
    // communicated globally). Initialise to -1 so we can use max
    // to extract the right one via MPI_Allreduce.
    // NOTE: We communicate these globally to facilitate distribution
    //       of refinement pattern. While the data itself can be
    //       sent point-to-point for non-halo elements,
    //       mesh refinement information also needs to be sent for
    //       halo elements which aren't known yet.
    unsigned n_base_element = Base_mesh_element_pt.size();
    Vector<int> old_domain_for_base_element_local(n_base_element, -1);
    Vector<int> new_domain_for_base_element_local(n_base_element, -1);

    // Loop over all non-halo elements on current processor and identify roots
    // -------------------------------------------------------------------
    // All leaf elements in associated tree (must!) get moved together
    //----------------------------------------------------------------
    unsigned count_non_halo_el = 0;
    // Get the number of submeshs, if there are no submeshes, then
    // increase the counter so that the loop below also work for the only
    // one mesh in the problem
    unsigned n_mesh = nsub_mesh();
    if (n_mesh == 0)
    {
      n_mesh = 1;
    }
    // We need to know if there are structure meshes (with elements) as
    // part of the problem in order to perform (or not) the proper
    // communications
    bool are_there_structured_meshes = false;
    // Go for the nonhalo elements only in the TreeBaseMeshes
    for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
    {
      // Only work with structured meshes
      TriangleMeshBase* sub_mesh_pt =
        dynamic_cast<TriangleMeshBase*>(mesh_pt(i_mesh));
      if (!(sub_mesh_pt != 0))
      {
        const unsigned nele = mesh_pt(i_mesh)->nelement();
        if (nele > 0)
        {
          // Change the flag to indicate that there are structured meshes
          // (with elements, because we may have meshes with face
          // elements and therefore zero elements at this point)
          are_there_structured_meshes = true;
        }

        for (unsigned e = 0; e < nele; e++)
        {
          GeneralisedElement* el_pt = mesh_pt(i_mesh)->element_pt(e);
          if (!el_pt->is_halo())
          {
            // New non-halo: Where is this element supposed to go to?
            //-------------------------------------------------------
            unsigned target_domain =
              target_domain_for_local_non_halo_element[count_non_halo_el];

            // Bump up counter for non-halo elements
            count_non_halo_el++;

            // Is it a root element? (It is, trivially, if it's not refineable)
            //------------------------------------------------------------------
            RefineableElement* ref_el_pt =
              dynamic_cast<RefineableElement*>(el_pt);
            if (ref_el_pt == 0)
            {
              // Not refineable so add element itself
              element_for_processor[target_domain].push_back(el_pt);

              // Number of elements associated with this root/base
              // element (just the element itself)
              nelement_batch_for_processor[target_domain].push_back(1);

              // This is the unique base/root element number in unrefined mesh
              unsigned element_number_in_base_mesh =
                Base_mesh_element_number_plus_one[el_pt];
#ifdef PARANOID
              if (element_number_in_base_mesh == 0)
              {
                throw OomphLibError("Base_mesh_element_number_plus_one[...]=0",
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif
              element_number_in_base_mesh -= 1;
              base_element_for_element_batch_for_processor[target_domain]
                .push_back(element_number_in_base_mesh);

              /// Where do I come from, where do I go to?
              old_domain_for_base_element_local[element_number_in_base_mesh] =
                my_rank;
              new_domain_for_base_element_local[element_number_in_base_mesh] =
                target_domain;
            } // if (ref_el_pt==0)
            // It's not a root element so we package its leaves into a batch
            //--------------------------------------------------------------
            // of elements
            //------------
            else
            {
              // Get the root element
              RefineableElement* root_el_pt = ref_el_pt->root_element_pt();

              // Has this root been visited yet?
              if (!root_el_done[root_el_pt])
              {
                // Now we've done it
                root_el_done[root_el_pt] = true;

                // Unique number of root element in base mesh
                unsigned element_number_in_base_mesh =
                  Base_mesh_element_number_plus_one[root_el_pt];
#ifdef PARANOID
                if (element_number_in_base_mesh == 0)
                {
                  throw OomphLibError(
                    "Base_mesh_element_number_plus_one[...]=0",
                    OOMPH_CURRENT_FUNCTION,
                    OOMPH_EXCEPTION_LOCATION);
                }
#endif
                element_number_in_base_mesh -= 1;

                /// Where do I come from, where do I go to?
                old_domain_for_base_element_local[element_number_in_base_mesh] =
                  my_rank;
                new_domain_for_base_element_local[element_number_in_base_mesh] =
                  target_domain;

#ifdef PARANOID
                // Store target domain associated with this root element
                // (offset by one) to allow checking that all elements
                // with the same root move to the same processor
                target_plus_one_for_root[root_el_pt] = target_domain + 1;
#endif

                // Package all leaves into batch of elements
                Vector<Tree*> all_leaf_nodes_pt;
                root_el_pt->tree_pt()->stick_leaves_into_vector(
                  all_leaf_nodes_pt);

                // Number of leaves
                unsigned n_leaf = all_leaf_nodes_pt.size();

                // Number of elements associated with this root/base element
                // (all the leaves)
                nelement_batch_for_processor[target_domain].push_back(n_leaf);

                // Store the unique base/root element number in unrefined mesh
                base_element_for_element_batch_for_processor[target_domain]
                  .push_back(element_number_in_base_mesh);

                // Loop over leaves
                for (unsigned i_leaf = 0; i_leaf < n_leaf; i_leaf++)
                {
                  // Add element object at leaf
                  RefineableElement* leaf_el_pt =
                    all_leaf_nodes_pt[i_leaf]->object_pt();
                  element_for_processor[target_domain].push_back(leaf_el_pt);

                  // Monitor/update maximum refinement level
                  unsigned level = all_leaf_nodes_pt[i_leaf]->level();
                  if (level > max_refinement_level)
                  {
                    max_refinement_level = level;
                  }
                }
              }

#ifdef PARANOID
              // Root element has already been visited
              else
              {
                // We don't have to do anything with this element since it's
                // already been processed earlier, but check that it's scheduled
                // to go onto the same processor as its root.
                if ((target_plus_one_for_root[root_el_pt] - 1) != target_domain)
                {
                  std::ostringstream error_message;
                  error_message
                    << "All elements associated with same root must have "
                    << "same target. during load balancing\n";
                  throw OomphLibError(error_message.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
              }
#endif
            } // else if (ref_el_pt==0)
          } // if (!ele_pt->is_halo())
        } // for (e < nele)
      } // if (!(sub_mesh_pt!=0))
    } // for (i_mesh < n_mesh)

#ifdef PARANOID
    // Have we processed all target domains?
    if (target_domain_for_local_non_halo_element.size() != count_non_halo_el)
    {
      std::ostringstream error_message;
      error_message
        << "Have processed " << count_non_halo_el << " of "
        << target_domain_for_local_non_halo_element.size()
        << " target domains for local non-halo elelemts. \n "
        << "Very Odd -- we do (now) strip out the information for elements\n"
        << "that are removed in actions_before_distribute()...\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Determine max. refinement level and origin/destination scheme
    // -------------------------------------------------------------
    // for all root/base elements
    // --------------------------

    // Allreduce to work out max max refinement level across all processors
    max_refinement_level_overall = 0;

    // Only perform this communications if necessary (it means if there
    // are structured meshes as part of the problem)
    if (are_there_structured_meshes)
    {
      MPI_Allreduce(&max_refinement_level,
                    &max_refinement_level_overall,
                    1,
                    MPI_UNSIGNED,
                    MPI_MAX,
                    comm_pt->mpi_comm());
    } // if (are_there_structured_meshes)

    // Allreduce to tell everybody about the original and new domains
    // for root elements
    Vector<int> tmp_old_domain_for_base_element(n_base_element);

    // Only perform this communications if necessary (it means if there
    // are structured meshes as part of the problem)
    if (are_there_structured_meshes)
    {
      MPI_Allreduce(&old_domain_for_base_element_local[0],
                    &tmp_old_domain_for_base_element[0],
                    n_base_element,
                    MPI_INT,
                    MPI_MAX,
                    comm_pt->mpi_comm());
    } // if (are_there_structured_meshes)

    Vector<int> tmp_new_domain_for_base_element(n_base_element);
    // Only perform this communications if necessary (it means if there
    // are structured meshes as part of the problem)
    if (are_there_structured_meshes)
    {
      MPI_Allreduce(&new_domain_for_base_element_local[0],
                    &tmp_new_domain_for_base_element[0],
                    n_base_element,
                    MPI_INT,
                    MPI_MAX,
                    comm_pt->mpi_comm());
    } // if (are_there_structured_meshes)

    // Copy across (after optional sanity check)
    old_domain_for_base_element.resize(n_base_element);
    new_domain_for_base_element.resize(n_base_element);
    for (unsigned j = 0; j < n_base_element; j++)
    {
#ifdef PARANOID
      if (tmp_old_domain_for_base_element[j] == -1)
      {
        std::ostringstream error_message;
        error_message << "Old domain for base element " << j << ": "
                      << Base_mesh_element_pt[j]
                      << "or its incarnation as refineable el: "
                      << dynamic_cast<RefineableElement*>(
                           Base_mesh_element_pt[j])
                      << " which is of type "
                      << typeid(*Base_mesh_element_pt[j]).name()
                      << " does not\n"
                      << "appear to have been assigned by any processor\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      old_domain_for_base_element[j] = tmp_old_domain_for_base_element[j];
#ifdef PARANOID
      if (tmp_new_domain_for_base_element[j] == -1)
      {
        std::ostringstream error_message;
        error_message << "New domain for base element " << j
                      << "which is of type "
                      << typeid(*Base_mesh_element_pt[j]).name()
                      << " does not\n"
                      << "appear to have been assigned by any processor\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      new_domain_for_base_element[j] = tmp_new_domain_for_base_element[j];
    }


    // Loop over all processors and accumulate data to be sent
    //--------------------------------------------------------
    send_data.clear();

    // Only do each node once (per processor!)
    Vector<std::map<Node*, bool>> node_done(n_proc);

    // Loop over all processors. NOTE: We include current processor
    // since we have to refine local elements too -- store their data
    // in same data structure as the one used for off-processor elements.
    for (int rank = 0; rank < n_proc; rank++)
    {
      // Set the offset for the current processor
      send_displacement[rank] = send_data.size();

#ifdef PARANOID
      // Check that total number of elements processed matches those
      // in individual batches
      unsigned total_nel = element_for_processor[rank].size();
#endif

      // Counter for number of elements
      unsigned el_count = 0;

      // How many baches are there for current rank?
      unsigned nbatch = nelement_batch_for_processor[rank].size();

      // Add to vector of doubles to save on number of comms
      send_data.push_back(double(nbatch));

      // Loop over batches of elemnts associated with same root
      for (unsigned b = 0; b < nbatch; b++)
      {
        // How many elements are to be sent in this batch?
        unsigned nel = nelement_batch_for_processor[rank][b];

        // Get the unique number of the root element in unrefined mesh for
        // all the elements in this batch
        unsigned base_el_no =
          base_element_for_element_batch_for_processor[rank][b];

        // Add unsigneds to send data to minimise number of
        // communications
        send_data.push_back(double(nel));
        send_data.push_back(double(base_el_no));

        // Loop over batch of elements
        for (unsigned e = 0; e < nel; e++)
        {
          // Get element
          GeneralisedElement* el_pt = element_for_processor[rank][el_count];

          // FE?
          FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(el_pt);
          if (fe_pt != 0)
          {
            // Loop over nodes
            unsigned nnod = fe_pt->nnode();
            for (unsigned j = 0; j < nnod; j++)
            {
              Node* nod_pt = fe_pt->node_pt(j);

              // Reconstruct the nodal values/position from the node's
              // possible hanging node representation to be on the safe side
              unsigned n_value = nod_pt->nvalue();
              unsigned nt = nod_pt->ntstorage();
              Vector<double> values(n_value);
              unsigned n_dim = nod_pt->ndim();
              Vector<double> position(n_dim);

              // Loop over all history values
              for (unsigned t = 0; t < nt; t++)
              {
                nod_pt->value(t, values);
                for (unsigned i = 0; i < n_value; i++)
                {
                  nod_pt->set_value(t, i, values[i]);
                }
                nod_pt->position(t, position);
                for (unsigned i = 0; i < n_dim; i++)
                {
                  nod_pt->x(t, i) = position[i];
                }
              }


              // Has the node already been done for current rank?
              if (!node_done[rank][nod_pt])
              {
                // Now it has been done
                node_done[rank][nod_pt] = true;

                // Store number of values (as double) to allow for resizing
                // before read (req'd in case we store data that
                // got introduced by attaching FaceElements to bulk)
                send_data.push_back(double(n_value));

                // Check if it's a boundary node
                BoundaryNodeBase* bnod_pt =
                  dynamic_cast<BoundaryNodeBase*>(nod_pt);

                // Not a boundary node
                if (bnod_pt == 0)
                {
#ifdef PARANOID
                  // Record status for checking
                  send_data.push_back(double(0));
#endif
                }
                // Yes it's a boundary node
                else
                {
#ifdef PARANOID
                  // Record status for checking
                  send_data.push_back(double(1));
#endif
                  // Get pointer to the map of indices associated with
                  // additional values created by face elements
                  std::map<unsigned, unsigned>* map_pt =
                    bnod_pt->index_of_first_value_assigned_by_face_element_pt();

                  // No additional values created
                  if (map_pt == 0)
                  {
                    send_data.push_back(double(0));
                  }
                  // Created additional values
                  else
                  {
                    // How many?
                    send_data.push_back(double(map_pt->size()));

                    // Loop over entries in map and add to send data
                    for (std::map<unsigned, unsigned>::iterator p =
                           map_pt->begin();
                         p != map_pt->end();
                         p++)
                    {
                      send_data.push_back(double((*p).first));
                      send_data.push_back(double((*p).second));
                    }
                  }
                }

                // Add the actual values
                nod_pt->add_values_to_vector(send_data);
              }
            }
          }

          // Now add internal data
          el_pt->add_internal_data_values_to_vector(send_data);

          // Bump up counter in long vector of elements
          el_count++;
        }
      }


#ifdef PARANOID
      // Check that total number of elements matches the total of those
      // in batches
      if (total_nel != el_count)
      {
        std::ostringstream error_message;
        error_message
          << "total_nel: " << total_nel << " "
          << " doesn't match total number of elements sent in batch: "
          << el_count << "\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Find the number of data added to the vector
      send_n[rank] = send_data.size() - send_displacement[rank];
    }
  }


  //==========================================================================
  /// Get flat-packed refinement pattern for each root element in current
  /// mesh (labeled by unique number of root element in unrefined base mesh).
  /// The vector stored for each root element contains the following
  /// information:
  /// - First entry: Number of tree nodes (not just leaves!) in refinement
  ///   tree emanating from this root [Zero if root element is not refineable]
  /// - Loop over all refinement levels
  ///   - Loop over all tree nodes (not just leaves!)
  ///     - If associated element exists when the mesh has been refined to
  ///       this level (either because it has been refined to this level or
  ///       because it's less refined): 1
  ///       - If the element is to be refined: 1; else: 0
  ///     - else (element doesn't exist when mesh is refined to this level
  ///       (because it's more refined): 0
  ///     .
  ///   .
  /// .
  //==========================================================================
  void Problem::get_flat_packed_refinement_pattern_for_load_balancing(
    const Vector<unsigned>& old_domain_for_base_element,
    const Vector<unsigned>& new_domain_for_base_element,
    const unsigned& max_refinement_level_overall,
    std::map<unsigned, Vector<unsigned>>& flat_packed_refinement_info_for_root)
  {
    // Map to store whether the root element has been visited yet
    std::map<RefineableElement*, bool> root_el_done;

    // Get the number of submeshs, if there are no submeshes, then
    // increase the counter so that the loop below also work for the only
    // one mesh in the problem
    unsigned n_mesh = nsub_mesh();
    if (n_mesh == 0)
    {
      n_mesh = 1;
    }
    // Go for the nonhalo elements only in the TreeBaseMeshes
    for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
    {
      // Only work with structured
      TriangleMeshBase* sub_mesh_pt =
        dynamic_cast<TriangleMeshBase*>(mesh_pt(i_mesh));
      if (!(sub_mesh_pt != 0))
      {
        const unsigned nele_submesh = mesh_pt(i_mesh)->nelement();
        for (unsigned e = 0; e < nele_submesh; e++)
        {
          // Get pointer to element
          GeneralisedElement* el_pt = mesh_pt(i_mesh)->element_pt(e);

          // Ignore halos
          if (!el_pt->is_halo())
          {
            // Is it refineable? No!
            RefineableElement* ref_el_pt =
              dynamic_cast<RefineableElement*>(el_pt);
            if (ref_el_pt == 0)
            {
              // The element is not refineable - stick a zero in refinement_info
              // indicating that there are no tree nodes following
              unsigned e = Base_mesh_element_number_plus_one[el_pt];
#ifdef PARANOID
              if (e == 0)
              {
                throw OomphLibError("Base_mesh_element_number_plus_one[...]=0",
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif
              e -= 1;
              flat_packed_refinement_info_for_root[e].push_back(0);
            }
            // Refineable
            else
            {
              // Get the root element
              RefineableElement* root_el_pt = ref_el_pt->root_element_pt();

              // Has this root been visited yet?
              if (!root_el_done[root_el_pt])
              {
                // Get unique number of root element in base mesh
                unsigned root_element_number =
                  Base_mesh_element_number_plus_one[root_el_pt];

#ifdef PARANOID
                if (root_element_number == 0)
                {
                  throw OomphLibError(
                    "Base_mesh_element_number_plus_one[...]=0",
                    OOMPH_CURRENT_FUNCTION,
                    OOMPH_EXCEPTION_LOCATION);
                }
#endif
                root_element_number -= 1;

                // Get all the nodes associated with this root element
                Vector<Tree*> all_tree_nodes_pt;
                root_el_pt->tree_pt()->stick_all_tree_nodes_into_vector(
                  all_tree_nodes_pt);

                // How many tree nodes are there?
                unsigned n_tree_nodes = all_tree_nodes_pt.size();
                flat_packed_refinement_info_for_root[root_element_number]
                  .push_back(n_tree_nodes);

                // Loop over all levels
                for (unsigned current_level = 0;
                     current_level < max_refinement_level_overall;
                     current_level++)
                {
                  // Loop over all tree nodes
                  for (unsigned e = 0; e < n_tree_nodes; e++)
                  {
                    // What's the level of this tree node?
                    unsigned level = all_tree_nodes_pt[e]->level();

                    // Element exists at this refinement level of the mesh
                    // if it's at this level or it's at a lower level and a leaf
                    if ((level == current_level) ||
                        ((level < current_level) &&
                         (all_tree_nodes_pt[e]->is_leaf())))
                    {
                      flat_packed_refinement_info_for_root[root_element_number]
                        .push_back(1);

                      // If it's at this level, and not a leaf, then it will
                      // need to be refined in the new mesh
                      if ((level == current_level) &&
                          (!all_tree_nodes_pt[e]->is_leaf()))
                      {
                        flat_packed_refinement_info_for_root
                          [root_element_number]
                            .push_back(1);
                      }
                      // Element exists at this level and is a leaf so it
                      // doesn't have to be refined
                      else
                      {
                        flat_packed_refinement_info_for_root
                          [root_element_number]
                            .push_back(0);
                      }
                    }
                    // Element does not exist at this level so it doesn't have
                    // to be refined
                    else
                    {
                      flat_packed_refinement_info_for_root[root_element_number]
                        .push_back(0);
                    }
                  }
                }
                // Now we've done it
                root_el_done[root_el_pt] = true;
              }
            }

          } // if (!el_pt->is_halo())
        } // for (e < nele_submesh)
      } // if (!(sub_mesh_pt!=0))
    } // for (i_mesh < n_mesh)
  }

  //==========================================================================
  /// Load balance helper routine:  Function performs max_level_overall
  /// successive refinements of the problem's mesh(es) using the following
  /// procdure: Given ID of root element, root_element_id, and current
  /// refinement level, level, the e-th entry in
  /// refinement_info_for_root_elements[root_element_id][level][e] is equal
  /// to 2 if the e-th element (using the enumeration when the mesh has been
  /// refined to the level-th level) is to be refined during the next
  /// refinement; it's 1 if it's not to be refined.
  //==========================================================================
  void Problem::refine_distributed_base_mesh(
    Vector<Vector<Vector<unsigned>>>& refinement_info_for_root_elements,
    const unsigned& max_level_overall)
  {
    // Loop over sub meshes
    unsigned n_sub_mesh = nsub_mesh();
    unsigned max_mesh = std::max(n_sub_mesh, unsigned(1));
    for (unsigned i_mesh = 0; i_mesh < max_mesh; i_mesh++)
    {
      // Choose the right mesh
      Mesh* my_mesh_pt = 0;
      if (n_sub_mesh == 0)
      {
        my_mesh_pt = mesh_pt();
      }
      else
      {
        my_mesh_pt = mesh_pt(i_mesh);
      }

      // Number of elements on this processor -- currently all elements
      // are "base" elements since the mesh hasn't been refined.
      unsigned n_el_on_this_proc = my_mesh_pt->nelement();

      // Storage for actual refinement pattern:
      // to_be_refined_on_this_proc[level][e] contains the element number
      // of the e-th element that is to refined at the level-th refinement level
      Vector<Vector<unsigned>> to_be_refined_on_this_proc(max_level_overall);

      // Count, at each level, the total number of elements in the mesh
      // (we can accumulate this because we know that elements are
      // enumerated tree by tree).
      Vector<unsigned> el_count_on_this_proc(max_level_overall, 0);

      // Loop over levels where refinement is taking place
      for (unsigned level = 0; level < max_level_overall; level++)
      {
        // Loop over roots = unrefined elements on this processor in order.
        // Note that this loops over the trees in unique order
        for (unsigned e = 0; e < n_el_on_this_proc; e++)
        {
          // Get the (root) element
          FiniteElement* el_pt = my_mesh_pt->finite_element_pt(e);

          // What is its unique number in the base mesh
          unsigned root_el_no = Base_mesh_element_number_plus_one[el_pt];
#ifdef PARANOID
          if (root_el_no == 0)
          {
            throw OomphLibError("Base_mesh_element_number_plus_one[...]=0",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif
          root_el_no -= 1;

          // Number of refinements to be performed starting from current
          // root element
          unsigned n_refinements =
            refinement_info_for_root_elements[root_el_no].size();

          // Perform refinement?
          if (level < n_refinements)
          {
            // Loop over elements at this level
            unsigned n_el =
              refinement_info_for_root_elements[root_el_no][level].size();
            for (unsigned ee = 0; ee < n_el; ee++)
            {
              // Refinement code 2: Element is to be refined at this
              // level
              if (refinement_info_for_root_elements[root_el_no][level][ee] == 2)
              {
                to_be_refined_on_this_proc[level].push_back(
                  el_count_on_this_proc[level]);
                el_count_on_this_proc[level]++;
              }
              // Refinement code 1: Element should not be refined at this
              // level -- keep going
              else if (refinement_info_for_root_elements[root_el_no][level]
                                                        [ee] == 1)
              {
                el_count_on_this_proc[level]++;
              }
            }
          }

        } // end of loop over elements on proc; all of which should be root
      }

      // Now do the actual refinement
      TreeBasedRefineableMeshBase* ref_mesh_pt =
        dynamic_cast<TreeBasedRefineableMeshBase*>(my_mesh_pt);
      if (ref_mesh_pt != 0)
      {
        ref_mesh_pt->refine_base_mesh(to_be_refined_on_this_proc);
      }
    }

    // Rebuild global mesh after refinement
    if (n_sub_mesh != 0)
    {
      // Rebuild the global mesh
      rebuild_global_mesh();
    }
  }


  //====================================================================
  /// Helper function to re-setup the Base_mesh enumeration
  /// (used during load balancing) after pruning.
  //====================================================================
  void Problem::setup_base_mesh_info_after_pruning()
  {
    // Storage for number of processors and current processor
    int n_proc = this->communicator_pt()->nproc();
    int my_rank = this->communicator_pt()->my_rank();

    // Loop over sub meshes
    unsigned n_sub_mesh = nsub_mesh();
    unsigned max_mesh = std::max(n_sub_mesh, unsigned(1));
    for (unsigned i_mesh = 0; i_mesh < max_mesh; i_mesh++)
    {
      // Choose the right mesh
      Mesh* my_mesh_pt = 0;
      if (n_sub_mesh == 0)
      {
        my_mesh_pt = mesh_pt();
      }
      else
      {
        my_mesh_pt = mesh_pt(i_mesh);
      }

      // Only work with structured meshes
      TriangleMeshBase* sub_mesh_pt =
        dynamic_cast<TriangleMeshBase*>(my_mesh_pt);
      if (!(sub_mesh_pt != 0))
      {
        // Storage for number of data to be sent to each processor
        Vector<int> send_n(n_proc, 0);

        // Storage for all values to be sent to all processors
        Vector<unsigned> send_data;

        // Start location within send_data for data to be sent to each processor
        Vector<int> send_displacement(n_proc, 0);

        // Loop over all processors
        for (int rank = 0; rank < n_proc; rank++)
        {
          // Set the offset for the current processor
          send_displacement[rank] = send_data.size();

          // Don't bother to do anything if the processor in the loop is the
          // current processor
          if (rank != my_rank)
          {
            // Get root haloed elements with that processor
            Vector<GeneralisedElement*> root_haloed_elements_pt =
              my_mesh_pt->root_haloed_element_pt(rank);
            unsigned nel = root_haloed_elements_pt.size();

            // Store element numbers for send
            for (unsigned e = 0; e < nel; e++)
            {
              GeneralisedElement* el_pt = root_haloed_elements_pt[e];
              send_data.push_back(Base_mesh_element_number_plus_one[el_pt]);
            }
          }

          // Find the number of data added to the vector
          send_n[rank] = send_data.size() - send_displacement[rank];
        }

        // Storage for the number of data to be received from each processor
        Vector<int> receive_n(n_proc, 0);

        // Now send numbers of data to be sent between all processors
        MPI_Alltoall(&send_n[0],
                     1,
                     MPI_INT,
                     &receive_n[0],
                     1,
                     MPI_INT,
                     this->communicator_pt()->mpi_comm());

        // We now prepare the data to be received
        // by working out the displacements from the received data
        Vector<int> receive_displacement(n_proc, 0);
        int receive_data_count = 0;
        for (int rank = 0; rank < n_proc; ++rank)
        {
          // Displacement is number of data received so far
          receive_displacement[rank] = receive_data_count;
          receive_data_count += receive_n[rank];
        }

        // Now resize the receive buffer for all data from all processors
        // Make sure that it has a size of at least one
        if (receive_data_count == 0)
        {
          ++receive_data_count;
        }
        Vector<unsigned> receive_data(receive_data_count);

        // Make sure that the send buffer has size at least one
        // so that we don't get a segmentation fault
        if (send_data.size() == 0)
        {
          send_data.resize(1);
        }

        // Now send the data between all the processors
        MPI_Alltoallv(&send_data[0],
                      &send_n[0],
                      &send_displacement[0],
                      MPI_UNSIGNED,
                      &receive_data[0],
                      &receive_n[0],
                      &receive_displacement[0],
                      MPI_UNSIGNED,
                      this->communicator_pt()->mpi_comm());

        // Now use the received data to update the halo element numbers in
        // base mesh
        for (int send_rank = 0; send_rank < n_proc; send_rank++)
        {
          // Don't bother to do anything for the processor corresponding to the
          // current processor or if no data were received from this processor
          if ((send_rank != my_rank) && (receive_n[send_rank] != 0))
          {
            // Counter for the data within the large array
            unsigned count = receive_displacement[send_rank];

            // Get root halo elements with that processor
            Vector<GeneralisedElement*> root_halo_elements_pt =
              my_mesh_pt->root_halo_element_pt(send_rank);
            unsigned nel = root_halo_elements_pt.size();

            // Read in element numbers
            for (unsigned e = 0; e < nel; e++)
            {
              GeneralisedElement* el_pt = root_halo_elements_pt[e];
              unsigned el_number_plus_one = receive_data[count++];
              Base_mesh_element_number_plus_one[el_pt] = el_number_plus_one;
              Base_mesh_element_pt[el_number_plus_one - 1] = el_pt;
            }
          }

        } // End of data is received

      } // if (!(sub_mesh_pt!=0))

    } // for (i_mesh<max_mesh)
  }

#endif

  /// Instantiation of public flag to allow suppression of warning
  /// messages re reading in unstructured meshes during restart.
  bool Problem::Suppress_warning_about_actions_before_read_unstructured_meshes =
    false;


} // namespace oomph
