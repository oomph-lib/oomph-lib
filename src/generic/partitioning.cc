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
#include <float.h>

#include "partitioning.h"
#include "mesh.h"
#include "refineable_mesh.h"
// Include to fill in additional_setup_shared_node_scheme() function
#include "refineable_mesh.template.cc"

#ifdef OOMPH_TRANSITION_TO_VERSION_3

// for the new METIS API, need to use symbols defined in the standard header
// which aren't available in the current frozen (old) version of METIS
// Version 3 will (presumably) have this header in the include path as standard
#include "metis.h"

#endif

namespace oomph
{
  //====================================================================
  /// Namespace for METIS graph partitioning routines
  //====================================================================
  namespace METIS
  {
    /// \short Default function that translates spatial
    /// error into weight for METIS partitioning (unit weight regardless
    /// of input).
    void default_error_to_weight_fct(const double& spatial_error,
                                     const double& max_error,
                                     const double& min_error,
                                     int& weight)
    {
      weight = 1;
    }

    /// \short Function pointer to to function that translates spatial
    /// error into weight for METIS partitioning.
    ErrorToWeightFctPt Error_to_weight_fct_pt = &default_error_to_weight_fct;

  } // namespace METIS


  //==================================================================
  /// Partition mesh uniformly by dividing elements
  /// equally over the partitions, in the order
  /// in which they are returned by problem.
  /// On return, element_domain[ielem] contains the number
  /// of the domain [0,1,...,ndomain-1] to which
  /// element ielem has been assigned.
  //==================================================================
  void METIS::uniform_partition_mesh(Problem* problem_pt,
                                     const unsigned& ndomain,
                                     Vector<unsigned>& element_domain)
  {
    // Number of elements
    unsigned nelem = problem_pt->mesh_pt()->nelement();

#ifdef PARANOID
    if (nelem != element_domain.size())
    {
      std::ostringstream error_stream;
      error_stream << "element_domain Vector has wrong length " << nelem << " "
                   << element_domain.size() << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // Uniform partitioning
    unsigned nel_per_domain = int(float(nelem) / float(ndomain));
    for (unsigned ielem = 0; ielem < nelem; ielem++)
    {
      unsigned idomain = unsigned(float(ielem) / float(nel_per_domain));
      element_domain[ielem] = idomain;
    }
  }


  //==================================================================
  /// Use METIS to assign each element to a domain.
  /// On return, element_domain[ielem] contains the number
  /// of the domain [0,1,...,ndomain-1] to which
  /// element ielem has been assigned.
  /// - objective=0: minimise edgecut.
  /// - objective=1: minimise total communications volume.
  /// .
  /// Partioning is based on dual graph of mesh.
  //==================================================================
  void METIS::partition_mesh(Problem* problem_pt,
                             const unsigned& ndomain,
                             const unsigned& objective,
                             Vector<unsigned>& element_domain)
  {
    // Global mesh
    Mesh* mesh_pt = problem_pt->mesh_pt();

    // Number of elements
    unsigned nelem = mesh_pt->nelement();

#ifdef PARANOID
    if (nelem != element_domain.size())
    {
      std::ostringstream error_stream;
      error_stream << "element_domain Vector has wrong length " << nelem << " "
                   << element_domain.size() << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Setup dual graph
    //------------------

    // Start timer
    clock_t cpu_start = clock();

    // Container to collect all elements associated with given global eqn number
    std::map<unsigned, std::set<unsigned>> elements_connected_with_global_eqn;

    // Container for all unique global eqn numbers
    std::set<unsigned> all_global_eqns;

    // Loop over all elements
    for (unsigned e = 0; e < nelem; e++)
    {
      GeneralisedElement* el_pt = mesh_pt->element_pt(e);

      // Add all global eqn numbers
      unsigned ndof = el_pt->ndof();
      for (unsigned j = 0; j < ndof; j++)
      {
        // Get global eqn number
        unsigned eqn_number = el_pt->eqn_number(j);
        elements_connected_with_global_eqn[eqn_number].insert(e);
        all_global_eqns.insert(eqn_number);
      }
    }

    // Now reverse the lookup scheme to find out all elements
    // that are connected because they share the same global eqn
    Vector<std::set<unsigned>> connected_elements(nelem);

    // Counter for total number of entries in connected_elements structure
    unsigned count = 0;

    // Loop over all global eqns
    for (std::set<unsigned>::iterator it = all_global_eqns.begin();
         it != all_global_eqns.end();
         it++)
    {
      // Get set of elements connected with this data item
      std::set<unsigned> elements = elements_connected_with_global_eqn[*it];

      // Double loop over connnected elements: Everybody's connected to
      // everybody
      for (std::set<unsigned>::iterator it1 = elements.begin();
           it1 != elements.end();
           it1++)
      {
        for (std::set<unsigned>::iterator it2 = elements.begin();
             it2 != elements.end();
             it2++)
        {
          if ((*it1) != (*it2))
          {
            connected_elements[(*it1)].insert(*it2);
          }
        }
      }
    }


    // Now convert into C-style packed array for interface with METIS
    int* xadj = new int[nelem + 1];
    Vector<int> adjacency_vector;

    // Reserve (too much) space
    adjacency_vector.reserve(count);

    // Initialise counters
    unsigned ientry = 0;

    // Loop over all elements
    for (unsigned e = 0; e < nelem; e++)
    {
      // First entry for current element
      xadj[e] = ientry;

      // Loop over elements that are connected to current element
      typedef std::set<unsigned>::iterator IT;
      for (IT it = connected_elements[e].begin();
           it != connected_elements[e].end();
           it++)
      {
        // Copy into adjacency array
        adjacency_vector.push_back(*it);

        // We've just made another entry
        ientry++;
      }

      // Entry after last entry for current element:
      xadj[e + 1] = ientry;
    }

    // End timer
    clock_t cpu_end = clock();

    // Doc
    double cpu0 = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;
    oomph_info
      << "CPU time for setup of METIS data structures            [nelem="
      << nelem << "]: " << cpu0 << " sec" << std::endl;


    // Call METIS graph partitioner
    //-----------------------------

    // Start timer
    cpu_start = clock();

    // Number of vertices in graph
    int nvertex = nelem;

    // No vertex weights
    int* vwgt = 0;

    // No edge weights
    int* adjwgt = 0;

    // Flag indicating that graph isn't weighted: 0; vertex weights only: 2
    // Note that wgtflag==2 requires nodal weights to be stored in vwgt.
    int wgtflag = 0;

    // Use C-style numbering (first array entry is zero)
    int numflag = 0;

    // Number of desired partitions
    int nparts = ndomain;

    // Use default options
    int* options = new int[10];
    options[0] = 0;

#ifdef OOMPH_TRANSITION_TO_VERSION_3
    switch (objective)
    {
      case 0:
        // Edge-cut minimization
        options[0] = METIS_OBJTYPE_CUT;
        break;

      case 1:
        // communication volume minimisation
        options[0] = METIS_OBJTYPE_VOL;
        break;

      default:
        std::ostringstream error_stream;
        error_stream << "Wrong objective for METIS. objective = " << objective
                     << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Number of cut edges in graph
    int* edgecut = new int[nelem];

    // Array containing the partition information
    int* part = new int[nelem];

    // Can we get an error estimate?

    unsigned n_mesh = problem_pt->nsub_mesh();

    if (n_mesh == 0)
    {
      RefineableMeshBase* mmesh_pt = dynamic_cast<RefineableMeshBase*>(mesh_pt);
      if (mmesh_pt != 0)
      {
        // Bias distribution?
        if (Error_to_weight_fct_pt != &default_error_to_weight_fct)
        {
          oomph_info
            << "Biasing element distribution via spatial error estimate\n";

          // Adjust flag and provide storage for weights
          wgtflag = 2;
          vwgt = new int[nelem];

          // Get error for all elements
          Vector<double> elemental_error(nelem);
          mmesh_pt->spatial_error_estimator_pt()->get_element_errors(
            mesh_pt, elemental_error);

          double max_error =
            *(std::max_element(elemental_error.begin(), elemental_error.end()));
          double min_error =
            *(std::min_element(elemental_error.begin(), elemental_error.end()));

          // Bias weights
          int weight = 1;
          for (unsigned e = 0; e < nelem; e++)
          {
            // Translate error into weight
            Error_to_weight_fct_pt(
              elemental_error[e], max_error, min_error, weight);
            vwgt[e] = weight;
          }
        }
      }
    }
    else // There are submeshes
    {
      // Are any of the submeshes refineable?
      bool refineable_submesh_exists = false;
      // Vector to store "start and end point" for loops in submeshes
      Vector<unsigned> loop_helper(n_mesh + 1);
      loop_helper[0] = 0;

      // Loop over submeshes
      for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
      {
        // Store the end of the loop
        loop_helper[i_mesh + 1] =
          problem_pt->mesh_pt(i_mesh)->nelement() + loop_helper[i_mesh];

        RefineableMeshBase* mmesh_pt =
          dynamic_cast<RefineableMeshBase*>(problem_pt->mesh_pt(i_mesh));
        if (mmesh_pt != 0)
        {
          refineable_submesh_exists = true;
        }
      }

      // If a refineable submesh exists
      if (refineable_submesh_exists)
      {
        // Bias distribution?
        if (Error_to_weight_fct_pt != &default_error_to_weight_fct)
        {
          oomph_info
            << "Biasing element distribution via spatial error estimate\n";

          // Adjust flag and provide storage for weights
          wgtflag = 2;
          vwgt = new int[nelem];

          // Loop over submeshes
          for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
          {
            RefineableMeshBase* mmesh_pt =
              dynamic_cast<RefineableMeshBase*>(problem_pt->mesh_pt(i_mesh));
            if (mmesh_pt != 0)
            {
              // Get error for all elements
              unsigned nsub_elem =
                loop_helper[i_mesh + 1] - loop_helper[i_mesh];
              Vector<double> elemental_error(nsub_elem);
              mmesh_pt->spatial_error_estimator_pt()->get_element_errors(
                problem_pt->mesh_pt(i_mesh), elemental_error);

              double max_error = *(std::max_element(elemental_error.begin(),
                                                    elemental_error.end()));
              double min_error = *(std::min_element(elemental_error.begin(),
                                                    elemental_error.end()));

              // Bias weights
              int weight = 1;
              unsigned start = loop_helper[i_mesh];
              unsigned end = loop_helper[i_mesh + 1];
              for (unsigned e = start; e < end; e++)
              {
                unsigned error_index = e - start;
                // Translate error into weight
                Error_to_weight_fct_pt(
                  elemental_error[error_index], max_error, min_error, weight);
                vwgt[e] = weight;
              }
            }
            else // This mesh is not refineable
            {
              // There's no error estimator, so use the default weight
              int weight = 1;
              unsigned start = loop_helper[i_mesh];
              unsigned end = loop_helper[i_mesh + 1];
              for (unsigned e = start; e < end; e++)
              {
                vwgt[e] = weight;
              }
            }
          }
        }
      }
    }

#ifdef OOMPH_TRANSITION_TO_VERSION_3

    // Call partitioner
    METIS_PartGraphKway(&nvertex,
                        xadj,
                        &adjacency_vector[0],
                        vwgt,
                        adjwgt,
                        &wgtflag,
                        &numflag,
                        &nparts,
                        options,
                        edgecut,
                        part);
#else
    // original code to delete in version 3

    // Call partitioner
    if (objective == 0)
    {
      // Partition with the objective of minimising the edge cut
      METIS_PartGraphKway(&nvertex,
                          xadj,
                          &adjacency_vector[0],
                          vwgt,
                          adjwgt,
                          &wgtflag,
                          &numflag,
                          &nparts,
                          options,
                          edgecut,
                          part);
    }
    else if (objective == 1)
    {
      // Partition with the objective of minimising the total communication
      // volume
      METIS_PartGraphVKway(&nvertex,
                           xadj,
                           &adjacency_vector[0],
                           vwgt,
                           adjwgt,
                           &wgtflag,
                           &numflag,
                           &nparts,
                           options,
                           edgecut,
                           part);
    }
    else
    {
      std::ostringstream error_stream;
      error_stream << "Wrong objective for METIS. objective = " << objective
                   << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

#ifdef PARANOID
    std::vector<bool> done(nparts, false);
#endif

    // Copy across
    for (unsigned e = 0; e < nelem; e++)
    {
      element_domain[e] = part[e];
#ifdef PARANOID
      done[part[e]] = true;
#endif
    }


#ifdef PARANOID
    // Check
    std::ostringstream error_stream;
    bool shout = false;
    for (int p = 0; p < nparts; p++)
    {
      if (!done[p])
      {
        shout = true;
        error_stream << "No elements on processor " << p
                     << "when trying to partition " << nelem << "elements over "
                     << nparts << " processors!\n";
      }
    }
    if (shout)
    {
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // End timer
    cpu_end = clock();

    // Doc
    double cpu1 = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;
    oomph_info
      << "CPU time for METIS mesh partitioning                   [nelem="
      << nelem << "]: " << cpu1 << " sec" << std::endl;


    // Cleanup
    delete[] xadj;
    delete[] part;
    delete[] edgecut;
    delete[] options;
  }


#ifdef OOMPH_HAS_MPI


  //==================================================================
  /// \short Use METIS to assign each element in an already-distributed mesh
  /// to a domain. On return, element_domain_on_this_proc[e] contains the number
  /// of the domain [0,1,...,ndomain-1] to which non-halo element e on THE
  /// CURRENT PROCESSOR ONLY has been assigned. The order of the non-halo
  /// elements is the same as in the Problem's mesh, with the halo
  /// elements being skipped.
  /// Objective:
  /// - objective=0: minimise edgecut.
  /// - objective=1: minimise total communications volume.
  /// .
  /// The partioning is based on the dof graph of the complete mesh by
  /// taking into
  /// account which global equation numbers are affected by each element and
  /// connecting elements which affect the same global equation number.
  /// Partitioning is done such that all elements associated with the
  /// same tree root move together. Non-refineable elements are
  /// treated as their own root elements. If the optional boolean
  /// flag is set to true (it defaults to false) each processor
  /// assigns a dumb-but-repeatable equidistribution of its non-halo
  /// elements over the domains and outputs the input that would have
  /// gone into METIS in the file metis_input_for_validation.dat
  //==================================================================
  void METIS::partition_distributed_mesh(
    Problem* problem_pt,
    const unsigned& objective,
    Vector<unsigned>& element_domain_on_this_proc,
    const bool& bypass_metis)
  {
    // Start timer
    clock_t cpu_start = clock();

    // Communicator
    OomphCommunicator* comm_pt = problem_pt->communicator_pt();

    // Number of processors / domains
    unsigned n_proc = comm_pt->nproc();
    unsigned my_rank = comm_pt->my_rank();

    // Global mesh
    Mesh* mesh_pt = problem_pt->mesh_pt();

    // Total number of elements (halo and nonhalo) on this proc
    unsigned n_elem = mesh_pt->nelement();

    // Get elemental assembly times
    Vector<double> elemental_assembly_time =
      problem_pt->elemental_assembly_time();

#ifdef PARANOID
    unsigned n = elemental_assembly_time.size();
    if ((n != 0) && (n != n_elem))
    {
      std::ostringstream error_stream;
      error_stream << "Number of elements doesn't match the \n"
                   << "number of elemental assembly times: " << n_elem << " "
                   << n << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Can we base load balancing on assembly times?
    bool can_load_balance_on_assembly_times = false;
    if (elemental_assembly_time.size() != 0)
    {
      can_load_balance_on_assembly_times = true;
    }

    // Storage for global eqn numbers on current processor
    std::set<unsigned> global_eqns_on_this_proc;

    // Storage for pointers to root elements that are connected with given
    // eqn number -- assembled on local processor
    std::map<unsigned, std::set<GeneralisedElement*>>
      root_elements_connected_with_global_eqn_on_this_proc;

    // Storage for long sequence of equation numbers as encountered
    // by the root elements on this processor
    Vector<unsigned> eqn_numbers_with_root_elements_on_this_proc;

    // Reserve number of elements x average/estimate (?) for number of dofs
    // per element
    eqn_numbers_with_root_elements_on_this_proc.reserve(n_elem * 9);

    // Storage for the number of eqn numbers associated with each
    // root element on this processors -- once this and the previous
    // container have been collected from all processors we're
    // able to reconstruct which root element (in the nominal "global" mesh)
    // is connected with which global equations
    Vector<unsigned> number_of_dofs_for_root_element;
    number_of_dofs_for_root_element.reserve(n_elem);

    // Ditto for number of "leaf" elements connected with each root
    Vector<unsigned> number_of_non_halo_elements_for_root_element;
    number_of_non_halo_elements_for_root_element.reserve(n_elem);

    // Ditto for total assembly time of "leaf" elements connected with each root
    Vector<double> total_assembly_time_for_root_element;
    total_assembly_time_for_root_element.reserve(n_elem);

    // Map storing the number of the root elements on this processor
    // (offset by one to bypass the zero default).
    std::map<GeneralisedElement*, unsigned> root_el_number_plus_one;

    // Loop over non-halo elements on this processor
    int number_of_root_elements = 0;
    unsigned number_of_non_halo_elements = 0;
    for (unsigned e = 0; e < n_elem; e++)
    {
      double el_assembly_time = 0.0;
      GeneralisedElement* el_pt = mesh_pt->element_pt(e);
      if (!el_pt->is_halo())
      {
        if (can_load_balance_on_assembly_times)
        {
          el_assembly_time = elemental_assembly_time[e];
        }

        // Get the associated root element which is either...
        GeneralisedElement* root_el_pt = 0;
        RefineableElement* ref_el_pt = dynamic_cast<RefineableElement*>(el_pt);
        if (ref_el_pt != 0)
        {
          //...the actual root element
          root_el_pt = ref_el_pt->root_element_pt();
        }
        // ...or the element itself
        else
        {
          root_el_pt = el_pt;
        }

        // Have we already encountered this root element?
        // (offset of one to bypass the default return of zero)
        bool already_encountered = false;
        unsigned root_el_number = root_el_number_plus_one[root_el_pt];
        if (root_el_number_plus_one[root_el_pt] == 0)
        {
          // This is a new one
          already_encountered = false;

          // Give it a number
          number_of_root_elements++;
          root_el_number_plus_one[root_el_pt] = number_of_root_elements;

          // Remove offset
          root_el_number = number_of_root_elements - 1;
        }
        else
        {
          // We've already visited this one before...
          already_encountered = true;

          // Remove offset
          root_el_number -= 1;
        }


        // Get global equation numbers of actual element
        unsigned n_dof = el_pt->ndof();
        for (unsigned i = 0; i < n_dof; i++)
        {
          unsigned eqn_no = el_pt->eqn_number(i);

          // Record which root elements are connected with this eqn number
          root_elements_connected_with_global_eqn_on_this_proc[eqn_no].insert(
            root_el_pt);

          // Record all global eqn numbers on this processor
          global_eqns_on_this_proc.insert(eqn_no);

          // Add eqn number of the current element to the long sequence
          // of eqn numbers
          eqn_numbers_with_root_elements_on_this_proc.push_back(eqn_no);
        }

        // Now record how many equations are associated with the current
        // non-halo element
        if (already_encountered)
        {
          number_of_dofs_for_root_element[root_el_number] += n_dof;
          number_of_non_halo_elements_for_root_element[root_el_number]++;
          total_assembly_time_for_root_element[root_el_number] +=
            el_assembly_time;
        }
        else
        {
          number_of_dofs_for_root_element.push_back(n_dof);
          number_of_non_halo_elements_for_root_element.push_back(1);
          total_assembly_time_for_root_element.push_back(el_assembly_time);
        }

        // Bump up number of non-halos
        number_of_non_halo_elements++;
      }
    }

    // Tell everybody how many root elements
    // are on each processor
    unsigned root_processor = 0;
    Vector<int> number_of_root_elements_on_each_proc(n_proc, 0);
    MPI_Allgather(&number_of_root_elements,
                  1,
                  MPI_INT,
                  &number_of_root_elements_on_each_proc[0],
                  1,
                  MPI_INT,
                  comm_pt->mpi_comm());


    // In the big sequence of concatenated root elements (enumerated
    // individually on the various processors) where do the root elements from a
    // given processor start? Also figure out how many root elements there are
    // in total by summing up their numbers
    Vector<int> start_index(n_proc, 0);
    unsigned total_number_of_root_elements = 0;
    for (unsigned i_proc = 0; i_proc < n_proc; i_proc++)
    {
      total_number_of_root_elements +=
        number_of_root_elements_on_each_proc[i_proc];
      if (i_proc != 0)
      {
        start_index[i_proc] = total_number_of_root_elements -
                              number_of_root_elements_on_each_proc[i_proc];
      }
      else
      {
        start_index[0] = 0;
      }
    }


    // How many global equations are held on this processor?
    int n_eqns_on_this_proc =
      eqn_numbers_with_root_elements_on_this_proc.size();

    // Gather this information for all processors:
    // n_eqns_on_each_proc[iproc] now contains the number of global
    // equations held on processor iproc.
    Vector<int> n_eqns_on_each_proc(n_proc, 0);
    MPI_Allgather(&n_eqns_on_this_proc,
                  1,
                  MPI_INT,
                  &n_eqns_on_each_proc[0],
                  1,
                  MPI_INT,
                  comm_pt->mpi_comm());


    // In the big sequence of equation numbers from the root elements
    // (enumerated individually on the various processors) where do the
    // equation numbers associated with the root elements from a given
    // processor start? Also figure out how long the sequence of equation
    // numbers is
    Vector<int> start_eqns_index(n_proc, 0);
    unsigned total_n_eqn = 0;
    for (unsigned i_proc = 0; i_proc < n_proc; i_proc++)
    {
      total_n_eqn += n_eqns_on_each_proc[i_proc];
      if (i_proc != 0)
      {
        start_eqns_index[i_proc] = total_n_eqn - n_eqns_on_each_proc[i_proc];
      }
      else
      {
        start_eqns_index[0] = 0;
      }
    }


    // Big vector that contains the number of dofs for each root element
    // (concatenated in processor-by-processor order)
    Vector<unsigned> number_of_dofs_for_global_root_element(
      total_number_of_root_elements);
    // Create at least one entry so we don't get a seg fault below
    if (number_of_dofs_for_root_element.size() == 0)
    {
      number_of_dofs_for_root_element.resize(1);
    }
    MPI_Gatherv(
      &number_of_dofs_for_root_element[0], // pointer to first entry in
                                           // vector to be gathered on root
      number_of_root_elements, // Number of entries to be sent
                               // from current processor
      MPI_UNSIGNED,
      &number_of_dofs_for_global_root_element[0], // Target -- this will
                                                  // store the concatenated
                                                  // vectors sent from
                                                  // everywhere
      &number_of_root_elements_on_each_proc[0], // Pointer to
                                                // vector containing
                                                // the length of the
                                                // vectors received
                                                // from elsewhere
      &start_index[0], // "offset" for storage of vector received
                       // from various processors in the global
                       // concatenated vector stored on root
      MPI_UNSIGNED,
      root_processor,
      comm_pt->mpi_comm());


    // ditto for number of non-halo elements associated with root element
    Vector<unsigned> number_of_non_halo_elements_for_global_root_element(
      total_number_of_root_elements);

    // Create at least one entry so we don't get a seg fault below
    if (number_of_non_halo_elements_for_root_element.size() == 0)
    {
      number_of_non_halo_elements_for_root_element.resize(1);
    }
    MPI_Gatherv(&number_of_non_halo_elements_for_root_element[0],
                // pointer to first entry in
                // vector to be gathered on root
                number_of_root_elements, // Number of entries to be sent
                                         // from current processor
                MPI_UNSIGNED,
                &number_of_non_halo_elements_for_global_root_element[0],
                // Target -- this will
                // store the concatenated
                // vectors sent from
                // everywhere
                &number_of_root_elements_on_each_proc[0], // Pointer to
                                                          // vector containing
                                                          // the length of the
                                                          // vectors received
                                                          // from elsewhere
                &start_index[0], // "offset" for storage of vector received
                                 // from various processors in the global
                                 // concatenated vector stored on root
                MPI_UNSIGNED,
                root_processor,
                comm_pt->mpi_comm());


    // ditto for assembly times elements associated with root element
    Vector<double> total_assembly_time_for_global_root_element(
      total_number_of_root_elements);

    // Create at least one entry so we don't get a seg fault below
    if (total_assembly_time_for_root_element.size() == 0)
    {
      total_assembly_time_for_root_element.resize(1);
    }
    MPI_Gatherv(&total_assembly_time_for_root_element[0],
                // pointer to first entry in
                // vector to be gathered on root
                number_of_root_elements, // Number of entries to be sent
                                         // from current processor
                MPI_DOUBLE,
                &total_assembly_time_for_global_root_element[0],
                // Target -- this will
                // store the concatenated
                // vectors sent from
                // everywhere
                &number_of_root_elements_on_each_proc[0], // Pointer to
                                                          // vector containing
                                                          // the length of the
                                                          // vectors received
                                                          // from elsewhere
                &start_index[0], // "offset" for storage of vector received
                                 // from various processors in the global
                                 // concatenated vector stored on root
                MPI_DOUBLE,
                root_processor,
                comm_pt->mpi_comm());


    // Big vector to store the long sequence of global equation numbers
    // associated with the long sequence of root elements
    Vector<unsigned> eqn_numbers_with_root_elements(total_n_eqn);

    // Create at least one entry so we don't get a seg fault below
    if (eqn_numbers_with_root_elements_on_this_proc.size() == 0)
    {
      eqn_numbers_with_root_elements_on_this_proc.resize(1);
    }
    MPI_Gatherv(&eqn_numbers_with_root_elements_on_this_proc[0],
                n_eqns_on_this_proc,
                MPI_UNSIGNED,
                &eqn_numbers_with_root_elements[0],
                &n_eqns_on_each_proc[0],
                &start_eqns_index[0],
                MPI_UNSIGNED,
                root_processor,
                comm_pt->mpi_comm());

    // Doc
    clock_t cpu_end = clock();

    double cpu0 = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;
    oomph_info
      << "CPU time for global setup of METIS data structures [nroot_elem="
      << total_number_of_root_elements << "]: " << cpu0 << " sec" << std::endl;


    // Now the root processor has gathered all the data needed to establish
    // the root element connectivity (as in the serial case) so use METIS
    // to determine "partitioning" for non-uniformly refined mesh
    //----------------------------------------------------------------------

    // Vector to store target domain for each of the root elements (concatenated
    // in processor-by-processor order)
    Vector<unsigned> root_element_domain(total_number_of_root_elements, 0);
    if (my_rank == root_processor) //--
    {
      // Start timer
      clock_t cpu_start = clock();

      // Repeat the steps used in the serial code: Storage for
      // the global equations (on root processor)
      std::set<unsigned> all_global_eqns_root_processor;

      // Set of root elements (as enumerated in the processor-by-processor
      // order) associated with given global equation number
      std::map<unsigned, std::set<unsigned>>
        root_elements_connected_with_global_eqn_on_root_processor;

      // Retrace the steps of the serial code: Who's connected with who
      unsigned count_all = 0;
      for (unsigned e = 0; e < total_number_of_root_elements; e++)
      {
        unsigned n_eqn_no = number_of_dofs_for_global_root_element[e];
        for (unsigned n = 0; n < n_eqn_no; n++)
        {
          unsigned eqn_no = eqn_numbers_with_root_elements[count_all];
          count_all++;
          root_elements_connected_with_global_eqn_on_root_processor[eqn_no]
            .insert(e);
          all_global_eqns_root_processor.insert(eqn_no);
        }
      }

      // Number of domains
      unsigned ndomain = n_proc;

      // Now reverse the lookup scheme to find out all root elements
      // that are connected because they share the same global eqn
      Vector<std::set<unsigned>> connected_root_elements(
        total_number_of_root_elements);

      // Counter for total number of entries in connected_root_elements
      // structure
      unsigned count = 0;

      // Loop over all global eqns
      for (std::set<unsigned>::iterator it =
             all_global_eqns_root_processor.begin();
           it != all_global_eqns_root_processor.end();
           it++)
      {
        // Get set of root elements connected with this data item
        std::set<unsigned> root_elements =
          root_elements_connected_with_global_eqn_on_root_processor[*it];

        // Double loop over connnected root elements: Everybody's connected to
        // everybody
        for (std::set<unsigned>::iterator it1 = root_elements.begin();
             it1 != root_elements.end();
             it1++)
        {
          for (std::set<unsigned>::iterator it2 = root_elements.begin();
               it2 != root_elements.end();
               it2++)
          {
            if ((*it1) != (*it2))
            {
              connected_root_elements[(*it1)].insert(*it2);
            }
          }
        }
      }

      // End timer
      clock_t cpu_end = clock();

      // Doc
      double cpu0b = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;
      oomph_info << "CPU time for setup of connected elements (load balance) "
                    "[nroot_elem="
                 << total_number_of_root_elements << "]: " << cpu0b << " sec"
                 << std::endl;

      // Now convert into C-style packed array for interface with METIS
      cpu_start = clock();
      int* xadj = new int[total_number_of_root_elements + 1];
      Vector<int> adjacency_vector;

      // Reserve (too much) space
      adjacency_vector.reserve(count);

      // Initialise counters
      unsigned ientry = 0;

      // Loop over all elements
      for (unsigned e = 0; e < total_number_of_root_elements; e++)
      {
        // First entry for current element
        xadj[e] = ientry;

        // Loop over elements that are connected to current element
        typedef std::set<unsigned>::iterator IT;
        for (IT it = connected_root_elements[e].begin();
             it != connected_root_elements[e].end();
             it++)
        {
          // Copy into adjacency array
          adjacency_vector.push_back(*it);

          // We've just made another entry
          ientry++;
        }

        // Entry after last entry for current element:
        xadj[e + 1] = ientry;
      }

      // End timer
      cpu_end = clock();

      // Doc
      double cpu0 = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;
      oomph_info << "CPU time for setup of METIS data structures (load "
                    "balance) [nroot_elem="
                 << total_number_of_root_elements << "]: " << cpu0 << " sec"
                 << std::endl;


      // Call METIS graph partitioner
      //-----------------------------

      // Start timer
      cpu_start = clock();

      // Number of vertices in graph
      int nvertex = total_number_of_root_elements;

      // No vertex weights
      int* vwgt = 0;

      // No edge weights
      int* adjwgt = 0;

      // Flag indicating that graph isn't weighted: 0; vertex weights only: 2
      // Note that wgtflag==2 requires nodal weights to be stored in vwgt.
      int wgtflag = 0;

      // Use C-style numbering (first array entry is zero)
      int numflag = 0;

      // Number of desired partitions
      int nparts = ndomain;

      // Use default options
      int* options = new int[10];
      options[0] = 0;

#ifdef OOMPH_TRANSITION_TO_VERSION_3
      switch (objective)
      {
        case 0:
          // Edge-cut minimization
          options[0] = METIS_OBJTYPE_CUT;
          break;

        case 1:
          // communication volume minimisation
          options[0] = METIS_OBJTYPE_VOL;
          break;

        default:
          std::ostringstream error_stream;
          error_stream << "Wrong objective for METIS. objective = " << objective
                       << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Number of cut edges in graph
      int* edgecut = new int[total_number_of_root_elements];

      // Array containing the partition information
      int* part = new int[total_number_of_root_elements];

      // Now bias distribution by giving each root element
      // a weight equal to the number of elements associated with it

      // Adjust flag and provide storage for weights
      wgtflag = 2;
      vwgt = new int[total_number_of_root_elements];


      // Load balance based on assembly times of all leaf
      // elements associated with root
      if (can_load_balance_on_assembly_times)
      {
        oomph_info << "Basing distribution on assembly times of elements\n";

        // Normalise
        double min_time = *(
          std::min_element(total_assembly_time_for_global_root_element.begin(),
                           total_assembly_time_for_global_root_element.end()));
#ifdef PARANOID
        if (min_time == 0.0)
        {
          std::ostringstream error_stream;
          error_stream << "Minimum assemble time for element is zero!\n";
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Bypass METIS (usually for validation) and use made-up but
        // repeatable timings
        if (bypass_metis)
        {
          for (unsigned e = 0; e < total_number_of_root_elements; e++)
          {
            vwgt[e] = e;
          }
        }
        else
        {
          for (unsigned e = 0; e < total_number_of_root_elements; e++)
          {
            // Use assembly times (relative to minimum) as weight
            vwgt[e] =
              int(total_assembly_time_for_global_root_element[e] / min_time);
          }
        }
      }
      // Load balanced based on number of leaf elements associated with
      // root
      else
      {
        oomph_info << "Basing distribution on number of elements\n";
        for (unsigned e = 0; e < total_number_of_root_elements; e++)
        {
          vwgt[e] = number_of_non_halo_elements_for_global_root_element[e];
        }
      }

      // Bypass METIS (usually for validation)
      if (bypass_metis)
      {
        // Simple repeatable partition: Equidistribute root element
        for (unsigned e = 0; e < total_number_of_root_elements; e++)
        {
          // Simple repeatable partition: Equidistribute elements on each
          // processor
          part[e] = (n_proc - 1) -
                    unsigned(double(e) / double(total_number_of_root_elements) *
                             double(n_proc));
        }

        oomph_info
          << "Bypassing METIS for validation purposes.\n"
          << "Appending input for metis in metis_input_for_validation.dat\n";
        std::ofstream outfile;
        outfile.open("metis_input_for_validation.dat", std::ios_base::app);

        // Dump out relevant input to metis
        for (unsigned e = 0; e < total_number_of_root_elements + 1; e++)
        {
          outfile << xadj[e] << std::endl;
        }
        unsigned n = adjacency_vector.size();
        for (unsigned i = 0; i < n; i++)
        {
          outfile << adjacency_vector[i] << std::endl;
        }
        for (unsigned e = 0; e < total_number_of_root_elements; e++)
        {
          outfile << vwgt[e] << std::endl;
        }
        outfile.close();
      }
      // Actually use METIS (good but not always repeatable!)
      else
      {
#ifdef OOMPH_TRANSITION_TO_VERSION_3

        METIS_PartGraphKway(&nvertex,
                            xadj,
                            &adjacency_vector[0],
                            vwgt,
                            adjwgt,
                            &wgtflag,
                            &numflag,
                            &nparts,
                            options,
                            edgecut,
                            part);
#else
        // for old version of METIS; these two functions have been merged
        // in the new METIS API

        if (objective == 0)
        {
          // Partition with the objective of minimising the edge cut
          METIS_PartGraphKway(&nvertex,
                              xadj,
                              &adjacency_vector[0],
                              vwgt,
                              adjwgt,
                              &wgtflag,
                              &numflag,
                              &nparts,
                              options,
                              edgecut,
                              part);
        }
        else if (objective == 1)
        {
          // Partition with the objective of minimising the total communication
          // volume
          METIS_PartGraphVKway(&nvertex,
                               xadj,
                               &adjacency_vector[0],
                               vwgt,
                               adjwgt,
                               &wgtflag,
                               &numflag,
                               &nparts,
                               options,
                               edgecut,
                               part);
        }
#endif
      }

      // Copy across
      Vector<unsigned> total_weight_on_proc(n_proc, 0);
      for (unsigned e = 0; e < total_number_of_root_elements; e++)
      {
        root_element_domain[e] = part[e];
        total_weight_on_proc[part[e]] += vwgt[e];
      }

      // Document success of partitioning
      for (unsigned j = 0; j < n_proc; j++)
      {
        oomph_info << "Total weight on proc " << j << " is "
                   << total_weight_on_proc[j] << std::endl;
      }

      // Doc
      double cpu1 = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;
      oomph_info << "CPU time for METIS mesh partitioning [nroot_elem="
                 << total_number_of_root_elements << "]: " << cpu1 << " sec"
                 << std::endl;

      // Cleanup
      delete[] xadj;
      delete[] part;
      delete[] vwgt;
      delete[] edgecut;
      delete[] options;
    }

    // Now scatter things back to processors: root_element_domain[] contains
    // the target domain for all elements (concatenated in processor-by
    // processor order on the root processor). Distribute this back
    // to the processors so that root_element_domain_on_this_proc[e] contains
    // the target domain for root element e (in whatever order the processor
    // decided to line up its root elements).
    cpu_start = clock();
    Vector<unsigned> root_element_domain_on_this_proc(number_of_root_elements);

    // Create at least one entry so we don't get a seg fault below
    if (root_element_domain_on_this_proc.size() == 0)
    {
      root_element_domain_on_this_proc.resize(1);
    }
    MPI_Scatterv(&root_element_domain[0],
                 &number_of_root_elements_on_each_proc[0],
                 &start_index[0],
                 MPI_UNSIGNED,
                 &root_element_domain_on_this_proc[0],
                 number_of_root_elements,
                 MPI_UNSIGNED,
                 root_processor,
                 comm_pt->mpi_comm());


    // Now translate back into target domain for the actual (non-root)
    // elements
    element_domain_on_this_proc.resize(number_of_non_halo_elements);
    unsigned count_non_halo = 0;
    for (unsigned e = 0; e < n_elem; e++)
    {
      GeneralisedElement* el_pt = mesh_pt->element_pt(e);
      if (!el_pt->is_halo())
      {
        // Get the associated root element which is either...
        GeneralisedElement* root_el_pt = 0;
        RefineableElement* ref_el_pt = dynamic_cast<RefineableElement*>(el_pt);
        if (ref_el_pt != 0)
        {
          //...the actual root element
          root_el_pt = ref_el_pt->root_element_pt();
        }
        // ...or the element itself
        else
        {
          root_el_pt = el_pt;
        }

        // Recover the root element number (offset by one)
        unsigned root_el_number = root_el_number_plus_one[root_el_pt] - 1;

        // Copy target domain across from root element
        element_domain_on_this_proc[count_non_halo] =
          root_element_domain_on_this_proc[root_el_number];

        // Bump up counter for non-root elements
        count_non_halo++;
      }
    }


#ifdef PARANOID
    if (count_non_halo != number_of_non_halo_elements)
    {
      std::ostringstream error_stream;
      error_stream << "Non-halo counts don't match: " << count_non_halo << " "
                   << number_of_non_halo_elements << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // End timer
    cpu_end = clock();

    // Doc
    double cpu2 = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;
    oomph_info << "CPU time for communication of partition to all processors "
                  "[nroot_elem="
               << total_number_of_root_elements << "]: " << cpu2 << " sec"
               << std::endl;
  }


#endif

} // namespace oomph
