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
// Domain partitioning

#ifndef OOMPH_PARTITIONING_HEADER
#define OOMPH_PARTITIONING_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif
#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

// ooomph-lib includes
#include "Vector.h"
#include "problem.h"


namespace oomph
{
  //==================================================================
  // Interfaces to METIS functions
  //==================================================================
  extern "C"
  {
    /// Metis graph partitioning function
    void METIS_PartGraphKway(
      int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);

#ifndef OOMPH_TRANSITION_TO_VERSION_3
    // function from old API which no longer exists in METIS 5.1,
    // remove this moving to oomph-lib version 3

    /// Metis graph partitioning function -- decomposes
    /// nodal graph based on minimum communication volume
    void METIS_PartGraphVKway(
      int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);

#endif
  }


  //====================================================================
  /// Namespace for METIS graph partitioning routines
  //====================================================================
  namespace METIS
  {
    /// Default function that translates spatial
    /// error into weight for METIS partitioning (unit weight regardless
    /// of input)
    extern void default_error_to_weight_fct(const double& spatial_error,
                                            const double& max_error,
                                            const double& min_error,
                                            int& weight);

    /// Typedef for function pointer to to function that translates
    /// spatial error into weight for METIS partitioning.
    typedef void (*ErrorToWeightFctPt)(const double& spatial_error,
                                       const double& max_error,
                                       const double& min_error,
                                       int& weight);

    /// Function pointer to to function that translates spatial
    /// error into weight for METIS partitioning.
    extern ErrorToWeightFctPt Error_to_weight_fct_pt;

    /// Partition mesh uniformly by dividing elements
    /// equally over the partitions, in the order
    /// in which they are returned by problem.
    /// On return, element_domain[ielem] contains the number
    /// of the domain [0,1,...,ndomain-1] to which
    /// element ielem has been assigned.
    extern void uniform_partition_mesh(Problem* problem_pt,
                                       const unsigned& ndomain,
                                       Vector<unsigned>& element_domain);


    /// Use METIS to assign each element to a domain.
    /// On return, element_domain[ielem] contains the number
    /// of the domain [0,1,...,ndomain-1] to which
    /// element ielem has been assigned.
    /// - objective=0: minimise edgecut.
    /// - objective=1: minimise total communications volume.
    /// .
    /// Partioning is based on nodal graph of mesh.
    extern void partition_mesh(Problem* problem_pt,
                               const unsigned& ndomain,
                               const unsigned& objective,
                               Vector<unsigned>& element_domain);


    /// Use METIS to assign each element to a domain.
    /// On return, element_domain[ielem] contains the number
    /// of the domain [0,1,...,ndomain-1] to which
    /// element ielem has been assigned.
    /// - objective=0: minimise edgecut.
    /// - objective=1: minimise total communications volume.
    /// .
    /// Partioning is based on nodal graph of mesh.
    extern void partition_mesh(OomphCommunicator* comm_pt,
                               Mesh* mesh_pt,
                               const unsigned& ndomain,
                               const unsigned& objective,
                               Vector<unsigned>& element_domain);

    //  /// Use METIS to assign each element to a domain.
    //  /// On return, element_domain[ielem] contains the number
    //  /// of the domain [0,1,...,ndomain-1] to which
    //  /// element ielem has been assigned.
    //  /// - objective=0: minimise edgecut.
    //  /// - objective=1: minimise total communications volume.
    //  /// .
    //  /// Partioning is based on "Data" graph of mesh.
    //  extern void partition_mesh_data(Problem* problem_pt,
    //                                  const unsigned& ndomain,
    //                                  const unsigned& objective,
    //                                 Vector<unsigned>& element_domain);

#ifdef OOMPH_HAS_MPI


    /// Use METIS to assign each element in an already-distributed mesh
    /// to a domain. On return, element_domain_on_this_proc[e] contains the
    /// number of the domain [0,1,...,ndomain-1] to which non-halo element e on
    /// THE CURRENT PROCESSOR ONLY has been assigned. The order of the non-halo
    /// elements is the same as in the Problem's mesh, with the halo
    /// elements being skipped.
    ///
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
    extern void partition_distributed_mesh(
      Problem* problem_pt,
      const unsigned& objective,
      Vector<unsigned>& element_domain_on_this_proc,
      const bool& bypass_metis = false);

#endif

  } // namespace METIS


} // namespace oomph

#endif
