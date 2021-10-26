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
// Non-templated multi-domain functions which act on more than one mesh
// and set up the storage and interaction between the two

// oomph-lib header
#include "multi_domain.h"
#include "multi_domain.template.cc"
#include "mesh.h"
#include "algebraic_elements.h"
#include "macro_element_node_update_element.h"
#include "Qelements.h"

namespace oomph
{
  //======================================================================
  // Namespace for "global" multi-domain functions
  //======================================================================
  namespace Multi_domain_functions
  {
    /// Output file to document the boundary coordinate
    /// along the mesh boundary of the bulk mesh during call to
    /// setup_bulk_elements_adjacent_to_face_mesh(...)
    std::ofstream Doc_boundary_coordinate_file;

    // Workspace for locate zeta methods
    //----------------------------------

    /// Boolean to indicate that failure in setup multi domain
    /// functions is acceptable; defaults to false. If set to true
    /// external element pointers are set to null for those elements
    /// for which external elements couldn't be located.
    bool Accept_failed_locate_zeta_in_setup_multi_domain_interaction = false;

    /// Dimension of zeta tuples (set by get_dim_helper) -- needed
    /// because we store the scalar coordinates in flat-packed form.
    unsigned Dim;

    /// Lookup scheme for whether a local element's integration point
    /// has had an external element assigned to it -- essentially boolean.
    /// External_element_located[e][ipt] = {0,1} if external element
    /// for ipt-th integration in local element e {has not, has} been found.
    /// Used locally to ensure that we're not searching for the same
    /// elements over and over again when we go around the spirals.
    Vector<Vector<unsigned>> External_element_located;

    /// Vector of flat-packed zeta coordinates for which the external
    /// element could not be found during current local search. These
    /// will be sent to the next processor in the ring-like parallel search.
    /// The zeta coordinates come in groups of Dim (scalar) coordinates.
    Vector<double> Flat_packed_zetas_not_found_locally;

    /// Vector of flat-packed zeta coordinates for which the external
    /// element could not be found on another processor and for which
    /// we're currently searching here. Whatever can't be found here,
    /// gets written into Flat_packed_zetas_not_found_locally and then
    /// passed on to the next processor during the ring-like parallel search.
    /// The zeta coordinates come in  groups of Dim (scalar) coordinates.
    Vector<double> Received_flat_packed_zetas_to_be_found;

    /// Proc_id_plus_one_of_external_element[i] contains the
    /// processor id (plus one) of the processor
    /// on which the i-th zeta coordinate tuple received from elsewhere
    /// (in the order in which these are stored in
    /// Received_flat_packed_zetas_to_be_found) was located; it's zero if
    /// it wasn't found during the current stage of the ring-like parallel
    /// search.
    Vector<int> Proc_id_plus_one_of_external_element;

    /// Vector to indicate (to another processor) whether a
    /// located element (that will have to represented as an external
    /// halo element on that processor) should be newly created on that
    /// processor (2), already exists on that processor (1), or
    /// is not on the current processor either (0).
    Vector<unsigned> Located_element_status;

    /// Vector of flat-packed local coordinates for zeta tuples
    /// that have been located
    Vector<double> Flat_packed_located_coordinates;

    /// Vector of flat-packed doubles to be communicated with
    /// other processors
    Vector<double> Flat_packed_doubles;

    /// Counter used when processing vector of flat-packed
    /// doubles -- this is really "private" data, declared here
    /// to avoid having to pass it (and the associated array)
    /// between the various helper functions
    unsigned Counter_for_flat_packed_doubles;

    /// Vector of flat-packed unsigneds to be communicated with
    /// other processors -- this is really "private" data, declared here
    /// to avoid having to pass the array between the various helper
    /// functions
    Vector<unsigned> Flat_packed_unsigneds;

#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION

    // Temporary vector of strings to enable full annotation of multi domain
    // comms (but keep alive because it would be such a bloody pain to
    // rewrite it if things ever go wrong again...)
    Vector<std::string> Flat_packed_unsigneds_string;


#endif

    /// Counter used when processing vector of flat-packed
    /// unsigneds -- this is really "private" data, declared here
    /// to avoid having to pass it (and the associated array)
    /// between the various helper functions
    unsigned Counter_for_flat_packed_unsigneds;

    // Other parameters
    //-----------------

    /// Boolean to indicate when to use the bulk element as the
    /// external element.  Defaults to false, you must have set up FaceElements
    /// properly first in order for it to work
    bool Use_bulk_element_as_external = false;

    /// Boolean to indicate if we're allowed to use halo elements
    /// as external elements. Can drastically reduce the number of
    /// external halo elements -- currently not aware of any problems
    /// therefore set to true by default but retention
    /// of this flag allows easy return to previous implementation.
    bool Allow_use_of_halo_elements_as_external_elements = true;

    /// Indicate whether we are allowed to use halo elements as
    /// external elements for projection, possibly only required in
    /// parallel unstructured mesh generation during the projection
    /// stage. Default set to true
    bool Allow_use_of_halo_elements_as_external_elements_for_projection = true;

    /// Boolean to indicate whether to doc timings or not.
    bool Doc_timings = false;

    /// Boolean to indicate whether to output basic info during
    ///        setup_multi_domain_interaction() routines
    bool Doc_stats = false;

    /// Boolean to indicate whether to output further info during
    ///        setup_multi_domain_interaction() routines
    bool Doc_full_stats = false;

#ifdef OOMPH_HAS_MPI

    // Functions for location method in multi-domain problems


    //========================================================================
    /// Send the zeta coordinates from the current process to
    /// the next process; receive from the previous process
    //========================================================================
    void send_and_receive_missing_zetas(Problem* problem_pt)
    {
      // MPI info
      MPI_Status status;
      MPI_Request request;

      // Storage for number of processors, current process and communicator
      int n_proc = problem_pt->communicator_pt()->nproc();
      int my_rank = problem_pt->communicator_pt()->my_rank();
      OomphCommunicator* comm_pt = problem_pt->communicator_pt();

      // Work out processors to send and receive from
      int send_to_proc = my_rank + 1;
      int recv_from_proc = my_rank - 1;
      if (send_to_proc == n_proc)
      {
        send_to_proc = 0;
      }
      if (recv_from_proc < 0)
      {
        recv_from_proc = n_proc - 1;
      }

      // Send the number  of flat-packed zetas that we couldn't find
      // locally to the next processor
      int n_missing_local_zetas = Flat_packed_zetas_not_found_locally.size();
      MPI_Isend(&n_missing_local_zetas,
                1,
                MPI_INT,
                send_to_proc,
                4,
                comm_pt->mpi_comm(),
                &request);

      // Receive the number of flat-packed zetas that couldn't be found
      // on the "previous" processor
      int count_zetas = 0;
      MPI_Recv(&count_zetas,
               1,
               MPI_INT,
               recv_from_proc,
               4,
               comm_pt->mpi_comm(),
               &status);

      MPI_Wait(&request, MPI_STATUS_IGNORE);

      // Send the vector of flat-packed zetas that we couldn't find
      // locally to the next processor
      if (n_missing_local_zetas != 0)
      {
        MPI_Isend(&Flat_packed_zetas_not_found_locally[0],
                  n_missing_local_zetas,
                  MPI_DOUBLE,
                  send_to_proc,
                  5,
                  comm_pt->mpi_comm(),
                  &request);
      }

      // Receive the vector of flat-packed zetas that couldn't be found
      // on the "previous" processor
      if (count_zetas != 0)
      {
        Received_flat_packed_zetas_to_be_found.resize(count_zetas);
        MPI_Recv(&Received_flat_packed_zetas_to_be_found[0],
                 count_zetas,
                 MPI_DOUBLE,
                 recv_from_proc,
                 5,
                 comm_pt->mpi_comm(),
                 &status);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      }
      else
      {
        Received_flat_packed_zetas_to_be_found.resize(0);
      }

      // Now we should have the Zeta arrays set up correctly
      // for the next round of locations
    }

    //========start of send_and_receive_located_info==========================
    /// Send location information from current process; Received location
    /// information from (current process + iproc) modulo (nproc)
    //========================================================================
    void send_and_receive_located_info(int& iproc,
                                       Mesh* const& external_mesh_pt,
                                       Problem* problem_pt)
    {
      // Set MPI info
      MPI_Status status;
      MPI_Request request;

      // Storage for number of processors, current process and communicator
      OomphCommunicator* comm_pt = problem_pt->communicator_pt();
      int n_proc = comm_pt->nproc();
      int my_rank = comm_pt->my_rank();

      // Prepare vectors to receive information
      Vector<double> received_double_values;
      Vector<unsigned> received_unsigned_values;
      Vector<double> received_located_coord;
      Vector<int> received_proc_id_plus_one_of_external_element;
      Vector<unsigned> received_located_element_status;

      // Communicate the located information back to the original process
      int orig_send_proc = my_rank - iproc;
      if (my_rank < iproc)
      {
        orig_send_proc = n_proc + orig_send_proc;
      }
      int orig_recv_proc = my_rank + iproc;
      if ((my_rank + iproc) >= n_proc)
      {
        orig_recv_proc = orig_recv_proc - n_proc;
      }

      // Send the double values associated with external halos
      //------------------------------------------------------
      unsigned send_count_double_values = Flat_packed_doubles.size();
      MPI_Isend(&send_count_double_values,
                1,
                MPI_UNSIGNED,
                orig_send_proc,
                1,
                comm_pt->mpi_comm(),
                &request);
      int receive_count_double_values = 0;
      MPI_Recv(&receive_count_double_values,
               1,
               MPI_INT,
               orig_recv_proc,
               1,
               comm_pt->mpi_comm(),
               &status);
      MPI_Wait(&request, MPI_STATUS_IGNORE);

      if (send_count_double_values != 0)
      {
        MPI_Isend(&Flat_packed_doubles[0],
                  send_count_double_values,
                  MPI_DOUBLE,
                  orig_send_proc,
                  2,
                  comm_pt->mpi_comm(),
                  &request);
      }
      if (receive_count_double_values != 0)
      {
        received_double_values.resize(receive_count_double_values);
        MPI_Recv(&received_double_values[0],
                 receive_count_double_values,
                 MPI_DOUBLE,
                 orig_recv_proc,
                 2,
                 comm_pt->mpi_comm(),
                 &status);
      }
      if (send_count_double_values != 0)
      {
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      }

      // Now send unsigned values associated with external halos
      //---------------------------------------------------------
      unsigned send_count_unsigned_values = Flat_packed_unsigneds.size();
      MPI_Isend(&send_count_unsigned_values,
                1,
                MPI_UNSIGNED,
                orig_send_proc,
                14,
                comm_pt->mpi_comm(),
                &request);

      int receive_count_unsigned_values = 0;
      MPI_Recv(&receive_count_unsigned_values,
               1,
               MPI_INT,
               orig_recv_proc,
               14,
               comm_pt->mpi_comm(),
               &status);

      MPI_Wait(&request, MPI_STATUS_IGNORE);

      if (send_count_unsigned_values != 0)
      {
        MPI_Isend(&Flat_packed_unsigneds[0],
                  send_count_unsigned_values,
                  MPI_UNSIGNED,
                  orig_send_proc,
                  15,
                  comm_pt->mpi_comm(),
                  &request);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        for (unsigned i = 0; i < send_count_unsigned_values; i++)
        {
          oomph_info << "Sent:" << i << " to orig_proc:" << orig_send_proc
                     << " " << Flat_packed_unsigneds_string[i] << ": "
                     << Flat_packed_unsigneds[i] << std::endl;
        }
#endif
      }
      if (receive_count_unsigned_values != 0)
      {
        received_unsigned_values.resize(receive_count_unsigned_values);
        MPI_Recv(&received_unsigned_values[0],
                 receive_count_unsigned_values,
                 MPI_UNSIGNED,
                 orig_recv_proc,
                 15,
                 comm_pt->mpi_comm(),
                 &status);
      }

      if (send_count_unsigned_values != 0)
      {
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      }

      // Send and receive the Located_element_status
      //--------------------------------------------
      int send_count = Received_flat_packed_zetas_to_be_found.size() / Dim;
      MPI_Isend(&send_count,
                1,
                MPI_INT,
                orig_send_proc,
                20,
                comm_pt->mpi_comm(),
                &request);
      int receive_count = 0;
      MPI_Recv(&receive_count,
               1,
               MPI_INT,
               orig_recv_proc,
               20,
               comm_pt->mpi_comm(),
               &status);
      MPI_Wait(&request, MPI_STATUS_IGNORE);

      if (send_count != 0)
      {
        MPI_Isend(&Located_element_status[0],
                  send_count,
                  MPI_UNSIGNED,
                  orig_send_proc,
                  3,
                  comm_pt->mpi_comm(),
                  &request);
      }

      if (receive_count != 0)
      {
        received_located_element_status.resize(receive_count);
        MPI_Recv(&received_located_element_status[0],
                 receive_count,
                 MPI_UNSIGNED,
                 orig_recv_proc,
                 3,
                 comm_pt->mpi_comm(),
                 &status);
      }
      if (send_count != 0)
      {
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      }

      // Send and receive Proc_id_plus_one_of_external_element array
      //------------------------------------------------------------
      if (send_count != 0)
      {
        MPI_Isend(&Proc_id_plus_one_of_external_element[0],
                  send_count,
                  MPI_INT,
                  orig_send_proc,
                  13,
                  comm_pt->mpi_comm(),
                  &request);
      }
      if (receive_count != 0)
      {
        received_proc_id_plus_one_of_external_element.resize(receive_count);
        MPI_Recv(&received_proc_id_plus_one_of_external_element[0],
                 receive_count,
                 MPI_INT,
                 orig_recv_proc,
                 13,
                 comm_pt->mpi_comm(),
                 &status);
      }
      if (send_count != 0)
      {
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      }


      // And finally the Flat_packed_located_coordinates array
      //------------------------------------------------------
      unsigned send_count_located_coord =
        Flat_packed_located_coordinates.size();
      MPI_Isend(&send_count_located_coord,
                1,
                MPI_UNSIGNED,
                orig_send_proc,
                4,
                comm_pt->mpi_comm(),
                &request);
      unsigned receive_count_located_coord = 0;
      MPI_Recv(&receive_count_located_coord,
               1,
               MPI_UNSIGNED,
               orig_recv_proc,
               4,
               comm_pt->mpi_comm(),
               &status);
      MPI_Wait(&request, MPI_STATUS_IGNORE);

      if (send_count_located_coord != 0)
      {
        MPI_Isend(&Flat_packed_located_coordinates[0],
                  send_count_located_coord,
                  MPI_DOUBLE,
                  orig_send_proc,
                  5,
                  comm_pt->mpi_comm(),
                  &request);
      }
      if (receive_count_located_coord != 0)
      {
        received_located_coord.resize(receive_count_located_coord);
        MPI_Recv(&received_located_coord[0],
                 receive_count_located_coord,
                 MPI_DOUBLE,
                 orig_recv_proc,
                 5,
                 comm_pt->mpi_comm(),
                 &status);
      }
      if (send_count_located_coord != 0)
      {
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      }

      // Copy across into original containers -- these can now
      //------------------------------------------------------
      // be processed by create_external_halo_elements() to generate
      //------------------------------------------------------------
      // external halo elements
      //------------------------
      Flat_packed_doubles.resize(receive_count_double_values);
      for (int ii = 0; ii < receive_count_double_values; ii++)
      {
        Flat_packed_doubles[ii] = received_double_values[ii];
      }
      Flat_packed_unsigneds.resize(receive_count_unsigned_values);
      for (int ii = 0; ii < receive_count_unsigned_values; ii++)
      {
        Flat_packed_unsigneds[ii] = received_unsigned_values[ii];
      }
      Proc_id_plus_one_of_external_element.resize(receive_count);
      Located_element_status.resize(receive_count);
      for (int ii = 0; ii < receive_count; ii++)
      {
        Proc_id_plus_one_of_external_element[ii] =
          received_proc_id_plus_one_of_external_element[ii];
        Located_element_status[ii] = received_located_element_status[ii];
      }
      Flat_packed_located_coordinates.resize(receive_count_located_coord);
      for (int ii = 0; ii < int(receive_count_located_coord); ii++)
      {
        Flat_packed_located_coordinates[ii] = received_located_coord[ii];
      }
    }


    //========start of add_external_haloed_node_to_storage====================
    /// Helper function to add external haloed nodes, including any masters
    //========================================================================
    void add_external_haloed_node_to_storage(int& iproc,
                                             Node* nod_pt,
                                             Problem* problem_pt,
                                             Mesh* const& external_mesh_pt,
                                             int& n_cont_inter_values)
    {
      // Add the node if required
      add_external_haloed_node_helper(
        iproc, nod_pt, problem_pt, external_mesh_pt, n_cont_inter_values);

      // Recursively add any master nodes (and their master nodes etc)
      recursively_add_masters_of_external_haloed_node(
        iproc, nod_pt, problem_pt, external_mesh_pt, n_cont_inter_values);
    }


    //========================================================================
    /// Recursively add any master nodes (and their master nodes etc) of
    /// external nodes
    //========================================================================
    void recursively_add_masters_of_external_haloed_node(
      int& iproc,
      Node* nod_pt,
      Problem* problem_pt,
      Mesh* const& external_mesh_pt,
      int& n_cont_inter_values)
    {
      // Loop over continuously interpolated values and add masters
      for (int i_cont = -1; i_cont < n_cont_inter_values; i_cont++)
      {
        if (nod_pt->is_hanging(i_cont))
        {
          // Indicate that this node is a hanging node so the other
          // process knows to create HangInfo and masters, etc.
          Flat_packed_unsigneds.push_back(1);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          Flat_packed_unsigneds_string.push_back("Is hanging");
#endif
          // If this is a hanging node then add all its masters as
          // external halo nodes if they have not yet been added
          HangInfo* hang_pt = nod_pt->hanging_pt(i_cont);
          // Loop over masters
          unsigned n_master = hang_pt->nmaster();

          // Indicate number of master nodes to add on other process
          Flat_packed_unsigneds.push_back(n_master);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          Flat_packed_unsigneds_string.push_back("nmaster");
#endif
          for (unsigned m = 0; m < n_master; m++)
          {
            Node* master_nod_pt = hang_pt->master_node_pt(m);

            // Call the helper function for master nodes
            add_external_haloed_master_node_helper(iproc,
                                                   master_nod_pt,
                                                   problem_pt,
                                                   external_mesh_pt,
                                                   n_cont_inter_values);

            // Indicate the weight of this master
            Flat_packed_doubles.push_back(hang_pt->master_weight(m));

            // Recursively add any master nodes (and their master nodes etc)
            recursively_add_masters_of_external_haloed_node(
              iproc,
              master_nod_pt,
              problem_pt,
              external_mesh_pt,
              n_cont_inter_values);
          }
        }
        else
        {
          // Indicate that it's not a hanging node in this variable
          Flat_packed_unsigneds.push_back(0);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          Flat_packed_unsigneds_string.push_back("Not hanging");
#endif
        }
      } // end loop over continously interpolated values
    }

    //==========start of add_external_haloed_node_helper======================
    /// Helper to add external haloed node that is not a master
    //========================================================================
    void add_external_haloed_node_helper(int& iproc,
                                         Node* nod_pt,
                                         Problem* problem_pt,
                                         Mesh* const& external_mesh_pt,
                                         int& n_cont_inter_values)
    {
      // Attempt to add this node as an external haloed node
      unsigned n_ext_haloed_nod =
        external_mesh_pt->nexternal_haloed_node(iproc);
      unsigned external_haloed_node_index =
        external_mesh_pt->add_external_haloed_node_pt(iproc, nod_pt);

      // If it was added then the new index should match the size of the storage
      if (external_haloed_node_index == n_ext_haloed_nod)
      {
        Flat_packed_unsigneds.push_back(1);

#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        std::stringstream junk;
        junk << "Node needs to be constructed [size="
             << Flat_packed_unsigneds.size() << "]; last entry: "
             << Flat_packed_unsigneds[Flat_packed_unsigneds.size() - 1];
        Flat_packed_unsigneds_string.push_back(junk.str());
#endif

        // This helper function gets all the required information for the
        // specified node and stores it into MPI-sendable information
        // so that a halo copy can be made on the receiving process
        get_required_nodal_information_helper(
          iproc, nod_pt, problem_pt, external_mesh_pt, n_cont_inter_values);
      }
      else // It was already added
      {
        Flat_packed_unsigneds.push_back(0);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        std::stringstream junk;
        junk << "Node was already added [size=" << Flat_packed_unsigneds.size()
             << "]; last entry: "
             << Flat_packed_unsigneds[Flat_packed_unsigneds.size() - 1];

        Flat_packed_unsigneds_string.push_back(junk.str());
#endif

        // This node is already an external haloed node, so tell
        // the other process its index in the equivalent external halo storage
        Flat_packed_unsigneds.push_back(external_haloed_node_index);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("external haloed node index");
#endif
      }
    }


    //==========start of add_external_haloed_master_node_helper===============
    /// Helper function to add external haloed node that is a master
    //========================================================================
    void add_external_haloed_master_node_helper(int& iproc,
                                                Node* master_nod_pt,
                                                Problem* problem_pt,
                                                Mesh* const& external_mesh_pt,
                                                int& n_cont_inter_values)
    {
      // Attempt to add node as an external haloed node
      unsigned n_ext_haloed_nod =
        external_mesh_pt->nexternal_haloed_node(iproc);
      unsigned external_haloed_node_index;
      external_haloed_node_index =
        external_mesh_pt->add_external_haloed_node_pt(iproc, master_nod_pt);

      // If it was added the returned index is the same as current storage size
      if (external_haloed_node_index == n_ext_haloed_nod)
      {
        // Indicate that this node needs to be constructed on
        // the other process
        Flat_packed_unsigneds.push_back(1);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back(
          "Node needs to be constructed[2]");
#endif

        // This gets all the required information for the specified
        // master node and stores it into MPI-sendable information
        // so that a halo copy can be made on the receiving process
        get_required_master_nodal_information_helper(iproc,
                                                     master_nod_pt,
                                                     problem_pt,
                                                     external_mesh_pt,
                                                     n_cont_inter_values);
      }
      else // It was already added
      {
        Flat_packed_unsigneds.push_back(0);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Node was already added[2]");
#endif

        // This node is already an external haloed node, so tell
        // the other process its index in the equivalent external halo storage
        Flat_packed_unsigneds.push_back(external_haloed_node_index);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("external haloed node index[2]");
#endif
      }
    }


    //========start of get_required_nodal_information_helper==================
    /// Helper function to get the required nodal information from an
    /// external haloed node so that a fully-functional external halo
    /// node (and therefore element) can be created on the receiving process
    //========================================================================
    void get_required_nodal_information_helper(int& iproc,
                                               Node* nod_pt,
                                               Problem* problem_pt,
                                               Mesh* const& external_mesh_pt,
                                               int& n_cont_inter_values)
    {
      // Tell the halo copy of this node how many values there are
      // [NB this may be different for nodes within the same element, e.g.
      //  when using Lagrange multipliers]
      unsigned n_val = nod_pt->nvalue();
      Flat_packed_unsigneds.push_back(n_val);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
      Flat_packed_unsigneds_string.push_back("Number of values");
#endif

      unsigned n_dim = nod_pt->ndim();
      TimeStepper* time_stepper_pt = nod_pt->time_stepper_pt();

      // Find the timestepper in the list of problem timesteppers
      bool found_timestepper = false;
      unsigned time_stepper_index;
      unsigned n_time_steppers = problem_pt->ntime_stepper();
      for (unsigned i = 0; i < n_time_steppers; i++)
      {
        if (time_stepper_pt == problem_pt->time_stepper_pt(i))
        {
          // Indicate the timestepper's index
          found_timestepper = true;
          time_stepper_index = i;
          break;
        }
      }

      if (found_timestepper)
      {
        Flat_packed_unsigneds.push_back(1);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Found timestepper");
#endif
        Flat_packed_unsigneds.push_back(time_stepper_index);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Timestepper index");
#endif
      }
      else
      {
        Flat_packed_unsigneds.push_back(0);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Not found timestepper");
#endif
      }

      // Default number of previous values to 1
      unsigned n_prev = 1;
      if (time_stepper_pt != 0)
      {
        // Add number of history values to n_prev
        n_prev = time_stepper_pt->ntstorage();
      }

      // Is the node on any boundaries?
      if (nod_pt->is_on_boundary())
      {
        Flat_packed_unsigneds.push_back(1);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Node is on boundary");
#endif

        // Loop over the boundaries of the external mesh
        Vector<unsigned> boundaries;
        unsigned n_bnd = external_mesh_pt->nboundary();
        for (unsigned i_bnd = 0; i_bnd < n_bnd; i_bnd++)
        {
          // Which boundaries (could be more than one) is it on?
          if (nod_pt->is_on_boundary(i_bnd))
          {
            boundaries.push_back(i_bnd);
          }
        }
        unsigned nb = boundaries.size();
        Flat_packed_unsigneds.push_back(nb);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        std::stringstream junk;
        junk << "Node is on " << nb << " boundaries";
        Flat_packed_unsigneds_string.push_back(junk.str());
#endif
        for (unsigned i = 0; i < nb; i++)
        {
          Flat_packed_unsigneds.push_back(boundaries[i]);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          std::stringstream junk;
          junk << "Node is on boundary " << boundaries[i] << " of " << n_bnd;
          Flat_packed_unsigneds_string.push_back(junk.str());
#endif
        }

        // Get pointer to the map of indices associated with
        // additional values created by face elements
        BoundaryNodeBase* bnod_pt = dynamic_cast<BoundaryNodeBase*>(nod_pt);
#ifdef PARANOID
        if (bnod_pt == 0)
        {
          throw OomphLibError("Failed to cast new node to boundary node\n",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        std::map<unsigned, unsigned>* map_pt =
          bnod_pt->index_of_first_value_assigned_by_face_element_pt();

        // No additional values created
        if (map_pt == 0)
        {
          Flat_packed_unsigneds.push_back(0);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          std::stringstream junk;
          Flat_packed_unsigneds_string.push_back(
            "No additional values were created by face element");
#endif
        }
        // Created additional values
        else
        {
          // How many?
          Flat_packed_unsigneds.push_back(map_pt->size());
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          std::stringstream junk;
          junk << "Map size " << map_pt->size() << n_bnd;
          Flat_packed_unsigneds_string.push_back(junk.str());
#endif
          // Loop over entries in map and add to send data
          for (std::map<unsigned, unsigned>::iterator p = map_pt->begin();
               p != map_pt->end();
               p++)
          {
            Flat_packed_unsigneds.push_back((*p).first);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
            std::stringstream junk;
            Flat_packed_unsigneds_string.push_back("Key of map entry");
#endif
            Flat_packed_unsigneds.push_back((*p).second);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
            Flat_packed_unsigneds_string.push_back("Value of map entry");
#endif
          }
        }
      }
      else
      {
        // Not on any boundary
        Flat_packed_unsigneds.push_back(0);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Node is not on any boundary");
#endif
      }

      // Is the Node algebraic?  If so, send its ref values and
      // an indication of its geometric objects if they are stored
      // in the algebraic mesh
      AlgebraicNode* alg_nod_pt = dynamic_cast<AlgebraicNode*>(nod_pt);
      if (alg_nod_pt != 0)
      {
        // The external mesh should be algebraic
        AlgebraicMesh* alg_mesh_pt =
          dynamic_cast<AlgebraicMesh*>(external_mesh_pt);

        // Get default node update function ID
        unsigned update_id = alg_nod_pt->node_update_fct_id();
        Flat_packed_unsigneds.push_back(update_id);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Alg Node update id");
#endif

        // Get reference values at default...
        unsigned n_ref_val = alg_nod_pt->nref_value();
        Flat_packed_unsigneds.push_back(n_ref_val);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Alg Node n ref values");
#endif
        for (unsigned i_ref_val = 0; i_ref_val < n_ref_val; i_ref_val++)
        {
          Flat_packed_doubles.push_back(alg_nod_pt->ref_value(i_ref_val));
        }

        // Access geometric objects at default...
        unsigned n_geom_obj = alg_nod_pt->ngeom_object();
        Flat_packed_unsigneds.push_back(n_geom_obj);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Alg Node n geom objects");
#endif
        for (unsigned i_geom = 0; i_geom < n_geom_obj; i_geom++)
        {
          GeomObject* geom_obj_pt = alg_nod_pt->geom_object_pt(i_geom);

          // Check this against the stored geometric objects in mesh
          unsigned n_geom_list = alg_mesh_pt->ngeom_object_list_pt();

          // Default found index to zero
          unsigned found_geom_object = 0;
          for (unsigned i_list = 0; i_list < n_geom_list; i_list++)
          {
            if (geom_obj_pt == alg_mesh_pt->geom_object_list_pt(i_list))
            {
              found_geom_object = i_list;
            }
          }
          Flat_packed_unsigneds.push_back(found_geom_object);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          Flat_packed_unsigneds_string.push_back("Found geom object");
#endif
        }
      }

      // If it is a MacroElementNodeUpdateNode, everything has been
      // dealt with by the new element already

      // Is it a SolidNode?
      SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
      if (solid_nod_pt != 0)
      {
        unsigned n_solid_val = solid_nod_pt->variable_position_pt()->nvalue();
        for (unsigned i_val = 0; i_val < n_solid_val; i_val++)
        {
          for (unsigned t = 0; t < n_prev; t++)
          {
            Flat_packed_doubles.push_back(
              solid_nod_pt->variable_position_pt()->value(t, i_val));
          }
        }
      }

      // Finally copy info required for all node types
      for (unsigned i_val = 0; i_val < n_val; i_val++)
      {
        for (unsigned t = 0; t < n_prev; t++)
        {
          Flat_packed_doubles.push_back(nod_pt->value(t, i_val));
        }
      }

      // Now do positions
      for (unsigned idim = 0; idim < n_dim; idim++)
      {
        for (unsigned t = 0; t < n_prev; t++)
        {
          Flat_packed_doubles.push_back(nod_pt->x(t, idim));
        }
      }
    }

    //=========start of get_required_master_nodal_information_helper==========
    /// Helper function to get the required master nodal information from an
    /// external haloed master node so that a fully-functional external halo
    /// master node (and possible element) can be created on the receiving
    /// process
    //========================================================================
    void get_required_master_nodal_information_helper(
      int& iproc,
      Node* master_nod_pt,
      Problem* problem_pt,
      Mesh* const& external_mesh_pt,
      int& n_cont_inter_values)
    {
      // Need to send over dimension, position type and number of values
      Flat_packed_unsigneds.push_back(master_nod_pt->ndim());
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
      Flat_packed_unsigneds_string.push_back("Master node ndim");
#endif
      Flat_packed_unsigneds.push_back(master_nod_pt->nposition_type());
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
      Flat_packed_unsigneds_string.push_back("Master node npos_type");
#endif
      Flat_packed_unsigneds.push_back(master_nod_pt->nvalue());
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
      Flat_packed_unsigneds_string.push_back("Master node nvalue");
#endif

      // If it's a solid node, also need to send lagrangian dim and type
      SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(master_nod_pt);
      if (solid_nod_pt != 0)
      {
        Flat_packed_unsigneds.push_back(solid_nod_pt->nlagrangian());
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Master solid node nlagr");
#endif
        Flat_packed_unsigneds.push_back(solid_nod_pt->nlagrangian_type());
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Master solid node nlagr_type");
#endif
      }

      unsigned n_dim = master_nod_pt->ndim();
      TimeStepper* time_stepper_pt = master_nod_pt->time_stepper_pt();

      // Find the timestepper in the list of problem timesteppers
      bool found_timestepper = false;
      unsigned time_stepper_index;
      unsigned n_time_steppers = problem_pt->ntime_stepper();
      for (unsigned i = 0; i < n_time_steppers; i++)
      {
        if (time_stepper_pt == problem_pt->time_stepper_pt(i))
        {
          // Indicate the timestepper's index
          // add 1 to the index so that 0 indicates no timestepper?
          found_timestepper = true;
          time_stepper_index = i;
          break;
        }
      }

      if (found_timestepper)
      {
        Flat_packed_unsigneds.push_back(1);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Master node Found timestepper");
#endif
        Flat_packed_unsigneds.push_back(time_stepper_index);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Master node Timestepper index");
#endif
      }
      else
      {
        Flat_packed_unsigneds.push_back(0);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back(
          "Master node Not found timestepper");
#endif
      }

      // Default number of previous values to 1
      unsigned n_prev = 1;
      if (time_stepper_pt != 0)
      {
        // Add number of history values to n_prev
        n_prev = time_stepper_pt->ntstorage();
      }

      // Is the node on any boundaries?
      if (master_nod_pt->is_on_boundary())
      {
        Flat_packed_unsigneds.push_back(1);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Master node is on boundary");
#endif
        // Loop over the boundaries of the external mesh
        Vector<unsigned> boundaries;
        unsigned n_bnd = external_mesh_pt->nboundary();
        for (unsigned i_bnd = 0; i_bnd < n_bnd; i_bnd++)
        {
          // Which boundaries (could be more than one) is it on?
          if (master_nod_pt->is_on_boundary(i_bnd))
          {
            boundaries.push_back(i_bnd);
          }
        }
        unsigned nb = boundaries.size();
        Flat_packed_unsigneds.push_back(nb);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        std::stringstream junk;
        junk << "Master node is on " << nb << " boundaries";
        Flat_packed_unsigneds_string.push_back(junk.str());
#endif
        for (unsigned i = 0; i < nb; i++)
        {
          Flat_packed_unsigneds.push_back(boundaries[i]);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          std::stringstream junk;
          junk << "Master noode is on boundary " << boundaries[i] << " of "
               << n_bnd;
          Flat_packed_unsigneds_string.push_back(junk.str());
#endif
        }

        // Get pointer to the map of indices associated with
        // additional values created by face elements
        BoundaryNodeBase* bnod_pt =
          dynamic_cast<BoundaryNodeBase*>(master_nod_pt);
#ifdef PARANOID
        if (bnod_pt == 0)
        {
          throw OomphLibError("Failed to cast new node to boundary node\n",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        std::map<unsigned, unsigned>* map_pt =
          bnod_pt->index_of_first_value_assigned_by_face_element_pt();

        // No additional values created
        if (map_pt == 0)
        {
          Flat_packed_unsigneds.push_back(0);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          std::stringstream junk;
          Flat_packed_unsigneds_string.push_back(
            "No additional values were created by face element for this master "
            "node");
#endif
        }
        // Created additional values
        else
        {
          // How many?
          Flat_packed_unsigneds.push_back(map_pt->size());
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          std::stringstream junk;
          junk << "Map size for master node " << map_pt->size() << n_bnd;
          Flat_packed_unsigneds_string.push_back(junk.str());
#endif
          // Loop over entries in map and add to send data
          for (std::map<unsigned, unsigned>::iterator p = map_pt->begin();
               p != map_pt->end();
               p++)
          {
            Flat_packed_unsigneds.push_back((*p).first);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
            std::stringstream junk;
            Flat_packed_unsigneds_string.push_back(
              "Key of map entry for master node");
#endif
            Flat_packed_unsigneds.push_back((*p).second);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
            Flat_packed_unsigneds_string.push_back(
              "Value of map entry for master node");
#endif
          }
        }
      }
      else
      {
        // Not on any boundary
        Flat_packed_unsigneds.push_back(0);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back(
          "Master node is not on any boundary");
#endif
      }

      // Is the Node algebraic?  If so, send its ref values and
      // an indication of its geometric objects if they are stored
      // in the algebraic mesh
      AlgebraicNode* alg_nod_pt = dynamic_cast<AlgebraicNode*>(master_nod_pt);
      if (alg_nod_pt != 0)
      {
        // The external mesh should be algebraic
        AlgebraicMesh* alg_mesh_pt =
          dynamic_cast<AlgebraicMesh*>(external_mesh_pt);

        // Get default node update function ID
        unsigned update_id = alg_nod_pt->node_update_fct_id();
        Flat_packed_unsigneds.push_back(update_id);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Master Alg Node update id");
#endif

        // Get reference values at default...
        unsigned n_ref_val = alg_nod_pt->nref_value();
        Flat_packed_unsigneds.push_back(n_ref_val);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Master Alg Node n ref values");
#endif
        for (unsigned i_ref_val = 0; i_ref_val < n_ref_val; i_ref_val++)
        {
          Flat_packed_doubles.push_back(alg_nod_pt->ref_value(i_ref_val));
        }

        // Access geometric objects at default...
        unsigned n_geom_obj = alg_nod_pt->ngeom_object();
        Flat_packed_unsigneds.push_back(n_geom_obj);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        Flat_packed_unsigneds_string.push_back(
          "Master Alg Node n geom objects");
#endif
        for (unsigned i_geom = 0; i_geom < n_geom_obj; i_geom++)
        {
          GeomObject* geom_obj_pt = alg_nod_pt->geom_object_pt(i_geom);
          // Check this against the stored geometric objects in mesh
          unsigned n_geom_list = alg_mesh_pt->ngeom_object_list_pt();
          // Default found index to zero
          unsigned found_geom_object = 0;
          for (unsigned i_list = 0; i_list < n_geom_list; i_list++)
          {
            if (geom_obj_pt == alg_mesh_pt->geom_object_list_pt(i_list))
            {
              found_geom_object = i_list;
            }
          }
          Flat_packed_unsigneds.push_back(found_geom_object);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          Flat_packed_unsigneds_string.push_back(
            "Master node Found geom object");
#endif
        }
      } // end AlgebraicNode check

      // Is it a MacroElementNodeUpdateNode?
      MacroElementNodeUpdateNode* macro_nod_pt =
        dynamic_cast<MacroElementNodeUpdateNode*>(master_nod_pt);
      if (macro_nod_pt != 0)
      {
        // Loop over current external haloed elements - has the element which
        // controls the node update for this node been added yet?
        GeneralisedElement* macro_node_update_el_pt =
          macro_nod_pt->node_update_element_pt();

        unsigned n_ext_haloed_el =
          external_mesh_pt->nexternal_haloed_element(iproc);
        unsigned external_haloed_el_index;
        external_haloed_el_index =
          external_mesh_pt->add_external_haloed_element_pt(
            iproc, macro_node_update_el_pt);

        // If it wasn't already added, we need to create a halo copy
        if (external_haloed_el_index == n_ext_haloed_el)
        {
          Flat_packed_unsigneds.push_back(1);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          Flat_packed_unsigneds_string.push_back(
            "Master Node needs to be constructed");
#endif
          // Cast to a finite elemnet
          FiniteElement* macro_node_update_finite_el_pt =
            dynamic_cast<FiniteElement*>(macro_node_update_el_pt);

          // We're using macro elements to update...
          MacroElementNodeUpdateMesh* macro_mesh_pt =
            dynamic_cast<MacroElementNodeUpdateMesh*>(external_mesh_pt);
          if (macro_mesh_pt != 0)
          {
            Flat_packed_unsigneds.push_back(1);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
            Flat_packed_unsigneds_string.push_back(
              "Mesh is macro element mesh");
#endif
            // Need to send the macro element number in the mesh across
            MacroElement* macro_el_pt =
              macro_node_update_finite_el_pt->macro_elem_pt();
            unsigned macro_el_num = macro_el_pt->macro_element_number();
            Flat_packed_unsigneds.push_back(macro_el_num);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
            Flat_packed_unsigneds_string.push_back("Number of macro element");
#endif
            // Also need to send
            // the lower left and upper right coordinates of the macro element
            QElementBase* q_el_pt =
              dynamic_cast<QElementBase*>(macro_node_update_el_pt);
            if (q_el_pt != 0)
            {
              // The macro element needs to be set first before
              // its lower left and upper right coordinates can be accessed
              // Now send the lower left and upper right coordinates
              unsigned el_dim = q_el_pt->dim();
              for (unsigned i_dim = 0; i_dim < el_dim; i_dim++)
              {
                Flat_packed_doubles.push_back(q_el_pt->s_macro_ll(i_dim));
                Flat_packed_doubles.push_back(q_el_pt->s_macro_ur(i_dim));
              }
            }
            else // Throw an error
            {
              std::ostringstream error_stream;
              error_stream << "You are using a MacroElement node update\n"
                           << "in a case with non-QElements. This has not\n"
                           << "yet been implemented.\n";
              throw OomphLibError(error_stream.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
          else // Not using macro elements for node update... umm, we're
               // already inside a loop over macro elements, so this
               // should never get here... an error should be thrown I suppose
          {
            Flat_packed_unsigneds.push_back(0);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
            Flat_packed_unsigneds_string.push_back(
              "Mesh is not a macro element mesh");
#endif
          }

          // This element needs to be fully functioning on the other
          // process, so send all the information required to create it
          unsigned n_node = macro_node_update_finite_el_pt->nnode();
          for (unsigned j = 0; j < n_node; j++)
          {
            Node* new_nod_pt = macro_node_update_finite_el_pt->node_pt(j);
            add_external_haloed_node_to_storage(iproc,
                                                new_nod_pt,
                                                problem_pt,
                                                external_mesh_pt,
                                                n_cont_inter_values);
          }
        }
        else // The external haloed element already exists
        {
          Flat_packed_unsigneds.push_back(0);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          Flat_packed_unsigneds_string.push_back(
            "External haloed element already exists");
#endif
          Flat_packed_unsigneds.push_back(external_haloed_el_index);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          Flat_packed_unsigneds_string.push_back(
            "Index of existing external haloed element");
#endif
        }

      } // end of MacroElementNodeUpdateNode check

      // Is it a SolidNode?
      if (solid_nod_pt != 0)
      {
        unsigned n_val = solid_nod_pt->variable_position_pt()->nvalue();
        for (unsigned i_val = 0; i_val < n_val; i_val++)
        {
          for (unsigned t = 0; t < n_prev; t++)
          {
            Flat_packed_doubles.push_back(
              solid_nod_pt->variable_position_pt()->value(t, i_val));
          }
        }
      }

      // Finally copy info required for all node types

      // Halo copy needs to know all the history values
      unsigned n_val = master_nod_pt->nvalue();
      for (unsigned i_val = 0; i_val < n_val; i_val++)
      {
        for (unsigned t = 0; t < n_prev; t++)
        {
          Flat_packed_doubles.push_back(master_nod_pt->value(t, i_val));
        }
      }

      // Now do positions
      for (unsigned idim = 0; idim < n_dim; idim++)
      {
        for (unsigned t = 0; t < n_prev; t++)
        {
          Flat_packed_doubles.push_back(master_nod_pt->x(t, idim));
        }
      }
    }


    //=======start of add_external_halo_node_helper===========================
    /// Helper functiono to add external halo node that is not a master
    //========================================================================
    void add_external_halo_node_helper(Node*& new_nod_pt,
                                       Mesh* const& external_mesh_pt,
                                       unsigned& loc_p,
                                       unsigned& node_index,
                                       FiniteElement* const& new_el_pt,
                                       int& n_cont_inter_values,
                                       Problem* problem_pt)
    {
      // Given the node and the external mesh, and received information
      // about them from process loc_p, construct them on the current process
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
      oomph_info << "Rec:" << Counter_for_flat_packed_unsigneds
                 << " Bool: New node needs to be constructed "
                 << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
                 << std::endl;
#endif
      if (Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++] == 1)
      {
        // Construct a new node based upon sent information
        construct_new_external_halo_node_helper(new_nod_pt,
                                                loc_p,
                                                node_index,
                                                new_el_pt,
                                                external_mesh_pt,
                                                problem_pt);
      }
      else
      {
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        oomph_info << "Rec:" << Counter_for_flat_packed_unsigneds
                   << "  Index of existing external halo node "
                   << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
                   << std::endl;
#endif

        // Copy node from received location
        new_nod_pt = external_mesh_pt->external_halo_node_pt(
          loc_p, Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++]);

        new_el_pt->node_pt(node_index) = new_nod_pt;
      }
    }


    //========start of construct_new_external_halo_node_helper=================
    /// Helper function which constructs a new external halo node (on new
    /// element) with the required information sent from the haloed process
    //========================================================================
    void construct_new_external_halo_node_helper(
      Node*& new_nod_pt,
      unsigned& loc_p,
      unsigned& node_index,
      FiniteElement* const& new_el_pt,
      Mesh* const& external_mesh_pt,
      Problem* problem_pt)
    {
      // The first entry indicates the number of values at this new Node
      // (which may be different across the same element e.g. Lagrange
      // multipliers)
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
      oomph_info << "Rec:" << Counter_for_flat_packed_unsigneds
                 << "  Number of values of external halo node "
                 << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
                 << std::endl;
#endif
      unsigned n_val =
        Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++];

      // Null TimeStepper for now
      TimeStepper* time_stepper_pt = 0;
      // Default number of previous values to 1
      unsigned n_prev = 1;

      // The next entry in Flat_packed_unsigneds indicates
      // if a timestepper is required for this halo node
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
      oomph_info << "Rec:" << Counter_for_flat_packed_unsigneds
                 << "  Timestepper req'd for node "
                 << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
                 << std::endl;
#endif
      if (Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++] == 1)
      {
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        oomph_info << "Rec:" << Counter_for_flat_packed_unsigneds
                   << "  Index of timestepper "
                   << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
                   << std::endl;
#endif
        // Index
        time_stepper_pt = problem_pt->time_stepper_pt(
          Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++]);

        // Check whether number of prev values is "sent" across
        n_prev = time_stepper_pt->ntstorage();
      }

      // If this node was on a boundary then it needs to
      // be on the same boundary here
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
      oomph_info << "Rec:" << Counter_for_flat_packed_unsigneds
                 << "  Is node on boundary? "
                 << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
                 << std::endl;
#endif
      if (Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++] == 1)
      {
        // Construct a new boundary node
        if (time_stepper_pt != 0)
        {
          new_nod_pt =
            new_el_pt->construct_boundary_node(node_index, time_stepper_pt);
        }
        else
        {
          new_nod_pt = new_el_pt->construct_boundary_node(node_index);
        }

        // How many boundaries does the node live on?
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        oomph_info << "Rec:" << Counter_for_flat_packed_unsigneds
                   << " Number of boundaries the node is on: "
                   << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
                   << std::endl;
#endif
        unsigned nb =
          Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++];
        for (unsigned i = 0; i < nb; i++)
        {
          // Boundary number
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          oomph_info << "Rec:" << Counter_for_flat_packed_unsigneds
                     << "  Node is on boundary "
                     << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
                     << std::endl;
#endif
          unsigned i_bnd =
            Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++];
          external_mesh_pt->add_boundary_node(i_bnd, new_nod_pt);
        }

        // Do we have additional values created by face elements?
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        oomph_info << "Rec:" << Counter_for_flat_packed_unsigneds
                   << " Number of additional values created by face element "
                   << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
                   << std::endl;
#endif
        unsigned n_entry =
          Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++];
        if (n_entry > 0)
        {
          // Create storage, if it doesn't already exist, for the map
          // that will contain the position of the first entry of
          // this face element's additional values,
          BoundaryNodeBase* bnew_nod_pt =
            dynamic_cast<BoundaryNodeBase*>(new_nod_pt);
#ifdef PARANOID
          if (bnew_nod_pt == 0)
          {
            throw OomphLibError("Failed to cast new node to boundary node\n",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif
          if (bnew_nod_pt->index_of_first_value_assigned_by_face_element_pt() ==
              0)
          {
            bnew_nod_pt->index_of_first_value_assigned_by_face_element_pt() =
              new std::map<unsigned, unsigned>;
          }

          // Get pointer to the map of indices associated with
          // additional values created by face elements
          std::map<unsigned, unsigned>* map_pt =
            bnew_nod_pt->index_of_first_value_assigned_by_face_element_pt();

          // Loop over number of entries in map
          for (unsigned i = 0; i < n_entry; i++)
          {
            // Read out pairs...

#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
            oomph_info
              << "Rec:" << Counter_for_flat_packed_unsigneds
              << " Key of map entry"
              << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
              << std::endl;
#endif
            unsigned first =
              Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++];

#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
            oomph_info
              << "Rec:" << Counter_for_flat_packed_unsigneds
              << " Value of map entry"
              << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
              << std::endl;
#endif
            unsigned second =
              Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++];

            // ...and assign
            (*map_pt)[first] = second;
          }
        }
      }
      else
      {
        // Construct an ordinary (non-boundary) node
        if (time_stepper_pt != 0)
        {
          new_nod_pt = new_el_pt->construct_node(node_index, time_stepper_pt);
        }
        else
        {
          new_nod_pt = new_el_pt->construct_node(node_index);
        }
      }

      // Node constructed: add to external halo nodes
      external_mesh_pt->add_external_halo_node_pt(loc_p, new_nod_pt);

      // Is the new constructed node Algebraic?
      AlgebraicNode* new_alg_nod_pt = dynamic_cast<AlgebraicNode*>(new_nod_pt);

      // If it is algebraic, its node update functions will
      // not yet have been set up properly
      if (new_alg_nod_pt != 0)
      {
        // The AlgebraicMesh is the external mesh
        AlgebraicMesh* alg_mesh_pt =
          dynamic_cast<AlgebraicMesh*>(external_mesh_pt);

        /// The first entry of All_alg_nodal_info contains
        /// the default node update id
        /// e.g. for the quarter circle there are
        /// "Upper_left_box", "Lower right box" etc...
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        oomph_info << "Rec:" << Counter_for_flat_packed_unsigneds
                   << "  Alg node update id "
                   << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
                   << std::endl;
#endif

        unsigned update_id =
          Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++];

        Vector<double> ref_value;

        // The size of this vector is in the next entry
        // of All_alg_nodal_info
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        oomph_info << "Rec:" << Counter_for_flat_packed_unsigneds
                   << "  Alg node # of ref values "
                   << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
                   << std::endl;
#endif
        unsigned n_ref_val =
          Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++];

        // The reference values themselves are in
        // All_alg_ref_value
        ref_value.resize(n_ref_val);
        for (unsigned i_ref = 0; i_ref < n_ref_val; i_ref++)
        {
          ref_value[i_ref] =
            Flat_packed_doubles[Counter_for_flat_packed_doubles++];
        }

        Vector<GeomObject*> geom_object_pt;
        /// again we need the size of this vector as it varies
        /// between meshes; we also need some indication
        /// as to which geometric object should be used...

        // The size of this vector is in the next entry
        // of All_alg_nodal_info
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
        oomph_info << "Rec:" << Counter_for_flat_packed_unsigneds
                   << "  Alg node # of geom objects "
                   << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
                   << std::endl;
#endif
        unsigned n_geom_obj =
          Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++];

        // The remaining indices are in the rest of
        // All_alg_nodal_info
        geom_object_pt.resize(n_geom_obj);
        for (unsigned i_geom = 0; i_geom < n_geom_obj; i_geom++)
        {
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
          oomph_info << "Rec:" << Counter_for_flat_packed_unsigneds
                     << "  Alg node: geom object index "
                     << Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]
                     << std::endl;
#endif
          unsigned geom_index =
            Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds++];
          // This index indicates which of the AlgebraicMesh's
          // stored geometric objects should be used
          // (0 is a null pointer; everything else should have
          //  been filled in by the specific Mesh).  If it
          // hasn't been filled in then the update_node_update
          // call should fix it
          geom_object_pt[i_geom] = alg_mesh_pt->geom_object_list_pt(geom_index);
        }

        /// For the received update_id, ref_value, geom_object
        /// call add_node_update_info
        new_alg_nod_pt->add_node_update_info(
          update_id, alg_mesh_pt, geom_object_pt, ref_value);

        /// Now call update_node_update
        alg_mesh_pt->update_node_update(new_alg_nod_pt);
      }

      // Is the node a MacroElementNodeUpdateNode?
      MacroElementNodeUpdateNode* macro_nod_pt =
        dynamic_cast<MacroElementNodeUpdateNode*>(new_nod_pt);

      if (macro_nod_pt != 0)
      {
        // Need to call set_node_update_info; this requires
        // a Vector<GeomObject*> (taken from the mesh)
        Vector<GeomObject*> geom_object_vector_pt;

        // Access the required geom objects from the
        // MacroElementNodeUpdateMesh
        MacroElementNodeUpdateMesh* macro_mesh_pt =
          dynamic_cast<MacroElementNodeUpdateMesh*>(external_mesh_pt);
        geom_object_vector_pt = macro_mesh_pt->geom_object_vector_pt();

        // Get local coordinate of node in new element
        Vector<double> s_in_macro_node_update_element;
        new_el_pt->local_coordinate_of_node(node_index,
                                            s_in_macro_node_update_element);

        // Set node update info for this node
        macro_nod_pt->set_node_update_info(
          new_el_pt, s_in_macro_node_update_element, geom_object_vector_pt);
      }

      // Is the new node a SolidNode?
      SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(new_nod_pt);
      if (solid_nod_pt != 0)
      {
        unsigned n_solid_val = solid_nod_pt->variable_position_pt()->nvalue();
        for (unsigned i_val = 0; i_val < n_solid_val; i_val++)
        {
          for (unsigned t = 0; t < n_prev; t++)
          {
            solid_nod_pt->variable_position_pt()->set_value(
              t, i_val, Flat_packed_doubles[Counter_for_flat_packed_doubles++]);
          }
        }
      }

      // If there are additional values, resize the node
      unsigned n_new_val = new_nod_pt->nvalue();
      if (n_val > n_new_val)
      {
        new_nod_pt->resize(n_val);
      }

      // Get copied history values
      //  unsigned n_val=new_nod_pt->nvalue();
      for (unsigned i_val = 0; i_val < n_val; i_val++)
      {
        for (unsigned t = 0; t < n_prev; t++)
        {
          new_nod_pt->set_value(
            t, i_val, Flat_packed_doubles[Counter_for_flat_packed_doubles++]);
        }
      }

      // Get copied history values for positions
      unsigned n_dim = new_nod_pt->ndim();
      for (unsigned idim = 0; idim < n_dim; idim++)
      {
        for (unsigned t = 0; t < n_prev; t++)
        {
          // Copy to coordinate
          new_nod_pt->x(t, idim) =
            Flat_packed_doubles[Counter_for_flat_packed_doubles++];
        }
      }
    }


    //=====================================================================
    /// Locate zeta for current set of missing coordinates; vector-based version
    //=====================================================================
    void locate_zeta_for_missing_coordinates(
      int& iproc,
      Mesh* const& external_mesh_pt,
      Problem* problem_pt,
      Vector<MeshAsGeomObject*>& mesh_geom_obj_pt)
    {
      // How many meshes are we dealing with?
      unsigned n_mesh = mesh_geom_obj_pt.size();

      // Storage for number of processors, current process and communicator
      OomphCommunicator* comm_pt = problem_pt->communicator_pt();
      int n_proc = comm_pt->nproc();
      int my_rank = comm_pt->my_rank();

      // Clear vectors containing data to be sent
      Flat_packed_doubles.resize(0);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
      Flat_packed_unsigneds_string.resize(0);
#endif
      Flat_packed_unsigneds.resize(0);
      Flat_packed_located_coordinates.resize(0);

      // Flush storage for zetas not found locally (when
      // processing the zeta coordinates received from "previous"
      // processor)
      Flat_packed_zetas_not_found_locally.resize(0);

      // Number of zeta tuples to be dealt with (includes padding!)
      unsigned n_zeta = Received_flat_packed_zetas_to_be_found.size() / Dim;

      // Create storage for the processor id (plus one) on which
      // the zetas stored in Flat_packed_zetas_not_found_locally[...]
      // were located. (Remains zero for padded entries).
      Proc_id_plus_one_of_external_element.resize(n_zeta, 0);

      // Create storage for the status of the (external halo) element associated
      // the zetas stored in Flat_packed_zetas_not_found_locally[...].
      // It either hasn't been found, already exists on the processor
      // that needs it, or needs to be newly created. (Remains Not_found
      // for padded entries).
      Located_element_status.resize(n_zeta, Not_found);

      // Counter for flat-packed array of external zeta coordinates
      unsigned count = 0;

      // Current mesh
      unsigned i_mesh = 0;

      // Loop over the zeta tuples that we received from elsewhere and
      // are trying to find here for current mesh
      for (unsigned i = 0; i < n_zeta; i++)
      {
        // Storage for global coordinates to be located
        Vector<double> x_global(Dim);

        // Loop to fill in coordinates
        for (unsigned ii = 0; ii < Dim; ii++)
        {
          x_global[ii] = Received_flat_packed_zetas_to_be_found[count];
          count++;
        }

        // Check if we've reached the end of the mesh
        bool reached_end_of_mesh = false;
        unsigned dbl_max_count = 0;
        for (unsigned ii = 0; ii < Dim; ii++)
        {
          if (x_global[ii] == DBL_MAX)
          {
            dbl_max_count++;
            reached_end_of_mesh = true;
          }
        }

        // Reached end of mesh
        if (reached_end_of_mesh)
        {
#ifdef PARANOID
          // Check if all coordinates were set to DBX_MAX
          if (dbl_max_count != Dim)
          {
            std::ostringstream error_stream;
            error_stream << "Appear to have reached end of mesh " << i_mesh
                         << " but only " << dbl_max_count << " out of " << Dim
                         << " zeta coordinates have been set to DBX_MAX\n";
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif
          // Indicate end of mesh in flat packed data
          for (unsigned ii = 0; ii < Dim; ii++)
          {
            Flat_packed_zetas_not_found_locally.push_back(DBL_MAX);
          }

          // Bump mesh counter
          i_mesh++;

          // Bail out if we're done
          if (i_mesh == n_mesh)
          {
            return;
          }
        }

        // Perform locate_zeta for these coordinates and current mesh
        GeomObject* sub_geom_obj_pt = 0;
        Vector<double> ss(Dim);
        if (!reached_end_of_mesh)
        {
          mesh_geom_obj_pt[i_mesh]->locate_zeta(x_global, sub_geom_obj_pt, ss);

          // Did the locate method work?
          if (sub_geom_obj_pt != 0)
          {
            // Get the source element - bulk or not?
            GeneralisedElement* source_el_pt = 0;
            if (!Use_bulk_element_as_external)
            {
              source_el_pt = dynamic_cast<FiniteElement*>(sub_geom_obj_pt);
            }
            else
            {
              FaceElement* face_el_pt =
                dynamic_cast<FaceElement*>(sub_geom_obj_pt);
              source_el_pt =
                dynamic_cast<FiniteElement*>(face_el_pt->bulk_element_pt());
            }

            // Check if the returned element is halo
            if (!source_el_pt->is_halo()) // cannot accept halo here
            {
              // The correct non-halo element has been located; this will become
              // an external haloed element on the current process, and an
              // external halo copy needs to be created on the current process
              // minus wherever we are in the "ring-loop"
              int halo_copy_proc = my_rank - iproc;

              // If iproc is bigger than my_rank then we've "gone through"
              // nproc-1
              if (my_rank < iproc)
              {
                halo_copy_proc = n_proc + halo_copy_proc;
              }

              // So, we found zeta on the current processor
              Proc_id_plus_one_of_external_element[i] = my_rank + 1;

              // This source element is an external halo on process
              // halo_copy_proc but it should only be added to the storage if it
              // hasn't been added already, and this information also needs to
              // be communicated over to the other process

              unsigned n_extern_haloed =
                external_mesh_pt->nexternal_haloed_element(halo_copy_proc);
              unsigned external_haloed_el_index =
                external_mesh_pt->add_external_haloed_element_pt(halo_copy_proc,
                                                                 source_el_pt);

              // If it was added to the storage then the returned index
              // will be the same as the (old) size of the storage
              if (external_haloed_el_index == n_extern_haloed)
              {
                // Set Located_element_status to say it
                // should be newly created
                Located_element_status[i] = New;

                // How many continuously interpolated values are there?
                int n_cont_inter_values = -1;
                if (dynamic_cast<RefineableElement*>(source_el_pt) != 0)
                {
                  n_cont_inter_values =
                    dynamic_cast<RefineableElement*>(source_el_pt)
                      ->ncont_interpolated_values();
                }

                // Since it is (externally) haloed from the current process,
                // the info required to create a new element in the equivalent
                // external halo layer on process halo_copy_proc needs to be
                // sent there

                // If we're using macro elements to update...
                MacroElementNodeUpdateMesh* macro_mesh_pt =
                  dynamic_cast<MacroElementNodeUpdateMesh*>(external_mesh_pt);
                if (macro_mesh_pt != 0)
                {
                  Flat_packed_unsigneds.push_back(1);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
                  Flat_packed_unsigneds_string.push_back(
                    "Mesh is macro element mesh[2]");
#endif
                  // Cast to finite element... this must work because it's
                  // a macroelement no update mesh
                  FiniteElement* source_finite_el_pt =
                    dynamic_cast<FiniteElement*>(source_el_pt);

                  MacroElement* macro_el_pt =
                    source_finite_el_pt->macro_elem_pt();
                  // Send the macro element number across
                  unsigned macro_el_num = macro_el_pt->macro_element_number();
                  Flat_packed_unsigneds.push_back(macro_el_num);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
                  Flat_packed_unsigneds_string.push_back(
                    "Number of macro element[2]");
#endif
                  // we need to send
                  // the lower left and upper right coordinates of the macro
                  QElementBase* q_el_pt =
                    dynamic_cast<QElementBase*>(source_el_pt);
                  if (q_el_pt != 0)
                  {
                    // The macro element needs to be set first before
                    // its lower left and upper right coordinates can be
                    // accessed Now send the lower left and upper right
                    // coordinates
                    unsigned el_dim = q_el_pt->dim();
                    for (unsigned i_dim = 0; i_dim < el_dim; i_dim++)
                    {
                      Flat_packed_doubles.push_back(q_el_pt->s_macro_ll(i_dim));
                      Flat_packed_doubles.push_back(q_el_pt->s_macro_ur(i_dim));
                    }
                  }
                  else // Throw an error
                  {
                    std::ostringstream error_stream;
                    error_stream
                      << "You are using a MacroElement node update\n"
                      << "in a case with non-QElements. This has not\n"
                      << "yet been implemented.\n";
                    throw OomphLibError(error_stream.str(),
                                        OOMPH_CURRENT_FUNCTION,
                                        OOMPH_EXCEPTION_LOCATION);
                  }
                }
                else // Not using macro elements to update
                {
                  Flat_packed_unsigneds.push_back(0);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
                  Flat_packed_unsigneds_string.push_back(
                    "Mesh is not a macro element mesh [2]");
#endif
                }


                // Cast to finite element... this must work because it's
                // a macroelement no update mesh
                FiniteElement* source_finite_el_pt =
                  dynamic_cast<FiniteElement*>(source_el_pt);
#ifdef PARANOID
                if (source_finite_el_pt == 0)
                {
                  throw OomphLibError(
                    "Unable to cast source function to finite element\n",
                    "Multi_domain_functions::locate_zeta_for_missing_"
                    "coordinates()",
                    OOMPH_EXCEPTION_LOCATION);
                }
#endif


                // Loop over the nodes of the new source element
                unsigned n_node = source_finite_el_pt->nnode();
                for (unsigned j = 0; j < n_node; j++)
                {
                  Node* nod_pt = source_finite_el_pt->node_pt(j);

                  // Add the node to the storage; this routine
                  // also takes care of any master nodes if the
                  // node is hanging
                  add_external_haloed_node_to_storage(halo_copy_proc,
                                                      nod_pt,
                                                      problem_pt,
                                                      external_mesh_pt,
                                                      n_cont_inter_values);
                }
              }
              else // it has already been added, so tell the other process
              {
                // Set Located_element_status to indicate an element has
                // already been added
                Located_element_status[i] = Exists;
                Flat_packed_unsigneds.push_back(external_haloed_el_index);
#ifdef ANNOTATE_MULTI_DOMAIN_COMMUNICATION
                Flat_packed_unsigneds_string.push_back(
                  "Index of existing external haloed element[2]");
#endif
              }

              // The coordinates returned by locate_zeta are also needed
              // in the setup of the source elements on the other process
              if (!Use_bulk_element_as_external)
              {
                for (unsigned ii = 0; ii < Dim; ii++)
                {
                  Flat_packed_located_coordinates.push_back(ss[ii]);
                }
              }
              else // translate the coordinates to the bulk element
              {
                // The translation is from Lagrangian to Eulerian
                FaceElement* face_el_pt =
                  dynamic_cast<FaceElement*>(sub_geom_obj_pt);
                // Get the dimension of the BulkElement
                unsigned bulk_el_dim =
                  dynamic_cast<FiniteElement*>(source_el_pt)->dim();
                Vector<double> s_trans(bulk_el_dim);
                face_el_pt->get_local_coordinate_in_bulk(ss, s_trans);
                for (unsigned ii = 0; ii < bulk_el_dim; ii++)
                {
                  Flat_packed_located_coordinates.push_back(s_trans[ii]);
                }
              }
            }
            else // halo, so search again until non-halo equivalent is located
            {
              // Add required information to arrays (as below)
              for (unsigned ii = 0; ii < Dim; ii++)
              {
                Flat_packed_zetas_not_found_locally.push_back(x_global[ii]);
              }
              // It wasn't found here
              Proc_id_plus_one_of_external_element[i] = 0;

              // Set Located_element_status to indicate not found
              Located_element_status[i] = Not_found;
            }
          }
          else // not successful this time (i.e. sub_geom_obj_pt==0), so
          // prepare for next process to try
          {
            // Add this global coordinate to the LOCAL zeta array
            for (unsigned ii = 0; ii < Dim; ii++)
            {
              Flat_packed_zetas_not_found_locally.push_back(x_global[ii]);
            }
            // It wasn't found here
            Proc_id_plus_one_of_external_element[i] = 0;

            // Set Located_element_status to indicate not found
            Located_element_status[i] = Not_found;
          }

        } // end of mesh not reached

      } // end of loop over flat-packed zeta tuples
    }


#endif


    //=====================================================================
    /// locate zeta for current set of "local" coordinates
    /// vector-based version
    //=====================================================================
    void locate_zeta_for_local_coordinates(
      const Vector<Mesh*>& mesh_pt,
      Mesh* const& external_mesh_pt,
      Vector<MeshAsGeomObject*>& mesh_geom_obj_pt,
      const unsigned& interaction_index)
    {
      // Flush storage for zetas not found locally
      Flat_packed_zetas_not_found_locally.resize(0);

      // Number of meshes
      unsigned n_mesh = mesh_pt.size();

#ifdef PARANOID
      if (mesh_geom_obj_pt.size() != n_mesh)
      {
        std::ostringstream error_stream;
        error_stream << "Sizes of mesh_geom_obj_pt [ "
                     << mesh_geom_obj_pt.size() << " ] and "
                     << "mesh_pt [ " << n_mesh << " ] don't match.\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Element counter
      unsigned e_count = 0;

      // Loop over meshes
      for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++)
      {
        // Number of local elements
        unsigned n_element = mesh_pt[i_mesh]->nelement();

        // Loop over this processor's elements
        for (unsigned e = 0; e < n_element; e++)
        {
          ElementWithExternalElement* el_pt =
            dynamic_cast<ElementWithExternalElement*>(
              mesh_pt[i_mesh]->element_pt(e));
#ifdef OOMPH_HAS_MPI
          // Only visit non-halo elements -- we're not setting up external
          // elements for on-halos!
          if (!el_pt->is_halo())
#endif
          {
            // Find number of Gauss points and element dimension
            unsigned n_intpt = el_pt->integral_pt()->nweight();
            unsigned el_dim = el_pt->dim();


#ifdef PARANOID
            if (el_dim != Dim)
            {
              std::ostringstream error_stream;
              error_stream << "Dimension of element " << el_dim
                           << " is not consitent with dimension assumed \n"
                           << " in multidomain namespace, " << Dim << std::endl;
              throw OomphLibError(error_stream.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
#endif

            // Set storage for local and global coordinates
            Vector<double> s_local(el_dim);
            Vector<double> x_global(el_dim);

            // Loop over integration points
            for (unsigned ipt = 0; ipt < n_intpt; ipt++)
            {
              // Has this integration point been done yet?
              if (External_element_located[e_count][ipt] == 0)
              {
                // Get local coordinates
                for (unsigned i = 0; i < el_dim; i++)
                {
                  s_local[i] = el_pt->integral_pt()->knot(ipt, i);
                }
                // Interpolate to global coordinates
                el_pt->interpolated_zeta(s_local, x_global);

                // Storage for geometric object and its local coordinates
                GeomObject* sub_geom_obj_pt = 0;
                Vector<double> s_ext(el_dim);
                mesh_geom_obj_pt[i_mesh]->locate_zeta(
                  x_global, sub_geom_obj_pt, s_ext);

                // Has the required element been located?
                if (sub_geom_obj_pt != 0)
                {
                  // The required element has been located
                  // The located coordinates have the same dimension as the bulk
                  GeneralisedElement* source_el_pt;
                  Vector<double> s_source(el_dim);

                  // Is the bulk element the actual external element?
                  if (!Use_bulk_element_as_external)
                  {
                    // Use the object directly (it must be a finite element)
                    source_el_pt =
                      dynamic_cast<FiniteElement*>(sub_geom_obj_pt);
                    s_source = s_ext;
                  }
                  else
                  {
                    // Cast to a FaceElement and use the bulk element
                    FaceElement* face_el_pt =
                      dynamic_cast<FaceElement*>(sub_geom_obj_pt);
                    source_el_pt = face_el_pt->bulk_element_pt();

                    // Need to resize the located coordinates to have the same
                    // dimension as the bulk element
                    s_source.resize(
                      dynamic_cast<FiniteElement*>(source_el_pt)->dim());

                    // Translate the returned local coords into the bulk element
                    face_el_pt->get_local_coordinate_in_bulk(s_ext, s_source);
                  }

                  // Check if it's a halo; if it is then the non-halo equivalent
                  // needs to be located from another processor (unless we
                  // accept halo elements as external elements)
#ifdef OOMPH_HAS_MPI
                  if (Allow_use_of_halo_elements_as_external_elements ||
                      (!source_el_pt->is_halo()))
#endif
                  {
                    // Need to cast to a FiniteElement
                    FiniteElement* source_finite_el_pt =
                      dynamic_cast<FiniteElement*>(source_el_pt);

                    // Set the external element pointer and local coordinates
                    el_pt->external_element_pt(interaction_index, ipt) =
                      source_finite_el_pt;
                    el_pt->external_element_local_coord(interaction_index,
                                                        ipt) = s_source;

                    // Set the lookup array to 1/true
                    External_element_located[e_count][ipt] = 1;
                  }
#ifdef OOMPH_HAS_MPI
                  // located element is halo and we're not accepting haloes
                  // obviously only makes sense in mpi mode...
                  else
                  {
                    // Add required information to arrays
                    for (unsigned i = 0; i < el_dim; i++)
                    {
                      Flat_packed_zetas_not_found_locally.push_back(
                        x_global[i]);
                    }
                  }
#endif
                }
                else
                {
                  // Search has failed then add the required information to the
                  // arrays which need to be sent to the other processors so
                  // that they can perform the locate_zeta

                  // Add this global coordinate to the LOCAL zeta array
                  for (unsigned i = 0; i < el_dim; i++)
                  {
                    Flat_packed_zetas_not_found_locally.push_back(x_global[i]);
                  }
                }
              }
            } // end loop over integration points
          } // end for halo

          // Bump up counter for all elements
          e_count++;

        } // end loop over local elements

        // Mark end of mesh data in flat packed array
        for (unsigned i = 0; i < Dim; i++)
        {
          Flat_packed_zetas_not_found_locally.push_back(DBL_MAX);
        }

      } // end of loop over meshes
    }


    //=====================================================================
    /// Helper function that computes the dimension of the elements within
    /// each of the specified meshes (and checks they are the same).
    /// Stores result in Dim.
    //=====================================================================
    void get_dim_helper(Problem* problem_pt,
                        Mesh* const& mesh_pt,
                        Mesh* const& external_mesh_pt)
    {
#ifdef OOMPH_HAS_MPI
      // Storage for number of processors, current process and communicator
      OomphCommunicator* comm_pt = problem_pt->communicator_pt();
#endif

      // Extract the element dimensions from the first element of each mesh
      unsigned mesh_dim = 0;
      if (mesh_pt->nelement() > 0)
      {
        mesh_dim = dynamic_cast<FiniteElement*>(mesh_pt->element_pt(0))->dim();
      }
      unsigned external_mesh_dim = 0;
      if (external_mesh_pt->nelement() > 0)
      {
        external_mesh_dim =
          dynamic_cast<FiniteElement*>(external_mesh_pt->element_pt(0))->dim();
      }

      // Need to do an Allreduce
#ifdef OOMPH_HAS_MPI
      int n_proc = comm_pt->nproc();
      if (n_proc > 1)
      {
        unsigned mesh_dim_reduce;
        MPI_Allreduce(&mesh_dim,
                      &mesh_dim_reduce,
                      1,
                      MPI_UNSIGNED,
                      MPI_MAX,
                      comm_pt->mpi_comm());
        mesh_dim = mesh_dim_reduce;

        unsigned external_mesh_dim_reduce;
        MPI_Allreduce(&external_mesh_dim,
                      &external_mesh_dim_reduce,
                      1,
                      MPI_UNSIGNED,
                      MPI_MAX,
                      comm_pt->mpi_comm());
        external_mesh_dim = external_mesh_dim_reduce;
      }
#endif

      // Check the dimensions are the same!
      if (mesh_dim != external_mesh_dim)
      {
        std::ostringstream error_stream;
        error_stream << "The elements within the two meshes do not\n"
                     << "have the same dimension, so the multi-domain\n"
                     << "method will not work.\n"
                     << "For the mesh, dim=" << mesh_dim
                     << ", and the external mesh, dim=" << external_mesh_dim
                     << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Set dimension
      Dim = mesh_dim;
    }


    //=================================================================
    /// Helper function that clears all the information used
    /// during the external storage creation
    //=================================================================
    void clean_up()
    {
      Flat_packed_zetas_not_found_locally.clear();
      Received_flat_packed_zetas_to_be_found.clear();
      Proc_id_plus_one_of_external_element.clear();
      Located_element_status.clear();
      Flat_packed_located_coordinates.clear();
      Flat_packed_doubles.clear();
      Flat_packed_unsigneds.clear();
      External_element_located.clear();
    }

    /// Vector of zeta coordinates that we're currently trying to locate;
    /// used in sorting of bin entries in further_away() comparison function
    Vector<double> Zeta_coords_for_further_away_comparison;

    /// Comparison function for sorting entries in bin: Returns true if
    /// point identified by p1 (comprising pointer to finite element and
    /// vector of local coordinates within that element) is closer to
    /// Zeta_coords_for_further_away_comparison than p2
    bool first_closer_than_second(
      const std::pair<FiniteElement*, Vector<double>>& p1,
      const std::pair<FiniteElement*, Vector<double>>& p2)
    {
      // First point
      FiniteElement* el_pt = p1.first;
      Vector<double> s(p1.second);
      Vector<double> zeta(Dim);
      el_pt->interpolated_zeta(s, zeta);
      double dist_squared1 = 0.0;
      for (unsigned i = 0; i < Dim; i++)
      {
        dist_squared1 +=
          (zeta[i] - Zeta_coords_for_further_away_comparison[i]) *
          (zeta[i] - Zeta_coords_for_further_away_comparison[i]);
      }

      // Second point
      el_pt = p2.first;
      s = p2.second;
      el_pt->interpolated_zeta(s, zeta);
      double dist_squared2 = 0.0;
      for (unsigned i = 0; i < Dim; i++)
      {
        dist_squared2 +=
          (zeta[i] - Zeta_coords_for_further_away_comparison[i]) *
          (zeta[i] - Zeta_coords_for_further_away_comparison[i]);
      }

      // Which one is further
      if (dist_squared1 < dist_squared2)
      {
        return true;
      }
      else
      {
        return false;
      }
    }
  } // namespace Multi_domain_functions

} // namespace oomph
