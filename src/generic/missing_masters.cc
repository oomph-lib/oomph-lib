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
// Non-templated missing masters functions which help to reconcile the hanging
// status of nodes in the halo layer of distributed meshes

// oomph-lib header
#include "missing_masters.h"
#include "missing_masters.template.cc"
#include "mesh.h"
#include "algebraic_elements.h"
#include "macro_element_node_update_element.h"
#include "Qelements.h"

namespace oomph
{
  //======================================================================
  // Namespace for "global" missing-master-locating functions
  //======================================================================
  namespace Missing_masters_functions
  {
    // Workspace for locate zeta methods
    //----------------------------------

#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION

    // Temporary vector of strings to enable full annotation of multi domain
    // comms (but keep alive because it would be such a bloody pain to
    // rewrite it if things ever go wrong again...)
    // This is left over from the multi-domain stuff and should work
    // in the same way, but it has not been tested.
    Vector<std::string> Flat_packed_unsigneds_string;

#endif

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

    //========start of add_external_haloed_node_to_storage====================
    /// Helper function to add external haloed nodes, including any masters
    //========================================================================
    void add_external_haloed_node_to_storage(int& iproc,
                                             Node* nod_pt,
                                             Mesh* const& mesh_pt,
                                             int& n_cont_inter_values,
                                             Vector<unsigned>& send_unsigneds,
                                             Vector<double>& send_doubles)
    {
      // Add the node if required
      add_external_haloed_node_helper(iproc,
                                      nod_pt,
                                      mesh_pt,
                                      n_cont_inter_values,
                                      send_unsigneds,
                                      send_doubles);

      // Recursively add any master nodes (and their master nodes etc)
      recursively_add_masters_of_external_haloed_node(iproc,
                                                      nod_pt,
                                                      mesh_pt,
                                                      n_cont_inter_values,
                                                      send_unsigneds,
                                                      send_doubles);
    }


    //========================================================================
    /// Recursively add any master nodes (and their master nodes etc) of
    /// external nodes
    //========================================================================
    void recursively_add_masters_of_external_haloed_node(
      int& iproc,
      Node* nod_pt,
      Mesh* const& mesh_pt,
      int& n_cont_inter_values,
      Vector<unsigned>& send_unsigneds,
      Vector<double>& send_doubles)
    {
      // Loop over continuously interpolated values and add masters
      for (int i_cont = -1; i_cont < n_cont_inter_values; i_cont++)
      {
        if (nod_pt->is_hanging(i_cont))
        {
          // Indicate that this node is a hanging node so the other
          // process knows to create HangInfo and masters, etc.
          send_unsigneds.push_back(1);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          Flat_packed_unsigneds_string.push_back("Is hanging");
#endif
          // If this is a hanging node then add all its masters as
          // external halo nodes if they have not yet been added
          HangInfo* hang_pt = nod_pt->hanging_pt(i_cont);
          // Loop over masters
          unsigned n_master = hang_pt->nmaster();

          // Indicate number of master nodes to add on other process
          send_unsigneds.push_back(n_master);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          Flat_packed_unsigneds_string.push_back("nmaster");
#endif
          for (unsigned m = 0; m < n_master; m++)
          {
            Node* master_nod_pt = hang_pt->master_node_pt(m);

            // Call the helper function for master nodes
            add_external_haloed_master_node_helper(iproc,
                                                   master_nod_pt,
                                                   mesh_pt,
                                                   n_cont_inter_values,
                                                   send_unsigneds,
                                                   send_doubles);

            // Indicate the weight of this master
            send_doubles.push_back(hang_pt->master_weight(m));

            // Recursively add any master nodes (and their master nodes etc)
            recursively_add_masters_of_external_haloed_node(iproc,
                                                            master_nod_pt,
                                                            mesh_pt,
                                                            n_cont_inter_values,
                                                            send_unsigneds,
                                                            send_doubles);
          }
        }
        else
        {
          // Indicate that it's not a hanging node in this variable
          send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
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
                                         Mesh* const& mesh_pt,
                                         int& n_cont_inter_values,
                                         Vector<unsigned>& send_unsigneds,
                                         Vector<double>& send_doubles)
    {
      // Check to see if this haloed node already exists in (internal)
      // haloed node storage with any processor
      bool found_internally = false;
      unsigned shared_node_index = 0;

      // Get vector of all shared nodes with processor iproc
      Vector<Node*> shared_node_pt;
      mesh_pt->get_shared_node_pt(iproc, shared_node_pt);

      // Search the internal haloed storage for this node
      Vector<Node*>::iterator it =
        std::find(shared_node_pt.begin(), shared_node_pt.end(), nod_pt);

      // Check if the node was found in shared storage
      if (it != shared_node_pt.end())
      {
        // Node found in (internal) haloed storage
        found_internally = true;
        // Store the index in this storage
        shared_node_index = it - shared_node_pt.begin();
      }

      /// / Slow search version without additional access function in Mesh class
      /// /Search the internal shared node storage for this node
      // for(unsigned i=0; i<mesh_pt->nshared_node(iproc); i++)
      // {
      //  if(nod_pt == mesh_pt->shared_node_pt(iproc,i))
      //   {
      //    //Node found in (internal) shared storage
      //    found_internally = true;
      //    shared_node_index = i;
      //    break;
      //   }
      // }

      // If we've found the node internally
      if (found_internally)
      {
        // Indicate that this node doesn not need to be constructed on
        // the other process
        send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        std::stringstream junk;
        junk << "Node was already added [size=" << send_unsigneds.size()
             << "]; last entry: " << send_unsigneds[send_unsigneds.size() - 1];

        Flat_packed_unsigneds_string.push_back(junk.str());
#endif

        // This node is already shared with processor iproc, so tell the other
        // processor its index in the shared node storage
        send_unsigneds.push_back(1);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("haloed node found internally");
#endif
        send_unsigneds.push_back(shared_node_index);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("(internal) haloed node index");
#endif
      }
      else
      {
        // Attempt to add this node as an external haloed node
        unsigned n_ext_haloed_nod = mesh_pt->nexternal_haloed_node(iproc);
        unsigned external_haloed_node_index;
        external_haloed_node_index =
          mesh_pt->add_external_haloed_node_pt(iproc, nod_pt);

        // If it was added then the new index should match the size of the
        // storage
        if (external_haloed_node_index == n_ext_haloed_nod)
        {
          // Indicate that this node needs to be constructed on
          // the other process
          send_unsigneds.push_back(1);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          std::stringstream junk;
          junk << "Node needs to be constructed [size=" << send_unsigneds.size()
               << "]; last entry: "
               << send_unsigneds[send_unsigneds.size() - 1];
          Flat_packed_unsigneds_string.push_back(junk.str());
#endif

          // This helper function gets all the required information for the
          // specified node and stores it into MPI-sendable information
          // so that a halo copy can be made on the receiving process
          get_required_nodal_information_helper(iproc,
                                                nod_pt,
                                                mesh_pt,
                                                n_cont_inter_values,
                                                send_unsigneds,
                                                send_doubles);
        }
        else // It was already added
        {
          // Indicate that this node doesn not need to be constructed on
          // the other process
          send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          std::stringstream junk;
          junk << "Node was already added [size=" << send_unsigneds.size()
               << "]; last entry: "
               << send_unsigneds[send_unsigneds.size() - 1];

          Flat_packed_unsigneds_string.push_back(junk.str());
#endif

          // This node is already an external haloed node, so tell
          // the other process its index in the equivalent external halo storage
          send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          Flat_packed_unsigneds_string.push_back(
            "haloed node found externally");
#endif
          send_unsigneds.push_back(external_haloed_node_index);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          Flat_packed_unsigneds_string.push_back("external haloed node index");
#endif
        }
      }
    }


    //==========start of add_external_haloed_master_node_helper===============
    /// Helper function to add external haloed node that is a master
    //========================================================================
    void add_external_haloed_master_node_helper(
      int& iproc,
      Node* master_nod_pt,
      Mesh* const& mesh_pt,
      int& n_cont_inter_values,
      Vector<unsigned>& send_unsigneds,
      Vector<double>& send_doubles)
    {
      // Check to see if this haloed node already exists in (internal)
      // haloed node storage with any processor
      bool found_internally = false;
      unsigned shared_node_index = 0;

      // Get vector of all shared nodes with processor iproc
      Vector<Node*> shared_node_pt;
      mesh_pt->get_shared_node_pt(iproc, shared_node_pt);

      // Search the internal haloed storage for this node
      Vector<Node*>::iterator it =
        std::find(shared_node_pt.begin(), shared_node_pt.end(), master_nod_pt);

      // Check if the node was found in shared storage
      if (it != shared_node_pt.end())
      {
        // Node found in (internal) haloed storage
        found_internally = true;
        // Store the index in this storage
        shared_node_index = it - shared_node_pt.begin();
      }

      /// / Slow search version without additional access function in Mesh class
      /// /Search the internal shared node storage for this node
      // for(unsigned i=0; i<mesh_pt->nshared_node(iproc); i++)
      // {
      //  if(master_nod_pt == mesh_pt->shared_node_pt(iproc,i))
      //   {
      //    //Node found in (internal) shared storage
      //    found_internally = true;
      //    shared_node_index = i;
      //    break;
      //   }
      // }

      // If we've found the node internally
      if (found_internally)
      {
        // Indicate that this node doesn not need to be constructed on
        // the other process
        send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        std::stringstream junk;
        junk << "Node was already added [size=" << send_unsigneds.size()
             << "]; last entry: " << send_unsigneds[send_unsigneds.size() - 1];

        Flat_packed_unsigneds_string.push_back(junk.str());
#endif

        // This node is already shared with processor iproc, so tell the other
        // processor its index in the shared node storage
        send_unsigneds.push_back(1);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("haloed node found internally");
#endif
        send_unsigneds.push_back(shared_node_index);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("(internal) haloed node index");
#endif
      }
      else
      {
        // Attempt to add node as an external haloed node
        unsigned n_ext_haloed_nod = mesh_pt->nexternal_haloed_node(iproc);
        unsigned external_haloed_node_index;
        external_haloed_node_index =
          mesh_pt->add_external_haloed_node_pt(iproc, master_nod_pt);

        // If it was added the returned index is the same as current storage
        // size
        if (external_haloed_node_index == n_ext_haloed_nod)
        {
          // Indicate that this node needs to be constructed on
          // the other process
          send_unsigneds.push_back(1);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          Flat_packed_unsigneds_string.push_back(
            "Node needs to be constructed[2]");
#endif

          // This gets all the required information for the specified
          // master node and stores it into MPI-sendable information
          // so that a halo copy can be made on the receiving process
          get_required_master_nodal_information_helper(iproc,
                                                       master_nod_pt,
                                                       mesh_pt,
                                                       n_cont_inter_values,
                                                       send_unsigneds,
                                                       send_doubles);
        }
        else // It was already added
        {
          // Indicate that this node doesn not need to be constructed on
          // the other process
          send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          Flat_packed_unsigneds_string.push_back("Node was already added[2]");
#endif

          // This node is already an external haloed node, so tell
          // the other process its index in the equivalent external halo storage
          send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          Flat_packed_unsigneds_string.push_back(
            "haloed node found externally");
#endif
          send_unsigneds.push_back(external_haloed_node_index);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          Flat_packed_unsigneds_string.push_back(
            "external haloed node index[2]");
#endif
        }
      }
    }


    //========start of get_required_nodal_information_helper==================
    /// Helper function to get the required nodal information from an
    /// external haloed node so that a fully-functional external halo
    /// node (and therefore element) can be created on the receiving process
    //========================================================================
    void get_required_nodal_information_helper(int& iproc,
                                               Node* nod_pt,
                                               Mesh* const& mesh_pt,
                                               int& n_cont_inter_values,
                                               Vector<unsigned>& send_unsigneds,
                                               Vector<double>& send_doubles)
    {
      // Tell the halo copy of this node how many values there are
      // [NB this may be different for nodes within the same element, e.g.
      //  when using Lagrange multipliers]
      unsigned n_val = nod_pt->nvalue();
      send_unsigneds.push_back(n_val);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
      Flat_packed_unsigneds_string.push_back("Number of values");
#endif

      unsigned n_dim = nod_pt->ndim();
      TimeStepper* time_stepper_pt = nod_pt->time_stepper_pt();

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
        send_unsigneds.push_back(1);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Node is on boundary");
#endif

        // Loop over the boundaries of the external mesh
        Vector<unsigned> boundaries;
        unsigned n_bnd = mesh_pt->nboundary();
        for (unsigned i_bnd = 0; i_bnd < n_bnd; i_bnd++)
        {
          // Which boundaries (could be more than one) is it on?
          if (nod_pt->is_on_boundary(i_bnd))
          {
            boundaries.push_back(i_bnd);
          }
        }
        unsigned nb = boundaries.size();
        send_unsigneds.push_back(nb);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        std::stringstream junk;
        junk << "Node is on " << nb << " boundaries";
        Flat_packed_unsigneds_string.push_back(junk.str());
#endif
        for (unsigned i = 0; i < nb; i++)
        {
          send_unsigneds.push_back(boundaries[i]);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
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
          send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          std::stringstream junk;
          Flat_packed_unsigneds_string.push_back(
            "No additional values were created by face element");
#endif
        }
        // Created additional values
        else
        {
          // How many?
          send_unsigneds.push_back(map_pt->size());
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          std::stringstream junk;
          junk << "Map size " << map_pt->size() << n_bnd;
          Flat_packed_unsigneds_string.push_back(junk.str());
#endif
          // Loop over entries in map and add to send data
          for (std::map<unsigned, unsigned>::iterator p = map_pt->begin();
               p != map_pt->end();
               p++)
          {
            send_unsigneds.push_back((*p).first);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
            std::stringstream junk;
            Flat_packed_unsigneds_string.push_back("Key of map entry");
#endif
            send_unsigneds.push_back((*p).second);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
            Flat_packed_unsigneds_string.push_back("Value of map entry");
#endif
          }
        }
      }
      else
      {
        // Not on any boundary
        send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
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
        AlgebraicMesh* alg_mesh_pt = dynamic_cast<AlgebraicMesh*>(mesh_pt);

        // Get default node update function ID
        unsigned update_id = alg_nod_pt->node_update_fct_id();
        send_unsigneds.push_back(update_id);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Alg Node update id");
#endif

        // Get reference values at default...
        unsigned n_ref_val = alg_nod_pt->nref_value();
        send_unsigneds.push_back(n_ref_val);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Alg Node n ref values");
#endif
        for (unsigned i_ref_val = 0; i_ref_val < n_ref_val; i_ref_val++)
        {
          send_doubles.push_back(alg_nod_pt->ref_value(i_ref_val));
        }

        // Access geometric objects at default...
        unsigned n_geom_obj = alg_nod_pt->ngeom_object();
        send_unsigneds.push_back(n_geom_obj);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
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
          send_unsigneds.push_back(found_geom_object);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
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
            send_doubles.push_back(
              solid_nod_pt->variable_position_pt()->value(t, i_val));
          }
        }
      }

      // Finally copy info required for all node types
      for (unsigned i_val = 0; i_val < n_val; i_val++)
      {
        for (unsigned t = 0; t < n_prev; t++)
        {
          send_doubles.push_back(nod_pt->value(t, i_val));
        }
      }

      // Now do positions
      for (unsigned idim = 0; idim < n_dim; idim++)
      {
        for (unsigned t = 0; t < n_prev; t++)
        {
          send_doubles.push_back(nod_pt->x(t, idim));
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
      Mesh* const& mesh_pt,
      int& n_cont_inter_values,
      Vector<unsigned>& send_unsigneds,
      Vector<double>& send_doubles)
    {
      // Need to send over dimension, position type and number of values
      send_unsigneds.push_back(master_nod_pt->ndim());
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
      Flat_packed_unsigneds_string.push_back("Master node ndim");
#endif
      send_unsigneds.push_back(master_nod_pt->nposition_type());
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
      Flat_packed_unsigneds_string.push_back("Master node npos_type");
#endif
      send_unsigneds.push_back(master_nod_pt->nvalue());
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
      Flat_packed_unsigneds_string.push_back("Master node nvalue");
#endif
      if (master_nod_pt->is_halo())
      {
        send_unsigneds.push_back(master_nod_pt->non_halo_proc_ID());
      }
      else
      {
        send_unsigneds.push_back(mesh_pt->communicator_pt()->my_rank());
      }
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
      Flat_packed_unsigneds_string.push_back(
        "Master node non-halo processor ID");
#endif

      // If it's a solid node, also need to send lagrangian dim and type
      SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(master_nod_pt);
      if (solid_nod_pt != 0)
      {
        send_unsigneds.push_back(solid_nod_pt->nlagrangian());
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Master solid node nlagr");
#endif
        send_unsigneds.push_back(solid_nod_pt->nlagrangian_type());
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Master solid node nlagr_type");
#endif
      }

      unsigned n_dim = master_nod_pt->ndim();
      TimeStepper* time_stepper_pt = master_nod_pt->time_stepper_pt();

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
        send_unsigneds.push_back(1);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Master node is on boundary");
#endif
        // Loop over the boundaries of the external mesh
        Vector<unsigned> boundaries;
        unsigned n_bnd = mesh_pt->nboundary();
        for (unsigned i_bnd = 0; i_bnd < n_bnd; i_bnd++)
        {
          // Which boundaries (could be more than one) is it on?
          if (master_nod_pt->is_on_boundary(i_bnd))
          {
            boundaries.push_back(i_bnd);
          }
        }
        unsigned nb = boundaries.size();
        send_unsigneds.push_back(nb);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        std::stringstream junk;
        junk << "Master node is on " << nb << " boundaries";
        Flat_packed_unsigneds_string.push_back(junk.str());
#endif
        for (unsigned i = 0; i < nb; i++)
        {
          send_unsigneds.push_back(boundaries[i]);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
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
          send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
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
          send_unsigneds.push_back(map_pt->size());
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          std::stringstream junk;
          junk << "Map size for master node " << map_pt->size() << n_bnd;
          Flat_packed_unsigneds_string.push_back(junk.str());
#endif
          // Loop over entries in map and add to send data
          for (std::map<unsigned, unsigned>::iterator p = map_pt->begin();
               p != map_pt->end();
               p++)
          {
            send_unsigneds.push_back((*p).first);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
            std::stringstream junk;
            Flat_packed_unsigneds_string.push_back(
              "Key of map entry for master node");
#endif
            send_unsigneds.push_back((*p).second);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
            Flat_packed_unsigneds_string.push_back(
              "Value of map entry for master node");
#endif
          }
        }
      }
      else
      {
        // Not on any boundary
        send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
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
        AlgebraicMesh* alg_mesh_pt = dynamic_cast<AlgebraicMesh*>(mesh_pt);

        // Get default node update function ID
        unsigned update_id = alg_nod_pt->node_update_fct_id();
        send_unsigneds.push_back(update_id);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Master Alg Node update id");
#endif

        // Get reference values at default...
        unsigned n_ref_val = alg_nod_pt->nref_value();
        send_unsigneds.push_back(n_ref_val);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        Flat_packed_unsigneds_string.push_back("Master Alg Node n ref values");
#endif
        for (unsigned i_ref_val = 0; i_ref_val < n_ref_val; i_ref_val++)
        {
          send_doubles.push_back(alg_nod_pt->ref_value(i_ref_val));
        }

        // Access geometric objects at default...
        unsigned n_geom_obj = alg_nod_pt->ngeom_object();
        send_unsigneds.push_back(n_geom_obj);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
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
          send_unsigneds.push_back(found_geom_object);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
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
        oomph_info << "Adding external haloed master node: " << master_nod_pt
                   << " at " << master_nod_pt->x(0) << ", "
                   << master_nod_pt->x(1) << " ]" << std::endl;
        // Loop over current external haloed elements - has the element which
        // controls the node update for this node been added yet?
        GeneralisedElement* macro_node_update_el_pt =
          macro_nod_pt->node_update_element_pt();

        // BENFLAG: Check that the node's macro element node update element
        // actually contains the node
        oomph_info << "Master node's macro update element:" << std::endl;
        bool really_bad = true;
        FiniteElement* mac_el_pt =
          dynamic_cast<FiniteElement*>(macro_node_update_el_pt);
        for (unsigned j = 0; j < mac_el_pt->nnode(); j++)
        {
          oomph_info << mac_el_pt->node_pt(j) << ": [ "
                     << mac_el_pt->node_pt(j)->x(0) << ", "
                     << mac_el_pt->node_pt(j)->x(1) << " ] " << std::endl;
          if (mac_el_pt->node_pt(j) == master_nod_pt)
          {
            really_bad = false;
            // oomph_info << "Found it!" << std::endl;
          }
        }
        if (really_bad == true)
        {
          oomph_info << "This is REALLY BAD! The master node is not part of "
                        "its own update element..."
                     << std::endl;
        }

        // BENFLAG: Search internal storage for node update element
        Vector<GeneralisedElement*> int_haloed_el_pt(
          mesh_pt->haloed_element_pt(iproc));
        // unsigned n_int_haloed_el=int_haloed_el_pt.size();
        Vector<GeneralisedElement*>::iterator it =
          std::find(int_haloed_el_pt.begin(),
                    int_haloed_el_pt.end(),
                    macro_node_update_el_pt);
        if (it != int_haloed_el_pt.end())
        {
          // Found in internal haloed storage
          unsigned int_haloed_el_index = it - int_haloed_el_pt.begin();
          oomph_info << "Found internally at index " << int_haloed_el_index
                     << std::endl;
          // BENFLAG: Check index corresponds to correct element
          if ((mesh_pt->haloed_element_pt(iproc))[int_haloed_el_index] !=
              macro_node_update_el_pt)
          {
            oomph_info << "Found wrong index!!!" << std::endl;
            throw;
          }
          else
          {
            oomph_info << "index and proc are correct in internal storage."
                       << std::endl;
            oomph_info << "i.e. "
                       << (mesh_pt->haloed_element_pt(
                            iproc))[int_haloed_el_index]
                       << "==" << macro_node_update_el_pt << std::endl;
          }


          // Say haloed element already exists
          send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          Flat_packed_unsigneds_string.push_back(
            "External haloed element already exists");
#endif
          // Say found internally
          send_unsigneds.push_back(1);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          Flat_packed_unsigneds_string.push_back(
            "Haloed element found internally");
#endif
          // Say what the index is
          send_unsigneds.push_back(int_haloed_el_index);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          Flat_packed_unsigneds_string.push_back(
            "Index of existing internal haloed element");
#endif
          // BENFLAG:
          FiniteElement* tmp_el_pt = dynamic_cast<FiniteElement*>(
            (mesh_pt->haloed_element_pt(iproc))[int_haloed_el_index]);
          oomph_info << "Internal haloed element (" << tmp_el_pt
                     << ") already exists..." << std::endl;
          oomph_info << "on proc " << iproc << " at index "
                     << int_haloed_el_index << std::endl;
          for (unsigned j = 0; j < tmp_el_pt->nnode(); j++)
          {
            oomph_info << tmp_el_pt->node_pt(j) << "=="
                       << dynamic_cast<FiniteElement*>(
                            (mesh_pt->haloed_element_pt(
                              iproc))[int_haloed_el_index])
                            ->node_pt(j)
                       << " at [ " << tmp_el_pt->node_pt(j)->x(0) << ", "
                       << tmp_el_pt->node_pt(j)->x(1) << " ]" << std::endl;
          }

          // oomph_info << "now " << tmp_el_pt << "==" <<
          // (mesh_pt->haloed_element_pt(iproc))[int_haloed_el_index] << "==" <<
          // macro_node_update_el_pt << std::endl;

          //       //BENFLAG: Add to external storage too
          //       unsigned n_ext_haloed_el=mesh_pt->
          //        nexternal_haloed_element(iproc);
          //       unsigned external_haloed_el_index;
          //       external_haloed_el_index=mesh_pt->
          //        add_external_haloed_element_pt(iproc,macro_node_update_el_pt);
          //
          //       // If it was already added, say
          //       if (external_haloed_el_index!=n_ext_haloed_el)
          //        {
          //         oomph_info << "Element (" << tmp_el_pt << "==" <<
          //         macro_node_update_el_pt << "==" <<
          //         mesh_pt->external_haloed_element_pt(iproc,external_haloed_el_index)
          //         << ") also exists in external storage with proc " << iproc
          //         << " at index " << external_haloed_el_index << " of " <<
          //         n_ext_haloed_el << "..." << std::endl;
          //         // Say also exists in external storage
          //         send_unsigneds.push_back(1234);
          //#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          //         Flat_packed_unsigneds_string.push_back("Haloed element also
          //         found externally");
          //#endif
          //         // Say what the index is
          //         send_unsigneds.push_back(external_haloed_el_index);
          //#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          //         Flat_packed_unsigneds_string.push_back("Index of existing
          //         internal haloed element");
          //#endif
          //         oomph_info << "sent external_haloed_el_index = " <<
          //         external_haloed_el_index << std::endl;
          //
          //        }
          //       else
          //        {
          //         oomph_info << "Element didn't exist in external storage
          //         with proc " << iproc << "..." << std::endl;
          //         // Say doesn't exists in external storage
          //         send_unsigneds.push_back(0);
          //#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          //         Flat_packed_unsigneds_string.push_back("Haloed element not
          //         also found externally");
          //#endif
          //        }
        }
        else
        {
          unsigned n_ext_haloed_el = mesh_pt->nexternal_haloed_element(iproc);
          unsigned external_haloed_el_index;
          external_haloed_el_index = mesh_pt->add_external_haloed_element_pt(
            iproc, macro_node_update_el_pt);

          // If it wasn't already added, we need to create a halo copy
          if (external_haloed_el_index == n_ext_haloed_el)
          {
            send_unsigneds.push_back(1);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
            Flat_packed_unsigneds_string.push_back(
              "Master Node needs to be constructed");
#endif
            // Cast to a finite elemnet
            FiniteElement* macro_node_update_finite_el_pt =
              dynamic_cast<FiniteElement*>(macro_node_update_el_pt);

            // We're using macro elements to update...
            MacroElementNodeUpdateMesh* macro_mesh_pt =
              dynamic_cast<MacroElementNodeUpdateMesh*>(mesh_pt);
            if (macro_mesh_pt != 0)
            {
              send_unsigneds.push_back(1);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
              Flat_packed_unsigneds_string.push_back(
                "Mesh is macro element mesh");
#endif
              // Need to send the macro element number in the mesh across
              MacroElement* macro_el_pt =
                macro_node_update_finite_el_pt->macro_elem_pt();
              unsigned macro_el_num = macro_el_pt->macro_element_number();
              send_unsigneds.push_back(macro_el_num);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
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
                  send_doubles.push_back(q_el_pt->s_macro_ll(i_dim));
                  send_doubles.push_back(q_el_pt->s_macro_ur(i_dim));
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
              send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
              Flat_packed_unsigneds_string.push_back(
                "Mesh is not a macro element mesh");
#endif
            }

            // If the element is p-refineable we need to send the p-order so
            // that the halo version can be constructed correctly
            PRefineableElement* p_refineable_el_pt =
              dynamic_cast<PRefineableElement*>(macro_node_update_finite_el_pt);
            if (p_refineable_el_pt != 0)
            {
              send_unsigneds.push_back(1);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
              Flat_packed_unsigneds_string.push_back("Element is p-refineable");
#endif
              // Send p-order of macro element node update element
              unsigned p_order = p_refineable_el_pt->p_order();
              send_unsigneds.push_back(p_order);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
              Flat_packed_unsigneds_string.push_back("p-order of element");
#endif
            }
            else
            {
              send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
              Flat_packed_unsigneds_string.push_back(
                "Element is not p-refineable");
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
                                                  mesh_pt,
                                                  n_cont_inter_values,
                                                  send_unsigneds,
                                                  send_doubles);
            }
          }
          else // The external haloed element already exists
          {
            // Say haloed element already exists
            send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
            Flat_packed_unsigneds_string.push_back(
              "External haloed element already exists");
#endif
            // Say found externally
            send_unsigneds.push_back(0);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
            Flat_packed_unsigneds_string.push_back(
              "Haloed element found externally");
#endif
            send_unsigneds.push_back(external_haloed_el_index);
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
            Flat_packed_unsigneds_string.push_back(
              "Index of existing external haloed element");
#endif
            // BENFLAG:
            oomph_info << "External haloed element already exists..."
                       << std::endl;
            FiniteElement* tmp_el_pt =
              dynamic_cast<FiniteElement*>(mesh_pt->external_haloed_element_pt(
                iproc, external_haloed_el_index));
            oomph_info << "on proc " << iproc << " at index "
                       << external_haloed_el_index << std::endl;
            for (unsigned j = 0; j < tmp_el_pt->nnode(); j++)
            {
              oomph_info << tmp_el_pt->node_pt(j) << " at [ "
                         << tmp_el_pt->node_pt(j)->x(0) << ", "
                         << tmp_el_pt->node_pt(j)->x(1) << " ]" << std::endl;
            }
          }
        } // End of case where not found internally

      } // end of MacroElementNodeUpdateNode check

      // Is it a SolidNode?
      if (solid_nod_pt != 0)
      {
        unsigned n_val = solid_nod_pt->variable_position_pt()->nvalue();
        for (unsigned i_val = 0; i_val < n_val; i_val++)
        {
          for (unsigned t = 0; t < n_prev; t++)
          {
            send_doubles.push_back(
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
          send_doubles.push_back(master_nod_pt->value(t, i_val));
        }
      }

      // Now do positions
      for (unsigned idim = 0; idim < n_dim; idim++)
      {
        for (unsigned t = 0; t < n_prev; t++)
        {
          send_doubles.push_back(master_nod_pt->x(t, idim));
        }
      }
    }


    //=======start of add_external_halo_node_helper===========================
    /// Helper functiono to add external halo node that is not a master
    //========================================================================
    void add_external_halo_node_helper(Node*& new_nod_pt,
                                       Mesh* const& mesh_pt,
                                       unsigned& loc_p,
                                       unsigned& node_index,
                                       FiniteElement* const& new_el_pt,
                                       int& n_cont_inter_values,
                                       unsigned& counter_for_recv_unsigneds,
                                       Vector<unsigned>& recv_unsigneds,
                                       unsigned& counter_for_recv_doubles,
                                       Vector<double>& recv_doubles)
    {
      // Given the node and the external mesh, and received information
      // about them from process loc_p, construct them on the current process
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
      oomph_info << "Rec:" << counter_for_recv_unsigneds
                 << " Bool: New node needs to be constructed "
                 << recv_unsigneds[counter_for_recv_unsigneds] << std::endl;
#endif
      if (recv_unsigneds[counter_for_recv_unsigneds++] == 1)
      {
        // Construct a new node based upon sent information
        construct_new_external_halo_node_helper(new_nod_pt,
                                                loc_p,
                                                node_index,
                                                new_el_pt,
                                                mesh_pt,
                                                counter_for_recv_unsigneds,
                                                recv_unsigneds,
                                                counter_for_recv_doubles,
                                                recv_doubles);
      }
      else
      {
        // Need to check which storage we should copy this halo node from
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        oomph_info << "Rec:" << counter_for_recv_unsigneds
                   << "  Existing external halo node was found externally (0) "
                      "or internally (1): "
                   << recv_unsigneds[counter_for_recv_unsigneds] << std::endl;
#endif
        unsigned node_found_internally =
          recv_unsigneds[counter_for_recv_unsigneds++];
        if (node_found_internally)
        {
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          oomph_info << "Rec:" << counter_for_recv_unsigneds
                     << "  index of existing (internal) halo master node "
                     << recv_unsigneds[counter_for_recv_unsigneds] << std::endl;
#endif

          // Copy node from received location
          new_nod_pt = mesh_pt->shared_node_pt(
            loc_p, recv_unsigneds[counter_for_recv_unsigneds++]);

          new_el_pt->node_pt(node_index) = new_nod_pt;
        }
        else
        {
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          oomph_info << "Rec:" << counter_for_recv_unsigneds
                     << "  Index of existing external halo node "
                     << recv_unsigneds[counter_for_recv_unsigneds] << std::endl;
#endif

          // Copy node from received location
          new_nod_pt = mesh_pt->external_halo_node_pt(
            loc_p, recv_unsigneds[counter_for_recv_unsigneds++]);

          new_el_pt->node_pt(node_index) = new_nod_pt;
        }
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
      Mesh* const& mesh_pt,
      unsigned& counter_for_recv_unsigneds,
      Vector<unsigned>& recv_unsigneds,
      unsigned& counter_for_recv_doubles,
      Vector<double>& recv_doubles)
    {
      // The first entry indicates the number of values at this new Node
      // (which may be different across the same element e.g. Lagrange
      // multipliers)
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
      oomph_info << "Rec:" << counter_for_recv_unsigneds
                 << "  Number of values of external halo node "
                 << recv_unsigneds[counter_for_recv_unsigneds] << std::endl;
#endif
      unsigned n_val = recv_unsigneds[counter_for_recv_unsigneds++];

      // Null TimeStepper for now
      TimeStepper* time_stepper_pt = 0;
      // Default number of previous values to 1
      unsigned n_prev = 1;

      // Just take timestepper from a node
      // Let's use first node of first element since this must exist
      time_stepper_pt =
        mesh_pt->finite_element_pt(0)->node_pt(0)->time_stepper_pt();

      // If this node was on a boundary then it needs to
      // be on the same boundary here
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
      oomph_info << "Rec:" << counter_for_recv_unsigneds
                 << "  Is node on boundary? "
                 << recv_unsigneds[counter_for_recv_unsigneds] << std::endl;
#endif
      if (recv_unsigneds[counter_for_recv_unsigneds++] == 1)
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
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        oomph_info << "Rec:" << counter_for_recv_unsigneds
                   << " Number of boundaries the node is on: "
                   << recv_unsigneds[counter_for_recv_unsigneds] << std::endl;
#endif
        unsigned nb = recv_unsigneds[counter_for_recv_unsigneds++];
        for (unsigned i = 0; i < nb; i++)
        {
          // Boundary number
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          oomph_info << "Rec:" << counter_for_recv_unsigneds
                     << "  Node is on boundary "
                     << recv_unsigneds[counter_for_recv_unsigneds] << std::endl;
#endif
          unsigned i_bnd = recv_unsigneds[counter_for_recv_unsigneds++];
          mesh_pt->add_boundary_node(i_bnd, new_nod_pt);
        }

        // Do we have additional values created by face elements?
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        oomph_info << "Rec:" << counter_for_recv_unsigneds
                   << " Number of additional values created by face element "
                   << recv_unsigneds[counter_for_recv_unsigneds] << std::endl;
#endif
        unsigned n_entry = recv_unsigneds[counter_for_recv_unsigneds++];
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

#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
            oomph_info << "Rec:" << counter_for_recv_unsigneds
                       << " Key of map entry"
                       << recv_unsigneds[counter_for_recv_unsigneds]
                       << std::endl;
#endif
            unsigned first = recv_unsigneds[counter_for_recv_unsigneds++];

#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
            oomph_info << "Rec:" << counter_for_recv_unsigneds
                       << " Value of map entry"
                       << recv_unsigneds[counter_for_recv_unsigneds]
                       << std::endl;
#endif
            unsigned second = recv_unsigneds[counter_for_recv_unsigneds++];

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
      mesh_pt->add_external_halo_node_pt(loc_p, new_nod_pt);

      // Is the new constructed node Algebraic?
      AlgebraicNode* new_alg_nod_pt = dynamic_cast<AlgebraicNode*>(new_nod_pt);

      // If it is algebraic, its node update functions will
      // not yet have been set up properly
      if (new_alg_nod_pt != 0)
      {
        // The AlgebraicMesh is the external mesh
        AlgebraicMesh* alg_mesh_pt = dynamic_cast<AlgebraicMesh*>(mesh_pt);

        /// The first entry of All_alg_nodal_info contains
        /// the default node update id
        /// e.g. for the quarter circle there are
        /// "Upper_left_box", "Lower right box" etc...
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        oomph_info << "Rec:" << counter_for_recv_unsigneds
                   << "  Alg node update id "
                   << recv_unsigneds[counter_for_recv_unsigneds] << std::endl;
#endif

        unsigned update_id = recv_unsigneds[counter_for_recv_unsigneds++];

        Vector<double> ref_value;

        // The size of this vector is in the next entry
        // of All_alg_nodal_info
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        oomph_info << "Rec:" << counter_for_recv_unsigneds
                   << "  Alg node # of ref values "
                   << recv_unsigneds[counter_for_recv_unsigneds] << std::endl;
#endif
        unsigned n_ref_val = recv_unsigneds[counter_for_recv_unsigneds++];

        // The reference values themselves are in
        // All_alg_ref_value
        ref_value.resize(n_ref_val);
        for (unsigned i_ref = 0; i_ref < n_ref_val; i_ref++)
        {
          ref_value[i_ref] = recv_doubles[counter_for_recv_doubles++];
        }

        Vector<GeomObject*> geom_object_pt;
        /// again we need the size of this vector as it varies
        /// between meshes; we also need some indication
        /// as to which geometric object should be used...

        // The size of this vector is in the next entry
        // of All_alg_nodal_info
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
        oomph_info << "Rec:" << counter_for_recv_unsigneds
                   << "  Alg node # of geom objects "
                   << recv_unsigneds[counter_for_recv_unsigneds] << std::endl;
#endif
        unsigned n_geom_obj = recv_unsigneds[counter_for_recv_unsigneds++];

        // The remaining indices are in the rest of
        // All_alg_nodal_info
        geom_object_pt.resize(n_geom_obj);
        for (unsigned i_geom = 0; i_geom < n_geom_obj; i_geom++)
        {
#ifdef ANNOTATE_MISSING_MASTERS_COMMUNICATION
          oomph_info << "Rec:" << counter_for_recv_unsigneds
                     << "  Alg node: geom object index "
                     << recv_unsigneds[counter_for_recv_unsigneds] << std::endl;
#endif
          unsigned geom_index = recv_unsigneds[counter_for_recv_unsigneds++];
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
          dynamic_cast<MacroElementNodeUpdateMesh*>(mesh_pt);
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
              t, i_val, recv_doubles[counter_for_recv_doubles++]);
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
            t, i_val, recv_doubles[counter_for_recv_doubles++]);
        }
      }

      // Get copied history values for positions
      unsigned n_dim = new_nod_pt->ndim();
      for (unsigned idim = 0; idim < n_dim; idim++)
      {
        for (unsigned t = 0; t < n_prev; t++)
        {
          // Copy to coordinate
          new_nod_pt->x(t, idim) = recv_doubles[counter_for_recv_doubles++];
        }
      }
    }


#endif

  } // namespace Missing_masters_functions

} // namespace oomph
