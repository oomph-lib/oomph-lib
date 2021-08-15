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
// Non-inline and non-templated functions for RefineableQElement<1> class

#include <algorithm>

// oomph-lib headers
#include "mesh.h"
#include "algebraic_elements.h"
#include "macro_element_node_update_element.h"
#include "refineable_line_element.h"

namespace oomph
{
  //==========================================================================
  /// Setup static matrix for coincidence between son nodal points and father
  /// boundaries:
  ///
  /// Father_bound[nnode_1d](nnode_son,son_type) = {L/R/OMEGA}
  ///
  /// so that node nnode_son in element of type son_type lies on boundary/
  /// vertex Father_bound[nnode_1d](nnode_son,son_type) in its father
  /// element. If the node doesn't lie on a boundary the value is OMEGA.
  //==========================================================================
  void RefineableQElement<1>::setup_father_bounds()
  {
    using namespace BinaryTreeNames;

    // Find the number of nodes along a 1D edge (which is the number of nodes
    // in the element for a 1D element!)
    const unsigned n_node = nnode_1d();

    // Allocate space for the boundary information
    Father_bound[n_node].resize(n_node, 2);

    // Initialise: By default points are not on the boundary
    for (unsigned n = 0; n < n_node; n++)
    {
      for (unsigned ison = 0; ison < 2; ison++)
      {
        Father_bound[n_node](n, ison) = Tree::OMEGA;
      }
    }

    // Left-hand son:
    // --------------

    // L node (0) is the L node of the parent
    Father_bound[n_node](0, L) = L;

    // Other boundary is in the interior

    // Right-hand son:
    // ---------------

    // R node (n_node-1) is the R node of the parent
    Father_bound[n_node](n_node - 1, R) = R;

    // Other boundary is in the interior
  }

  //==========================================================================
  /// If a neighbouring element has already created a node at a position
  /// corresponding to the local fractional position within the present
  /// element, s_fraction, return a pointer to that node. If not, return
  /// NULL (0). If the node is periodic the flag is_periodic will be true.
  //==========================================================================
  Node* RefineableQElement<1>::node_created_by_neighbour(
    const Vector<double>& s_fraction, bool& is_periodic)
  {
    using namespace BinaryTreeNames;

    // Initialise edge of current element on which node lies
    int edge_in_current = OMEGA;

    // Determine the edge of the current element on which the node lies
    if (s_fraction[0] == 0.0)
    {
      edge_in_current = L;
    }
    if (s_fraction[0] == 1.0)
    {
      edge_in_current = R;
    }

    // If the node does not lie on an edge then there is no neighbour:
    // return NULL
    if (edge_in_current == OMEGA)
    {
      return 0;
    }

    // Allocate storage for edge in neighbouring element
    int edge_in_neighbour;

    // Allocate storage for difference in size between current and
    // neighbouring element
    int diff_level;

    // Allocate storage for local coordinate of node in neighbouring tree
    Vector<double> s_in_neighbour(1);

    // Allocate storage for flag indicating if the node is not in the same
    // binary tree
    bool in_neighbouring_tree;

    // Allocate storage for the pointer to the neighbouring element
    // (using its binary tree representation)
    BinaryTree* neighbour_pt;

    // Find pointer to neighbouring element along the edge in question and
    // calculate the local coordinate of the node within that element
    // s_in_neighbour
    neighbour_pt = binary_tree_pt()->gteq_edge_neighbour(edge_in_current,
                                                         s_in_neighbour,
                                                         edge_in_neighbour,
                                                         diff_level,
                                                         in_neighbouring_tree);

    // If a neighbour exists...
    if (neighbour_pt != 0)
    {
      // ...check whether its nodes have been created yet
      if (neighbour_pt->object_pt()->nodes_built())
      {
        // If they have, find the node in question in the neighbour
        Node* neighbour_node_pt =
          neighbour_pt->object_pt()->get_node_at_local_coordinate(
            s_in_neighbour);

        // If there is no node at this position, there is a problem, since in
        // a 1D element (whose nodes have already been built) there should
        // ALWAYS be a node at each edge of the element.
        if (neighbour_node_pt == 0)
        {
          std::string error_message =
            "Problem: an element claims to have had its nodes built, yet\n";
          error_message += "it is missing (a least) a node at its edge.\n";
          throw OomphLibError(
            error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
        // Otherwise, carry on
        else
        {
          // Now work out whether it's a periodic boundary (this is only
          // (possible if we have moved into a neighbouring tree)
          if (in_neighbouring_tree)
          {
            // Return whether the neighbour is periodic
            is_periodic = binary_tree_pt()->root_pt()->is_neighbour_periodic(
              edge_in_current);
          }
          // Return the pointer to the neighbouring node
          return neighbour_node_pt;
        }
      }
    }
    // Node not found, return null
    return 0;
  }

  //==========================================================================
  /// Build the element by doing the following:
  /// - Give it nodal positions (by establishing the pointers to its nodes).
  ///   In the process create new nodes where required (i.e. if they don't
  ///   already exist in the father element and have not already been created
  ///   whilst building new neighbouring elements). Node building involves
  ///   the following steps:
  ///   - Get the nodal position from the father element.
  ///   - Establish the time-history of the newly created nodal point
  ///     (its coordinates and previous values) consistent with the father's
  ///     history.
  ///   - Add the new node to the mesh itself.
  ///   - Doc newly created nodes in "new_nodes.dat" stored in the directory
  ///     of the DocInfo object (only if it's open!).
  ///   - NOTE: Unlike in higher-dimensional elements, in 1D it is impossible
  ///     for newly-created nodes to be on a mesh boundary, since any boundary
  ///     nodes must exist in the initial (coarse) mesh. Therefore it is not
  ///     necessary to add any nodes to the mesh's boundary node storage
  ///     schemes, and we always create normal "bulk" nodes.
  /// - Once the element has a full complement of nodes, excute the element-
  ///   specific further_build() (empty by default -- must be overloaded for
  ///   specific elements). This deals with any build operations that are not
  ///   included in the generic process outlined above. For instance, in
  ///   Crouzeix Raviart elements we need to initialise the internal pressure
  ///   pressure values in a manner consistent with the pressure distribution
  ///   in the father element.
  //==========================================================================
  void RefineableQElement<1>::build(Mesh*& mesh_pt,
                                    Vector<Node*>& new_node_pt,
                                    bool& was_already_built,
                                    std::ofstream& new_nodes_file)
  {
    using namespace BinaryTreeNames;

    // Find the number of nodes along a 1D edge (which is the number of nodes
    // in the element for a 1D element!)
    const unsigned n_node = nnode_1d();

    // Check whether static father_bound needs to be created
    if (Father_bound[n_node].nrow() == 0)
    {
      setup_father_bounds();
    }

    // Pointer to the current element's father (in binary tree impersonation)
    BinaryTree* father_pt =
      dynamic_cast<BinaryTree*>(binary_tree_pt()->father_pt());

    // What type of son is the current element?
    // Ask its binary tree representation...
    const int son_type = Tree_pt->son_type();

    // Has somebody built the current element already?
    // Check this by determining whether or not the first node has been built

    // If this element has not already been built...
    if (!nodes_built())
    {
#ifdef PARANOID
      if (father_pt == 0)
      {
        std::string error_message =
          "Something fishy here: I have no father and yet \n";
        error_message += "I have no nodes. Who has created me then?!\n";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Indicate status
      was_already_built = false;

      // Return pointer to father element
      RefineableQElement<1>* father_el_pt =
        dynamic_cast<RefineableQElement<1>*>(father_pt->object_pt());

      // Timestepper should be the same for all nodes in father element.
      // Use it create timesteppers for new nodes.
      TimeStepper* time_stepper_pt =
        father_el_pt->node_pt(0)->time_stepper_pt();

      // Determine number of history values (including present)
      const unsigned ntstorage = time_stepper_pt->ntstorage();

      // Currently we can't handle the case of generalised coordinates
      // since we haven't established how they should be interpolated.
      // Buffer this case:
      if (father_el_pt->node_pt(0)->nposition_type() != 1)
      {
        throw OomphLibError("Can't handle generalised nodal positions (yet).",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Set up vertex coordinates in the father element:
      // ------------------------------------------------

      // Allocate storage for the vertex coordinates s_left and s_right of
      // the current element as viewed by this element's father (i.e. s_left[0]
      // stores the local coordinate within the father element at which the
      // node on the current element's left-hand edge is located. Likewise,
      // s_right[0] stores the local coordinate within the father element at
      // which the node on the current element's right-hand edge is located).
      Vector<double> s_left(1);
      Vector<double> s_right(1);

      // In order to set up the vertex coordinates we need to know which
      // type of son the current element is
      switch (son_type)
      {
        case L:
          s_left[0] = -1.0;
          s_right[0] = 0.0;
          break;

        case R:
          s_left[0] = 0.0;
          s_right[0] = 1.0;
          break;
      }

      // Pass macro element pointer on to sons and set coordinates in
      // macro element
      // hierher why can I see this?
      if (father_el_pt->Macro_elem_pt != 0)
      {
        set_macro_elem_pt(father_el_pt->Macro_elem_pt);

        s_macro_ll(0) =
          father_el_pt->s_macro_ll(0) +
          0.5 * (s_left[0] + 1.0) *
            (father_el_pt->s_macro_ur(0) - father_el_pt->s_macro_ll(0));
        s_macro_ur(0) =
          father_el_pt->s_macro_ll(0) +
          0.5 * (s_right[0] + 1.0) *
            (father_el_pt->s_macro_ur(0) - father_el_pt->s_macro_ll(0));
      }

      // If the father element hasn't been generated yet, we're stuck...
      if (father_el_pt->node_pt(0) == 0)
      {
        throw OomphLibError(
          "Trouble: father_el_pt->node_pt(0)==0\n Can't build son element!\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      // Otherwise, carry on
      else
      {
        // Allocate storage for the location of a node in terms of the local
        // coordinate of the father element
        Vector<double> s_in_father(1);

        // Allocate storage for the fractional position (in the current
        // element) of the node in the direction of s[0]
        Vector<double> s_fraction(1);

        // Loop over all nodes in the element
        for (unsigned n = 0; n < n_node; n++)
        {
          // Get the fractional position (in the current element) of the node
          // in the direction of s[0]
          s_fraction[0] = local_one_d_fraction_of_node(n, 0);

          // Evaluate the local coordinate of the node in the father element
          s_in_father[0] = s_left[0] + (s_right[0] - s_left[0]) * s_fraction[0];

          // Initialise flag: So far, this node hasn't been built or copied yet
          bool node_done = false;

          // Get the pointer to the node in the father (returns NULL if there
          // is no node at this position)
          Node* created_node_pt =
            father_el_pt->get_node_at_local_coordinate(s_in_father);

          // Does this node already exist in father element?
          // -----------------------------------------------

          // If it does...
          if (created_node_pt != 0)
          {
            // ...copy node across
            node_pt(n) = created_node_pt;

            // Make sure that we update the values at the node so that they
            // are consistent with the present representation. This is only
            // needed for mixed interpolation where the value at the father
            // could now become active.

            // Loop over all history values
            for (unsigned t = 0; t < ntstorage; t++)
            {
              // Get values from father element
              // Note: get_interpolated_values() sets Vector size itself
              Vector<double> prev_values;
              father_el_pt->get_interpolated_values(
                t, s_in_father, prev_values);

              // Find the minimum number of values (either those stored at the
              // node, or those returned by the function)
              unsigned n_val_at_node = created_node_pt->nvalue();
              unsigned n_val_from_function = prev_values.size();

              // Use the ternary conditional operator here
              unsigned n_var = n_val_at_node < n_val_from_function ?
                                 n_val_at_node :
                                 n_val_from_function;

              // Assign the values that we can
              for (unsigned k = 0; k < n_var; k++)
              {
                created_node_pt->set_value(t, k, prev_values[k]);
              }
            }

            // Indicate that node has been created by copy
            node_done = true;
          }

          // Node does not exist in father element but might already
          // -------------------------------------------------------
          // have been created by neighbouring elements
          // ------------------------------------------

          else
          {
            // Was the node created by one of its neighbours? Whether or not
            // the node lies on an edge can be calculated from the fractional
            // position.
            bool is_periodic = false;
            created_node_pt =
              node_created_by_neighbour(s_fraction, is_periodic);

            // If the node was so created...
            if (created_node_pt != 0)
            {
              // ...assign the pointer
              node_pt(n) = created_node_pt;

              // Indicate that node has been created by copy
              node_done = true;

              // In a 1D mesh there is no way that a periodic node (which must
              // be on a boundary) can exist without being part of the initial
              // (coarse) mesh. Therefore issue an error message if
              // node_created_by_neighbour(...) returns `is_periodic==true'.
#ifdef PARANOID
              if (is_periodic)
              {
                std::string error_message =
                  "node_created_by_neighbour returns a node which claims\n";
                error_message += "to be periodic. In a 1D mesh any periodic "
                                 "nodes must exist\n";
                error_message += "in the initial (coarse) mesh.";

                throw OomphLibError(error_message,
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif
            }
          }

          // Node does not exist in father element or neighbouring element
          // -------------------------------------------------------------

          // If the node has not been built anywhere ---> build it here
          if (!node_done)
          {
            // In a 1D mesh any node which lies on the boundary must exist in
            // the initial (coarse) mesh, so any newly-built nodes cannot be
            // boundary nodes. Therefore we always create a normal "bulk" node.

            // Create node and set the pointer to it from the element
            created_node_pt = construct_node(n, time_stepper_pt);

            // Add to vector of new nodes
            new_node_pt.push_back(created_node_pt);

            // Now we set the position and values at the newly created node.
            // In the first instance use macro element or FE representation
            // to create past and present nodal positions.
            // (THIS STEP SHOULD NOT BE SKIPPED FOR ALGEBRAIC ELEMENTS AS NOT
            // ALL OF THEM NECESSARILY IMPLEMENT NONTRIVIAL NODE UPDATE
            // FUNCTIONS. CALLING THE NODE UPDATE FOR SUCH ELEMENTS/NODES WILL
            // LEAVE THEIR NODAL POSITIONS WHERE THEY WERE (THIS IS APPROPRIATE
            // ONCE THEY HAVE BEEN GIVEN POSITIONS) BUT WILL NOT ASSIGN SENSIBLE
            // INITIAL POSITIONS!)

            // Loop over history values
            for (unsigned t = 0; t < ntstorage; t++)
            {
              // Allocate storage for the previous position of the node
              Vector<double> x_prev(1);

              // Get position from father element -- this uses the macro
              // element representation if appropriate.
              father_el_pt->get_x(t, s_in_father, x_prev);

              // Set the previous position of the new node
              created_node_pt->x(t, 0) = x_prev[0];

              // Allocate storage for the previous values at the node
              // NOTE: the size of this vector is equal to the number of values
              // (e.g. 3 velocity components and 1 pressure, say)
              // associated with each node and NOT the number of history values
              // which the node stores!
              Vector<double> prev_values;

              // Get values from father element
              // Note: get_interpolated_values() sets Vector size itself.
              father_el_pt->get_interpolated_values(
                t, s_in_father, prev_values);

              // Determine the number of values at the new node
              const unsigned n_value = created_node_pt->nvalue();

              // Loop over all values and set the previous values
              for (unsigned v = 0; v < n_value; v++)
              {
                created_node_pt->set_value(t, v, prev_values[v]);
              }
            } // End of loop over history values

            // Add new node to mesh
            mesh_pt->add_node_pt(created_node_pt);

          } // End of case where we build the node ourselves

          // Check if the element is an algebraic element
          AlgebraicElementBase* alg_el_pt =
            dynamic_cast<AlgebraicElementBase*>(this);

          // If it is, set up node position (past and present) from algebraic
          // node position (past and present) from algebraic node update
          // node update function. This over-writes previous assingments that
          // were made based on the macro-element/FE representation.
          // NOTE: YES, THIS NEEDS TO BE CALLED REPEATEDLY IF THE NODE IS A
          // MEMBER OF MULTIPLE ELEMENTS: THEY ALL ASSIGN THE SAME NODAL
          // POSITIONS BUT WE NEED TO ADD THE REMESH INFO FOR *ALL* ROOT
          // ELEMENTS!
          if (alg_el_pt != 0)
          {
            // Build algebraic node update info for new node. This sets up
            // the node update data for all node update functions that are
            // shared by all nodes in the father element.
            alg_el_pt->setup_algebraic_node_update(
              node_pt(n), s_in_father, father_el_pt);
          }

          // If we have built the node and we are documenting our progress,
          // write the (hopefully consistent position) to the outputfile
          if ((!node_done) && (new_nodes_file.is_open()))
          {
            new_nodes_file << node_pt(n)->x(0) << std::endl;
          }
        } // End of loop over all nodes in element


        // If the element is a MacroElementNodeUpdateElement, set the update
        // parameters for the current element's nodes -- all this needs is
        // the vector of (pointers to the) geometric objects that affect the
        // MacroElement-based node update. This is the same as that in the
        // father element
        MacroElementNodeUpdateElementBase* father_m_el_pt =
          dynamic_cast<MacroElementNodeUpdateElementBase*>(father_el_pt);
        if (father_m_el_pt != 0)
        {
          // Get Vector of geometric objects from father (construct Vector
          // via copy operation)
          Vector<GeomObject*> geom_object_pt(father_m_el_pt->geom_object_pt());

          // Cast current element to MacroElementNodeUpdateElement:
          MacroElementNodeUpdateElementBase* m_el_pt =
            dynamic_cast<MacroElementNodeUpdateElementBase*>(this);

#ifdef PARANOID
          if (m_el_pt == 0)
          {
            std::string error_message =
              "Failed to cast to MacroElementNodeUpdateElementBase*\n";
            error_message +=
              "Strange -- if the father is a MacroElementNodeUpdateElement\n";
            error_message += "the son should be too....\n";

            throw OomphLibError(
              error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
          }
#endif
          // Build update info by passing Vector of geometric objects:
          // This sets the current element to be the update element for all
          // of the element's nodes. This is reversed if the element is
          // ever un-refined in the father element's rebuild_from_sons()
          // function which overwrites this assignment to avoid nasty
          // segmentation faults that occur when a node tries to update
          // itself via an element that no longer exists...
          m_el_pt->set_node_update_info(geom_object_pt);
        }

#ifdef OOMPH_HAS_MPI
        // Pass on non-halo proc ID
        Non_halo_proc_ID =
          tree_pt()->father_pt()->object_pt()->non_halo_proc_ID();
#endif

        // Is the new element an ElementWithMovingNodes?
        ElementWithMovingNodes* aux_el_pt =
          dynamic_cast<ElementWithMovingNodes*>(this);

        // Pass down the information re the method for the evaluation
        // of the shape derivatives
        if (aux_el_pt != 0)
        {
          ElementWithMovingNodes* aux_father_el_pt =
            dynamic_cast<ElementWithMovingNodes*>(father_el_pt);

#ifdef PARANOID
          if (aux_father_el_pt == 0)
          {
            std::string error_message =
              "Failed to cast to ElementWithMovingNodes*\n";
            error_message +=
              "Strange -- if the son is a ElementWithMovingNodes\n";
            error_message += "the father should be too....\n";

            throw OomphLibError(
              error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
          }
#endif

          // If evaluating the residuals by finite differences in the father
          // continue to do so in the child
          if (aux_father_el_pt
                ->are_dresidual_dnodal_coordinates_always_evaluated_by_fd())
          {
            aux_el_pt
              ->enable_always_evaluate_dresidual_dnodal_coordinates_by_fd();
          }

          aux_el_pt->method_for_shape_derivs() =
            aux_father_el_pt->method_for_shape_derivs();

          // If bypassing the evaluation of fill_in_jacobian_from_geometric_data
          // continue to do so
          if (aux_father_el_pt
                ->is_fill_in_jacobian_from_geometric_data_bypassed())
          {
            aux_el_pt->enable_bypass_fill_in_jacobian_from_geometric_data();
          }
        }


        // Now do further build (if any)
        further_build();

      } // Sanity check: Father element has been generated

    } // End of if element has not already been built

    // If element has already been built, say so
    else
    {
      was_already_built = true;
    }
  }

  //==========================================================================
  ///  Print corner nodes, use colour (default "BLACK")
  //==========================================================================
  void RefineableQElement<1>::output_corners(std::ostream& outfile,
                                             const std::string& colour) const
  {
    // Allocate storage for local coordinate s
    Vector<double> s(1);

    // Allocate storage for global coordinate of element vertex
    Vector<double> vertex(1);

    outfile << "ZONE I=2,J=2, C=" << colour << std::endl;

    // Left-hand vertex
    s[0] = -1.0;
    get_x(s, vertex);
    outfile << vertex[0] << " " << Number << std::endl;

    // Right-hand vertex
    s[0] = 1.0;
    get_x(s, vertex);
    outfile << vertex[0] << " " << Number << std::endl;

    outfile << "TEXT  CS = GRID, X = " << vertex[0]
            << ", HU = GRID, H = 0.01, AN = MIDCENTER, T =\"" << Number << "\""
            << std::endl;
  }

  //==========================================================================
  /// Check inter-element continuity of
  /// - nodal positions
  /// - (nodally) interpolated function values
  //==========================================================================
  void RefineableQElement<1>::check_integrity(double& max_error)
  {
    using namespace BinaryTreeNames;

    // Calculate number of nodes
    const unsigned n_node = nnode_1d();

    // Number of timesteps (including present) for which continuity is
    // to be checked
    const unsigned n_time = 1;

    // Initialise errors
    max_error = 0.0;
    Vector<double> max_error_x(1, 0.0);
    double max_error_val = 0.0;

    // Initialise vector of element edges in the current element
    Vector<int> edge_in_current(2);
    edge_in_current[0] = L;
    edge_in_current[1] = R;

    // Loop over both edges
    for (unsigned edge_counter = 0; edge_counter < 2; edge_counter++)
    {
      // Allocate storage for the fractional position (in the current
      // element) of the node in the direction of s[0]
      Vector<double> s_fraction(1);

      // Allocate storage for the location of a node in terms of the local
      // coordinate of the current element
      Vector<double> s_in_current(1);

      // Allocate storage for the location of a node in terms of the local
      // coordinate of the neighbouring element
      Vector<double> s_in_neighbour(1);

      // Allocate storage for edge in neighbouring element
      int edge_in_neighbour;

      // Allocate storage for difference in size between current and
      // neighbouring element
      int diff_level;

      // Allocate storage for flag indicating if the node is not in the same
      // binary tree
      bool in_neighbouring_tree;

      // Calculate the local coordinate and the local coordinate as viewed
      // from the neighbour

      // Allocate storage for the pointer to the neighbouring element
      // (using its binary tree representation)
      BinaryTree* neighbour_pt;

      // Find pointer to neighbouring element along the edge in question and
      // calculate the local coordinate of the node within that element,
      // s_in_neighbour
      neighbour_pt =
        binary_tree_pt()->gteq_edge_neighbour(edge_in_current[edge_counter],
                                              s_in_neighbour,
                                              edge_in_neighbour,
                                              diff_level,
                                              in_neighbouring_tree);

      // If neighbour exists and has existing nodes
      if ((neighbour_pt != 0) && (neighbour_pt->object_pt()->nodes_built()))
      {
        // Need to exclude periodic nodes from this check.
        // There are only periodic nodes if we are in a neighbouring tree...
        bool is_periodic = false;
        if (in_neighbouring_tree)
        {
          // Is it periodic?
          is_periodic = tree_pt()->root_pt()->is_neighbour_periodic(
            edge_in_current[edge_counter]);
        }

        // Allocate storage for pointer to the local node
        Node* local_node_pt = 0;

        switch (edge_counter)
        {
          case 0:
            // Local fraction of node
            s_fraction[0] = 0.0;
            // Get pointer to local node
            local_node_pt = node_pt(0);
            break;

          case 1:
            // Local fraction of node
            s_fraction[0] = 1.0;
            // Get pointer to local node
            local_node_pt = node_pt(n_node - 1);
            break;
        }

        // Evaluate the local coordinate of the node in the current element
        s_in_current[0] = -1.0 + 2.0 * s_fraction[0];

        // NOTE: We have already calculated the local coordinate of the node
        // in the neighbouring element above when calling gteq_edge_neighbour()

        // Loop over timesteps
        for (unsigned t = 0; t < n_time; t++)
        {
          // Allocate storage for the nodal position in the neighbouring element
          Vector<double> x_in_neighbour(1);

          // Get the nodal position from the neighbouring element
          neighbour_pt->object_pt()->interpolated_x(
            t, s_in_neighbour, x_in_neighbour);

          // Check error only if the node is NOT periodic
          if (is_periodic == false)
          {
            // Find the spatial error
            double err = std::fabs(local_node_pt->x(t, 0) - x_in_neighbour[0]);

            // If it's bigger than our tolerance, say so
            if (err > 1e-9)
            {
              oomph_info << "errx " << err << " " << t << " "
                         << local_node_pt->x(t, 0) << " " << x_in_neighbour[0]
                         << std::endl;

              oomph_info << "at " << local_node_pt->x(0) << std::endl;
            }

            // If it's bigger than the previous max error, it is the
            // new max error!
            if (err > max_error_x[0])
            {
              max_error_x[0] = err;
            }
          }
          // Allocate storage for the nodal values in the neighbouring element
          Vector<double> values_in_neighbour;

          // Get the values from neighbouring element. NOTE: Number of values
          // gets set by routine (because in general we don't know how many
          // interpolated values a certain element has.
          neighbour_pt->object_pt()->get_interpolated_values(
            t, s_in_neighbour, values_in_neighbour);

          // Allocate storage for the nodal values in the current element
          Vector<double> values_in_current;

          // Get the values in current element
          get_interpolated_values(t, s_in_current, values_in_current);

          // Now figure out how many continuously interpolated values there are
          const unsigned num_val =
            neighbour_pt->object_pt()->ncont_interpolated_values();

          // Check error
          for (unsigned ival = 0; ival < num_val; ival++)
          {
            // Find the spatial error
            double err =
              std::fabs(values_in_current[ival] - values_in_neighbour[ival]);

            // If it's bigger than our tolerance, say so
            if (err > 1.0e-10)
            {
              oomph_info << local_node_pt->x(0) << "\n# "
                         << "erru " << err << " " << ival << " "
                         << get_node_number(local_node_pt) << " "
                         << values_in_current[ival] << " "
                         << values_in_neighbour[ival] << std::endl;
            }

            // If it's bigger than the previous max error, it is the
            // new max error!
            if (err > max_error_val)
            {
              max_error_val = err;
            }
          }
        } // End of loop over timesteps
      } // End of if neighbour exists and has existing nodes
    } // End of loop over edges

    // Update max_error if necessary
    max_error = max_error_x[0];
    if (max_error_val > max_error)
    {
      max_error = max_error_val;
    }

    // Output max_error information
    if (max_error > 1e-9)
    {
      oomph_info << "\n#------------------------------------ \n#Max error ";
      oomph_info << max_error_x[0] << " " << max_error_val << std::endl;
      oomph_info << "#------------------------------------ \n " << std::endl;
    }
  }

  //==========================================================================
  /// Static matrix for coincidence between son nodal points and father
  /// boundaries
  //==========================================================================
  std::map<unsigned, DenseMatrix<int>> RefineableQElement<1>::Father_bound;

} // namespace oomph
