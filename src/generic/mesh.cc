// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
// Non-inline member functions for general mesh classes

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

#include <algorithm>
#include <limits.h>
#include <typeinfo>


// oomph-lib headers
#include "oomph_utilities.h"
#include "mesh.h"
#include "problem.h"
#include "elastic_problems.h"
#include "refineable_mesh.h"
#include "triangle_mesh.h"
#include "shape.h"

namespace oomph
{
  //======================================================
  /// The Steady Timestepper
  //======================================================
  Steady<0> Mesh::Default_TimeStepper;


  //=======================================================================
  /// Static boolean flag to control warning about mesh level timesteppers
  //=======================================================================
  bool Mesh::Suppress_warning_about_empty_mesh_level_time_stepper_function =
    false;

  //=======================================================================
  /// Merge meshes.
  /// Note: This simply merges the meshes' elements and nodes (ignoring
  /// duplicates; no boundary information etc. is created).
  //=======================================================================
  void Mesh::merge_meshes(const Vector<Mesh*>& sub_mesh_pt)
  {
    // No boundary lookup scheme is set up for the combined mesh
    Lookup_for_elements_next_boundary_is_setup = false;

    // Number of submeshes
    unsigned nsub_mesh = sub_mesh_pt.size();

    // Initialise element, node and boundary counters for global mesh
    unsigned long n_element = 0;
    unsigned long n_node = 0;
    unsigned n_bound = 0;

    // Loop over submeshes and get total number of elements, nodes and
    // boundaries
    for (unsigned imesh = 0; imesh < nsub_mesh; imesh++)
    {
      n_element += sub_mesh_pt[imesh]->nelement();
      n_node += sub_mesh_pt[imesh]->nnode();
      n_bound += sub_mesh_pt[imesh]->nboundary();
    }

    // Reserve storage for element and node pointers
    Element_pt.clear();
    Element_pt.reserve(n_element);
    Node_pt.clear();
    Node_pt.reserve(n_node);

    // Resize vector of vectors of nodes
    Boundary_node_pt.clear();
    Boundary_node_pt.resize(n_bound);

    // Sets of pointers to elements and nodes (to exlude duplicates -- they
    // shouldn't occur anyway but if they do, they must only be added
    // once in the global mesh to avoid trouble in the timestepping)
    std::set<GeneralisedElement*> element_set_pt;
    std::set<Node*> node_set_pt;

    // Counter for total number of boundaries in all the submeshes
    unsigned ibound_global = 0;
    // Loop over the number of submeshes
    for (unsigned imesh = 0; imesh < nsub_mesh; imesh++)
    {
      // Loop over the elements of the submesh and add to vector
      // duplicates are ignored
      unsigned nel_before = 0;
      unsigned long n_element = sub_mesh_pt[imesh]->nelement();
      for (unsigned long e = 0; e < n_element; e++)
      {
        GeneralisedElement* el_pt = sub_mesh_pt[imesh]->element_pt(e);
        element_set_pt.insert(el_pt);
        // Was it a duplicate?
        unsigned nel_now = element_set_pt.size();
        if (nel_now == nel_before)
        {
          std::ostringstream warning_stream;
          warning_stream << "WARNING: " << std::endl
                         << "Element " << e << " in submesh " << imesh
                         << " is a duplicate \n and was ignored when assembling"
                         << " combined mesh." << std::endl;
          OomphLibWarning(warning_stream.str(),
                          "Mesh::Mesh(const Vector<Mesh*>&)",
                          OOMPH_EXCEPTION_LOCATION);
        }
        else
        {
          Element_pt.push_back(el_pt);
        }
        nel_before = nel_now;
      }

      // Loop over the nodes of the submesh and add to vector
      // duplicates are ignored
      unsigned nnod_before = 0;
      unsigned long n_node = sub_mesh_pt[imesh]->nnode();
      for (unsigned long n = 0; n < n_node; n++)
      {
        Node* nod_pt = sub_mesh_pt[imesh]->node_pt(n);
        node_set_pt.insert(nod_pt);
        // Was it a duplicate?
        unsigned nnod_now = node_set_pt.size();
        if (nnod_now == nnod_before)
        {
          std::ostringstream warning_stream;
          warning_stream
            << "WARNING: " << std::endl
            << "Node " << n << " in submesh " << imesh
            << " is a duplicate \n and was ignored when assembling "
            << "combined mesh." << std::endl;
          OomphLibWarning(warning_stream.str(),
                          "Mesh::Mesh(const Vector<Mesh*>&)",
                          OOMPH_EXCEPTION_LOCATION);
        }
        else
        {
          Node_pt.push_back(nod_pt);
        }
        nnod_before = nnod_now;
      }

      // Loop over the boundaries of the submesh
      unsigned n_bound = sub_mesh_pt[imesh]->nboundary();
      for (unsigned ibound = 0; ibound < n_bound; ibound++)
      {
        // Loop over the number of nodes on the boundary and add to the
        // global vector
        unsigned long n_bound_node = sub_mesh_pt[imesh]->nboundary_node(ibound);
        for (unsigned long n = 0; n < n_bound_node; n++)
        {
          Boundary_node_pt[ibound_global].push_back(
            sub_mesh_pt[imesh]->boundary_node_pt(ibound, n));
        }
        // Increase the number of the global boundary counter
        ibound_global++;
      }
    } // End of loop over submeshes
  }


  //========================================================
  /// Remove the information about nodes stored on the
  /// b-th boundary of the mesh
  //========================================================
  void Mesh::remove_boundary_nodes(const unsigned& b)
  {
    // Loop over all the nodes on the boundary and call
    // their remove_from_boundary function
    unsigned n_boundary_node = Boundary_node_pt[b].size();
    for (unsigned n = 0; n < n_boundary_node; n++)
    {
      boundary_node_pt(b, n)->remove_from_boundary(b);
    }
    // Clear the storage
    Boundary_node_pt[b].clear();
  }

  //=================================================================
  /// Remove all information about mesh boundaries
  //================================================================
  void Mesh::remove_boundary_nodes()
  {
    // Loop over each boundary call remove_boundary_nodes
    unsigned n_bound = Boundary_node_pt.size();
    for (unsigned b = 0; b < n_bound; b++)
    {
      remove_boundary_nodes(b);
    }
    // Clear the storage
    Boundary_node_pt.clear();
  }

  //============================================================
  /// Remove the node node_pt from the b-th boundary of the mesh
  /// This function also removes the information from the Node
  /// itself
  //===========================================================
  void Mesh::remove_boundary_node(const unsigned& b, Node* const& node_pt)
  {
    // Find the location of the node in the boundary
    Vector<Node*>::iterator it = std::find(
      Boundary_node_pt[b].begin(), Boundary_node_pt[b].end(), node_pt);
    // If the node is on this boundary
    if (it != Boundary_node_pt[b].end())
    {
      // Remove the node from the mesh's list of boundary nodes
      Boundary_node_pt[b].erase(it);
      // Now remove the node's boundary information
      node_pt->remove_from_boundary(b);
    }
    // If not do nothing
  }


  //========================================================
  /// Add the node node_pt to the b-th boundary of the mesh
  /// This function also sets the boundary information in the
  /// Node itself
  //=========================================================
  void Mesh::add_boundary_node(const unsigned& b, Node* const& node_pt)
  {
    // Tell the node that it's on boundary b.
    // At this point, if the node is not a BoundaryNode, the function
    // should throw an exception.
    node_pt->add_to_boundary(b);

    // Get the size of the Boundary_node_pt vector
    unsigned nbound_node = Boundary_node_pt[b].size();
    bool node_already_on_this_boundary = false;
    // Loop over the vector
    for (unsigned n = 0; n < nbound_node; n++)
    {
      // is the current node here already?
      if (node_pt == Boundary_node_pt[b][n])
      {
        node_already_on_this_boundary = true;
      }
    }

    // Add the base node pointer to the vector if it's not there already
    if (!node_already_on_this_boundary)
    {
      Boundary_node_pt[b].push_back(node_pt);
    }
  }


  //=======================================================
  /// Update nodal positions in response to changes in the domain shape.
  /// Uses the FiniteElement::get_x(...) function for FiniteElements
  /// and doesn't do anything for other element types.
  /// If a MacroElement pointer has been set for a FiniteElement,
  /// the MacroElement representation is used to update the
  /// nodal positions; if not get_x(...) uses the FE interpolation
  /// and thus leaves the nodal positions unchanged.
  /// Virtual, so it can be overloaded by specific meshes,
  /// such as AlgebraicMeshes or SpineMeshes.
  /// Generally, this function updates the position of all nodes
  /// in response to changes in the boundary position. For
  /// SolidNodes it only applies the update to those SolidNodes
  /// whose position is determined by the boundary position, unless
  /// the bool flag is set to true.
  //========================================================
  void Mesh::node_update(const bool& update_all_solid_nodes)
  {
    // Get the current time
    double t_start = TimingHelpers::timer();

#ifdef PARANOID
#ifdef OOMPH_HAS_MPI
    // Paranoid check to throw an error if node update is called for elements
    // with nonuniformly spaced nodes for which some masters are 'external'
    for (unsigned long n = 0; n < nnode(); n++)
    {
      Node* nod_pt = Node_pt[n];
      if (nod_pt->is_hanging())
      {
        // Loop over master nodes
        unsigned nmaster = nod_pt->hanging_pt()->nmaster();
        for (unsigned imaster = 0; imaster < nmaster; imaster++)
        {
          // Get pointer to master node
          Node* master_nod_pt = nod_pt->hanging_pt()->master_node_pt(imaster);

          // Get vector of all external halo nodes
          Vector<Node*> external_halo_node_pt;
          get_external_halo_node_pt(external_halo_node_pt);

          // Search the external halo storage for this node
          Vector<Node*>::iterator it = std::find(external_halo_node_pt.begin(),
                                                 external_halo_node_pt.end(),
                                                 master_nod_pt);

          // Check if the node was found
          if (it != external_halo_node_pt.end())
          {
            // Throw error becase node update won't work
            // It's ok to throw an error here because this function is
            // overloaded for Algebraic and MacroElementNodeUpdate
            // Meshes. This is only a problem for meshes of ordinary
            // nodes.
            std::ostringstream err_stream;

            err_stream << "Calling node_update() for a mesh which contains"
                       << std::endl
                       << "master nodes which live in the external storage."
                       << std::endl
                       << "These nodes do not belong to elements on this"
                       << std::endl
                       << "processor and therefore cannot be updated locally."
                       << std::endl;

            throw OomphLibError(err_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }
    }
    // If we get to here then none of the masters of any of the nodes in the
    // mesh live in the external storage, so we'll be fine if we carry on.
#endif
#endif

    /// Local and global (Eulerian) coordinate
    Vector<double> s;
    Vector<double> r;

    // NB: This repeats nodes a lot - surely it would be
    // quicker to modify it so that it only does each node once,
    // particularly in the update_all_solid_nodes=true case?
    // Create a map to indicate whether or not we've updated a node already
    std::map<Node*, bool> has_node_been_updated;

    // How many nodes are there?
    unsigned n_node = nnode();

    // Loop over all Nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // Get pointer to node
      Node* nod_pt = node_pt(n);

      // Initialise the boolean value associated with this node
      has_node_been_updated[nod_pt] = false;
    }

    // Loop over all elements
    unsigned nel = nelement();
    for (unsigned e = 0; e < nel; e++)
    {
      // Try to cast to FiniteElement
      FiniteElement* el_pt = dynamic_cast<FiniteElement*>(element_pt(e));

      // If it's a finite element we can proceed: FiniteElements have
      // nodes and a get_x() function
      if (el_pt != 0)
      {
        // Find out dimension of element = number of local coordinates
        unsigned ndim_el = el_pt->dim();
        s.resize(ndim_el);

        // Loop over nodal points
        unsigned n_node = el_pt->nnode();
        for (unsigned j = 0; j < n_node; j++)
        {
          // Get pointer to node
          Node* nod_pt = el_pt->node_pt(j);

          // Get spatial dimension of node
          unsigned ndim_node = nod_pt->ndim();
          r.resize(ndim_node);

          // For non-hanging nodes
          if (!(nod_pt->is_hanging()))
          {
            // If we've not dealt with this Node yet
            if (!has_node_been_updated[nod_pt])
            {
              // Get the position of the node
              el_pt->local_coordinate_of_node(j, s);

              // Get new position
              el_pt->get_x(s, r);

              // Try to cast to SolidNode
              SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);

              // Loop over coordinate directions
              for (unsigned i = 0; i < ndim_node; i++)
              {
                // It's a SolidNode:
                if (solid_node_pt != 0)
                {
                  // only do it if explicitly requested!
                  if (update_all_solid_nodes)
                  {
                    solid_node_pt->x(i) = r[i];
                  }
                }
                // Not a SolidNode: Definitely update
                else
                {
                  nod_pt->x(i) = r[i];
                }
              }

              // Indicate that we're done with this node, regardless of whether
              // or not it even needed updating
              has_node_been_updated[nod_pt] = true;
            } // if (!has_node_been_updated[nod_pt])
          } // if (!(nod_pt->is_hanging()))
        } // for (unsigned j=0;j<n_node;j++)
      } // if (el_pt!=0)
    } // for (unsigned e=0;e<nel;e++)

    // Now update the external halo nodes before we adjust the positions of the
    // hanging nodes incase any are masters of local nodes
#ifdef OOMPH_HAS_MPI
    // Loop over all external halo nodes with other processors
    // and update them
    for (std::map<unsigned, Vector<Node*>>::iterator it =
           External_halo_node_pt.begin();
         it != External_halo_node_pt.end();
         it++)
    {
      // Get vector of external halo nodes
      Vector<Node*> ext_halo_node_pt = (*it).second;
      unsigned nnod = ext_halo_node_pt.size();
      for (unsigned j = 0; j < nnod; j++)
      {
        ext_halo_node_pt[j]->node_update();
      }
    }
#endif

    // Now loop over hanging nodes and adjust their position
    // in line with their hanging node constraints
    for (unsigned long n = 0; n < n_node; n++)
    {
      Node* nod_pt = Node_pt[n];
      if (nod_pt->is_hanging())
      {
        // Get spatial dimension of node
        unsigned ndim_node = nod_pt->ndim();

        // Initialise
        for (unsigned i = 0; i < ndim_node; i++)
        {
          nod_pt->x(i) = 0.0;
        }

        // Loop over master nodes
        unsigned nmaster = nod_pt->hanging_pt()->nmaster();
        for (unsigned imaster = 0; imaster < nmaster; imaster++)
        {
          // Loop over directions
          for (unsigned i = 0; i < ndim_node; i++)
          {
            nod_pt->x(i) +=
              nod_pt->hanging_pt()->master_node_pt(imaster)->x(i) *
              nod_pt->hanging_pt()->master_weight(imaster);
          }
        }
      }
    }

    // Loop over all nodes again and execute auxiliary node update
    // function
    for (unsigned long n = 0; n < n_node; n++)
    {
      Node_pt[n]->perform_auxiliary_node_update_fct();
    }

    // Tell the user how long it's taken
    oomph_info << "Time taken to update nodal positions [sec]: "
               << TimingHelpers::timer() - t_start << std::endl;
  }


  //=======================================================
  /// Reorder nodes in the order in which they are
  /// encountered when stepping through the elements
  //========================================================
  void Mesh::reorder_nodes(const bool& use_old_ordering)
  {
    // Create storage for the reordered nodes
    Vector<Node*> reordering;

    // Get the reordered nodes (without altering the mesh's node vector)
    get_node_reordering(reordering, use_old_ordering);

    // Get the number of nodes in the mesh
    unsigned n_node = nnode();

    // Loop over all of the nodes
    for (unsigned i = 0; i < n_node; i++)
    {
      // Replace the Mesh's i-th node pointer with the reordered node pointer
      node_pt(i) = reordering[i];
    }
  } // End of reorder_nodes

  //=======================================================
  /// Get a vector of the nodes in the order in which they are encountered
  /// when stepping through the elements (similar to reorder_nodes() but
  /// without changing the mesh's node vector).
  //========================================================
  void Mesh::get_node_reordering(Vector<Node*>& reordering,
                                 const bool& use_old_ordering) const
  {
    // If the user wants to use the original order
    if (use_old_ordering)
    {
      // Setup map to check if nodes have been done yet
      std::map<Node*, bool> done;

      // Loop over all nodes
      unsigned nnod = nnode();

      // Initialise the vector
      reordering.assign(nnod, 0);

      // Return immediately if there are no nodes: Note assumption:
      // Either all the elements' nodes stored here or none. If only a subset
      // is stored in the Node_pt vector we'll get a range checking error below
      // (only if run with range checking, of course).
      if (nnod == 0)
      {
        // Return immediately
        return;
      }

      // Loop over the nodes in the mesh
      for (unsigned j = 0; j < nnod; j++)
      {
        // Indicate whether or not the node has been swapped
        done[node_pt(j)] = false;
      }

      // Initialise counter for number of nodes
      unsigned long count = 0;

      // Get the number of elements in the mesh
      unsigned nel = nelement();

      // Loop over all elements
      for (unsigned e = 0; e < nel; e++)
      {
        // Make sure the e-th element is a FiniteElement (or derived) class
        // object
        FiniteElement* el_pt =
          checked_dynamic_cast<FiniteElement*>(element_pt(e));

        // Get the number of nodes in this element
        unsigned nnod = el_pt->nnode();

        // Loop over nodes in element
        for (unsigned j = 0; j < nnod; j++)
        {
          // Get a pointer to the j-th node in the element
          Node* nod_pt = el_pt->node_pt(j);

          // Has node been done yet?
          if (!done[nod_pt])
          {
            // Insert into node vector. NOTE: If you get a seg fault/range
            // checking error here then you probably haven't added all the
            // elements' nodes to the Node_pt vector -- this is most likely to
            // arise in the case of meshes of face elements (though they usually
            // don't store the nodes at all so if you have any problems here
            // there's something unusual/not quite right in any case... For this
            // reason we don't range check here by default (not even under
            // paranoia) but force you turn on proper (costly) range checking to
            // track this down...
            reordering[count] = nod_pt;

            // Indicate that the node has been done
            done[nod_pt] = true;

            // Increase counter
            count++;
          }
        } // for (unsigned j=0;j<nnod;j++)
      } // for (unsigned e=0;e<nel;e++)

      // Sanity check
      if (count != nnod)
      {
        // Create an error message
        std::string error_message = "Trouble: Number of nodes hasn't stayed ";

        // Finish off the message
        error_message += "constant during reordering!\n";

        // Throw an error
        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
    else
    {
      // Copy node vector out
      unsigned n_node = nnode();

      // Resize the node ordering vector
      reordering.resize(n_node);

      // Loop over the nodes
      for (unsigned i = 0; i < n_node; i++)
      {
        // Assign the i-th node pointer entry
        reordering[i] = node_pt(i);
      }

      // Now sort the nodes lexicographically
      std::sort(reordering.begin(),
                reordering.end(),
                &NodeOrdering::node_global_position_comparison);
    } // if (use_old_ordering)
  } // End of get_node_reordering


  //========================================================
  /// Virtual Destructor to clean up all memory
  //========================================================
  Mesh::~Mesh()
  {
    // Free the nodes first
    // Loop over the nodes in reverse order
    unsigned long Node_pt_range = Node_pt.size();
    for (unsigned long i = Node_pt_range; i > 0; i--)
    {
      delete Node_pt[i - 1];
      Node_pt[i - 1] = 0;
    }

    // Free the elements
    // Loop over the elements in reverse order
    unsigned long Element_pt_range = Element_pt.size();
    for (unsigned long i = Element_pt_range; i > 0; i--)
    {
      delete Element_pt[i - 1];
      Element_pt[i - 1] = 0;
    }

    // Wipe the storage for all externally-based elements and delete halos
    delete_all_external_storage();
  }

  //========================================================
  /// Assign (global) equation numbers to the nodes
  //========================================================
  unsigned long Mesh::assign_global_eqn_numbers(Vector<double*>& Dof_pt)
  {
    // Find out the current number of equations
    unsigned long equation_number = Dof_pt.size();

    // Loop over the nodes and call their assigment functions
    unsigned long nnod = Node_pt.size();

    for (unsigned long i = 0; i < nnod; i++)
    {
      Node_pt[i]->assign_eqn_numbers(equation_number, Dof_pt);
    }

    // Loop over the elements and number their internals
    unsigned long nel = Element_pt.size();
    for (unsigned long i = 0; i < nel; i++)
    {
      Element_pt[i]->assign_internal_eqn_numbers(equation_number, Dof_pt);
    }

    // Return the total number of equations
    return (equation_number);
  }

  //========================================================
  /// Function to describe the dofs of the Mesh. The ostream
  /// specifies the output stream to which the description
  /// is written; the string stores the currently
  /// assembled output that is ultimately written to the
  /// output stream by Data::describe_dofs(...); it is typically
  /// built up incrementally as we descend through the
  /// call hierarchy of this function when called from
  /// Problem::describe_dofs(...)
  //========================================================
  void Mesh::describe_dofs(std::ostream& out,
                           const std::string& current_string) const
  {
    // Loop over the nodes and call their classification functions
    unsigned long nnod = Node_pt.size();
    for (unsigned long i = 0; i < nnod; i++)
    {
      std::stringstream conversion;
      conversion << " of Node " << i << current_string;
      std::string in(conversion.str());
      Node_pt[i]->describe_dofs(out, in);
    }

    // Loop over the elements and classify.
    unsigned long nel = Element_pt.size();
    for (unsigned long i = 0; i < nel; i++)
    {
      std::stringstream conversion;
      conversion << " in Element " << i << " [" << typeid(*Element_pt[i]).name()
                 << "] " << current_string;
      std::string in(conversion.str());
      Element_pt[i]->describe_dofs(out, in);
    }
  }

  //========================================================
  /// Function to describe the local dofs of the elements. The ostream
  /// specifies the output stream to which the description
  /// is written; the string stores the currently
  /// assembled output that is ultimately written to the
  /// output stream by Data::describe_dofs(...); it is typically
  /// built up incrementally as we descend through the
  /// call hierarchy of this function when called from
  /// Problem::describe_dofs(...)
  //========================================================
  void Mesh::describe_local_dofs(std::ostream& out,
                                 const std::string& current_string) const
  {
    // Now loop over the elements and describe local dofs
    unsigned long nel = Element_pt.size();
    for (unsigned long i = 0; i < nel; i++)
    {
      std::stringstream conversion;
      conversion << " in Element" << i << " [" << typeid(*Element_pt[i]).name()
                 << "] " << current_string;
      std::string in(conversion.str());
      Element_pt[i]->describe_local_dofs(out, in);
    }
  }


  //========================================================
  /// Assign local equation numbers in all elements
  //========================================================
  void Mesh::assign_local_eqn_numbers(const bool& store_local_dof_pt)
  {
    // Now loop over the elements and assign local equation numbers
    unsigned long Element_pt_range = Element_pt.size();
    for (unsigned long i = 0; i < Element_pt_range; i++)
    {
      Element_pt[i]->assign_local_eqn_numbers(store_local_dof_pt);
    }
  }

  //========================================================
  /// Self-test: Check elements and nodes. Return 0 for OK
  //========================================================
  unsigned Mesh::self_test()
  {
    // Initialise
    bool passed = true;

    // Check the mesh for repeated nodes (issues its own error message)
    if (0 != check_for_repeated_nodes()) passed = false;

    // hierher -- re-enable once problem with Hermite elements has been
    // resolved.
    //  // Check if there are any inverted elements
    //  bool mesh_has_inverted_elements=false;
    //  check_inverted_elements(mesh_has_inverted_elements);
    //  if (mesh_has_inverted_elements)
    //   {
    //    passed=false;
    //    oomph_info << "\n ERROR: Mesh has inverted elements\n"
    //               << "     Run Mesh::check_inverted_elements(...) with"
    //               << "     with output stream to find out which elements are"
    //               << "     inverted.\n";
    //   }

    // Loop over the elements, check for duplicates and do self test
    std::set<GeneralisedElement*> element_set_pt;
    unsigned long Element_pt_range = Element_pt.size();
    for (unsigned long i = 0; i < Element_pt_range; i++)
    {
      if (Element_pt[i]->self_test() != 0)
      {
        passed = false;
        oomph_info << "\n ERROR: Failed Element::self_test() for element i="
                   << i << std::endl;
      }
      // Add to set (which ignores duplicates):
      element_set_pt.insert(Element_pt[i]);
    }

    // Check for duplicates:
    if (element_set_pt.size() != Element_pt_range)
    {
      oomph_info << "ERROR:  " << Element_pt_range - element_set_pt.size()
                 << " duplicate elements were encountered in mesh!"
                 << std::endl;
      passed = false;
    }


    // Loop over the nodes, check for duplicates and do self test
    std::set<Node*> node_set_pt;
    unsigned long Node_pt_range = Node_pt.size();
    for (unsigned long i = 0; i < Node_pt_range; i++)
    {
      if (Node_pt[i]->self_test() != 0)
      {
        passed = false;
        oomph_info << "\n ERROR: Failed Node::self_test() for node i=" << i
                   << std::endl;
      }
      // Add to set (which ignores duplicates):
      node_set_pt.insert(Node_pt[i]);
    }

    // Check for duplicates:
    if (node_set_pt.size() != Node_pt_range)
    {
      oomph_info << "ERROR:  " << Node_pt_range - node_set_pt.size()
                 << " duplicate nodes were encountered in mesh!" << std::endl;
      passed = false;
    }

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


  //========================================================
  ///  Check for inverted elements and report outcome
  /// in boolean variable. This visits all elements at their
  /// integration points and checks if the Jacobian of the
  /// mapping between local and global coordinates is positive --
  /// using the same test that would be carried out (but only in PARANOID
  /// mode) during the assembly of the elements' Jacobian matrices.
  /// Inverted elements are output in inverted_element_file (if the
  /// stream is open).
  //========================================================
  void Mesh::check_inverted_elements(bool& mesh_has_inverted_elements,
                                     std::ofstream& inverted_element_file)
  {
    // Initialise flag
    mesh_has_inverted_elements = false;

    // Suppress output while checking for inverted elements
    bool backup =
      FiniteElement::Suppress_output_while_checking_for_inverted_elements;
    FiniteElement::Suppress_output_while_checking_for_inverted_elements = true;

    // Loop over all elements
    unsigned nelem = nelement();
    for (unsigned e = 0; e < nelem; e++)
    {
      FiniteElement* el_pt = finite_element_pt(e);

      // Only check for finite elements
      if (el_pt != 0)
      {
        // Find out number of nodes and local coordinates in the element
        unsigned n_node = el_pt->nnode();
        unsigned n_dim = el_pt->dim();
        unsigned ndim_node = el_pt->nodal_dimension();

        // Can't check Jacobian for elements in which nodal and elementa
        // dimensions don't match
        if (n_dim == ndim_node)
        {
          // Set up memory for the shape and test function and local coord
          Shape psi(n_node);
          DShape dpsidx(n_node, n_dim);
          Vector<double> s(n_dim);

          // Initialise element-level test
          bool is_inverted = false;

          unsigned n_intpt = el_pt->integral_pt()->nweight();
          for (unsigned ipt = 0; ipt < n_intpt; ipt++)
          {
            // Local coordinates
            for (unsigned i = 0; i < n_dim; i++)
            {
              s[i] = el_pt->integral_pt()->knot(ipt, i);
            }

            double J = 0;
            // Dummy assignment to keep gcc from complaining about
            // "set but unused".
            J += 0.0;
            try
            {
              // Call the derivatives of the shape functions and Jacobian
              J = el_pt->dshape_eulerian(s, psi, dpsidx);

              // If code is compiled without PARANOID setting, the
              // above call will simply return the negative Jacobian
              // without failing, so we need to check manually
#ifndef PARANOID
              try
              {
                el_pt->check_jacobian(J);
              }
              catch (OomphLibQuietException& error)
              {
                is_inverted = true;
              }
#endif
            }
            catch (OomphLibQuietException& error)
            {
              is_inverted = true;
            }
          }
          if (is_inverted)
          {
            mesh_has_inverted_elements = true;
            if (inverted_element_file.is_open())
            {
              el_pt->output(inverted_element_file);
            }
          }
        }
      }
    }
    // Reset
    FiniteElement::Suppress_output_while_checking_for_inverted_elements =
      backup;
  }


  //========================================================
  /// Nodes that have been marked as obsolete are removed
  /// from the mesh and the its boundaries. Returns vector
  /// of pointers to deleted nodes.
  //========================================================
  Vector<Node*> Mesh::prune_dead_nodes()
  {
    // Only copy the 'live' nodes across to new mesh
    //----------------------------------------------

    // New Vector of pointers to nodes
    Vector<Node*> new_node_pt;
    Vector<Node*> deleted_node_pt;

    // Loop over all nodes in mesh
    unsigned long n_node = nnode();
    for (unsigned long n = 0; n < n_node; n++)
    {
      // If the node still exists: Copy across
      if (!(Node_pt[n]->is_obsolete()))
      {
        new_node_pt.push_back(Node_pt[n]);
      }
      // Otherwise the Node is gone:
      // Delete it for good if it does not lie on a boundary
      // (if it lives on a boundary we have to remove it from
      // the boundary lookup schemes below)
      else
      {
        if (!(Node_pt[n]->is_on_boundary()))
        {
          deleted_node_pt.push_back(Node_pt[n]);
          delete Node_pt[n];
          Node_pt[n] = 0;
        }
      }
    }

    // Now update old vector by setting it equal to the new vector
    Node_pt = new_node_pt;


    // Boundaries
    //-----------

    // Only copy the 'live' nodes into new boundary node arrays
    //---------------------------------------------------------
    // Loop over the boundaries
    unsigned num_bound = nboundary();
    for (unsigned ibound = 0; ibound < num_bound; ibound++)
    {
      // New Vector of pointers to existent boundary nodes
      Vector<Node*> new_boundary_node_pt;

      // Loop over the boundary nodes
      unsigned long Nboundary_node = Boundary_node_pt[ibound].size();

      // Reserve contiguous memory for new vector of pointers
      // Must be equal in size to the number of nodes or less
      new_boundary_node_pt.reserve(Nboundary_node);

      for (unsigned long n = 0; n < Nboundary_node; n++)
      {
        // If node still exists: Copy across
        if (!(Boundary_node_pt[ibound][n]->is_obsolete()))
        {
          new_boundary_node_pt.push_back(Boundary_node_pt[ibound][n]);
        }
        // Otherwise Node is gone: Delete it for good
        else
        {
          // The node may lie on multiple boundaries, so remove the node
          // from the current boundary
          Boundary_node_pt[ibound][n]->remove_from_boundary(ibound);

          // Now if the node is no longer on any boundaries, delete it
          if (!Boundary_node_pt[ibound][n]->is_on_boundary())
          {
            deleted_node_pt.push_back(
              dynamic_cast<Node*>(Boundary_node_pt[ibound][n]));

            delete Boundary_node_pt[ibound][n];
          }
        }
      }

      // Update the old vector by setting it equal to the new vector
      Boundary_node_pt[ibound] = new_boundary_node_pt;

    } // End of loop over boundaries

    // Tell us who you deleted
    return deleted_node_pt;
  }


  //========================================================
  /// Output function for the mesh boundaries
  ///
  /// Loop over all boundaries and dump out the coordinates
  /// of the points on the boundary (in individual tecplot
  /// zones)
  //========================================================
  void Mesh::output_boundaries(std::ostream& outfile)
  {
    // Loop over the boundaries
    unsigned num_bound = nboundary();
    for (unsigned long ibound = 0; ibound < num_bound; ibound++)
    {
      unsigned nnod = Boundary_node_pt[ibound].size();
      if (nnod > 0)
      {
        outfile << "ZONE T=\"boundary" << ibound << "\"\n";

        for (unsigned inod = 0; inod < nnod; inod++)
        {
          Boundary_node_pt[ibound][inod]->output(outfile);
        }
      }
    }
  }


  //===================================================================
  /// Dump function for the mesh class.
  /// Loop over all nodes and elements and dump them
  //===================================================================
  void Mesh::dump(std::ofstream& dump_file, const bool& use_old_ordering) const
  {
    // Get a reordering of the nodes so that the dump file is in a standard
    // ordering regardless of the sequence of mesh refinements etc.
    Vector<Node*> reordering;
    this->get_node_reordering(reordering, use_old_ordering);

    // Find number of nodes
    unsigned long Node_pt_range = this->nnode();

    // Doc # of nodes
    dump_file << Node_pt_range << " # number of nodes " << std::endl;

    // Loop over all the nodes and dump their data
    for (unsigned nd = 0; nd < Node_pt_range; nd++)
    {
      reordering[nd]->dump(dump_file);
    }

    // Loop over elements and deal with internal data
    unsigned n_element = this->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      GeneralisedElement* el_pt = this->element_pt(e);
      unsigned n_internal = el_pt->ninternal_data();
      if (n_internal > 0)
      {
        dump_file << n_internal
                  << " # number of internal Data items in element " << e
                  << std::endl;
        for (unsigned i = 0; i < n_internal; i++)
        {
          el_pt->internal_data_pt(i)->dump(dump_file);
        }
      }
    }
  }


  //=======================================================
  /// Read solution from restart file
  //=======================================================
  void Mesh::read(std::ifstream& restart_file)
  {
    std::string input_string;

    // Reorder the nodes within the mesh's node vector
    // to establish a standard ordering regardless of the sequence
    // of mesh refinements etc
    this->reorder_nodes();

    // Read nodes

    // Find number of nodes
    unsigned long n_node = this->nnode();

    // Read line up to termination sign
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Check # of nodes:
    unsigned long check_n_node = atoi(input_string.c_str());
    if (check_n_node != n_node)
    {
      std::ostringstream error_stream;
      error_stream << "The number of nodes allocated " << n_node
                   << " is not the same as specified in the restart file "
                   << check_n_node << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Loop over the nodes
    for (unsigned long n = 0; n < n_node; n++)
    {
      /// Try to cast to elastic node
      SolidNode* el_node_pt = dynamic_cast<SolidNode*>(this->node_pt(n));
      if (el_node_pt != 0)
      {
        el_node_pt->read(restart_file);
      }
      else
      {
        this->node_pt(n)->read(restart_file);
      }
    }

    // Read internal data of elements:
    //--------------------------------
    // Loop over elements and deal with internal data
    unsigned n_element = this->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      GeneralisedElement* el_pt = this->element_pt(e);
      unsigned n_internal = el_pt->ninternal_data();
      if (n_internal > 0)
      {
        // Read line up to termination sign
        getline(restart_file, input_string, '#');

        // Ignore rest of line
        restart_file.ignore(80, '\n');

        // Check # of internals :
        unsigned long check_n_internal = atoi(input_string.c_str());
        if (check_n_internal != n_internal)
        {
          std::ostringstream error_stream;
          error_stream << "The number of internal data  " << n_internal
                       << " is not the same as specified in the restart file "
                       << check_n_internal << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        for (unsigned i = 0; i < n_internal; i++)
        {
          el_pt->internal_data_pt(i)->read(restart_file);
        }
      }
    }
  }


  //========================================================
  /// Output in paraview format into specified file.
  ///
  /// Breaks up each element into sub-elements for plotting
  /// purposes. We assume that all elements are of the same
  /// type (fct will break (in paranoid mode) if paraview
  /// output fcts of the elements are inconsistent).
  //========================================================
  void Mesh::output_paraview(std::ofstream& file_out,
                             const unsigned& nplot) const
  {
    // Change the scientific format so that E is used rather than e
    file_out.setf(std::ios_base::uppercase);

    // Decide how many elements there are to be plotted
    unsigned long number_of_elements = this->Element_pt.size();

    // Cast to finite element and return if cast fails.
    FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(0));

#ifdef PARANOID
    if (fe_pt == 0)
    {
      throw OomphLibError("Recast for FiniteElement failed for element 0!\n",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif


#ifdef PARANOID
    // Check if all elements have the same number of degrees of freedom,
    // if they don't, paraview will break
    unsigned el_zero_ndof = fe_pt->nscalar_paraview();
    for (unsigned i = 1; i < number_of_elements; i++)
    {
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));
      unsigned el_i_ndof = fe_pt->nscalar_paraview();
      if (el_zero_ndof != el_i_ndof)
      {
        std::stringstream error_stream;
        error_stream
          << "Element " << i << " has different number of degrees of freedom\n"
          << "than from previous elements, Paraview cannot handle this.\n"
          << "We suggest that the problem is broken up into submeshes instead."
          << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Make variables to hold the number of nodes and elements
    unsigned long number_of_nodes = 0;
    unsigned long total_number_of_elements = 0;

    // Loop over all the elements to find total number of plot points
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " failed" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      number_of_nodes += fe_pt->nplot_points_paraview(nplot);
      total_number_of_elements += fe_pt->nsub_elements_paraview(nplot);
    }


    // File Declaration
    //------------------

    // Insert the necessary lines plus header of file, and
    // number of nodes and elements
    file_out << "<?xml version=\"1.0\"?>\n"
             << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
             << "byte_order=\"LittleEndian\">\n"
             << "<UnstructuredGrid>\n"
             << "<Piece NumberOfPoints=\"" << number_of_nodes
             << "\" NumberOfCells=\"" << total_number_of_elements << "\">\n";


    // Point Data
    //-----------

    // Check the number of degrees of freedom
    unsigned ndof = fe_pt->nscalar_paraview();

    // Point data is going in here
    file_out << "<PointData ";

    // Insert just the first scalar name, since paraview reads everything
    // else after that as being of the same type. Get information from
    // first element.
    file_out << "Scalars=\"" << fe_pt->scalar_name_paraview(0) << "\">\n";

    // Loop over i scalar fields and j number of elements
    for (unsigned i = 0; i < ndof; i++)
    {
      file_out << "<DataArray type=\"Float32\" "
               << "Name=\"" << fe_pt->scalar_name_paraview(i) << "\" "
               << "format=\"ascii\""
               << ">\n";

      for (unsigned j = 0; j < number_of_elements; j++)
      {
        // Cast to FiniteElement and (in paranoid mode) check
        // if cast has failed.
        FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(j));

#ifdef PARANOID
        if (fe_pt == 0)
        {
          std::stringstream error_stream;
          error_stream << "Recast for element " << j << " failed" << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        fe_pt->scalar_value_paraview(file_out, i, nplot);
      }

      // Close of the DataArray
      file_out << "</DataArray>\n";
    }

    // Close off the PointData set
    file_out << "</PointData>\n";


    // Geometric Points
    //------------------

    file_out << "<Points>\n"
             << "<DataArray type=\"Float32\""
             << " NumberOfComponents=\""
             // This always has to be 3 for an unstructured grid
             << 3 << "\" "
             << "format=\"ascii\">\n";

    // Loop over all the elements to print their plot points
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " faild" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      fe_pt->output_paraview(file_out, nplot);
    }

    file_out << "</DataArray>\n"
             << "</Points>\n";


    // Cells
    //-------

    file_out
      << "<Cells>\n"
      << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    // Make counter for keeping track of all the local elements,
    // because Paraview requires global coordinates
    unsigned counter = 0;

    // Write connectivity with the local elements
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " faild" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      fe_pt->write_paraview_output_offset_information(file_out, nplot, counter);
    }

    file_out << "</DataArray>\n"
             << "<DataArray type=\"Int32\" "
             << "Name=\"offsets\" format=\"ascii\">\n";

    // Make variable that holds the current offset number
    unsigned offset_sum = 0;

    // Write the offset for the specific elements
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " failed" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      fe_pt->write_paraview_offsets(file_out, nplot, offset_sum);
    }

    file_out << "</DataArray>\n"
             << "<DataArray type=\"UInt8\" Name=\"types\">\n";

    // Loop over all elements to get the type that they have
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " failed" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      fe_pt->write_paraview_type(file_out, nplot);
    }

    file_out << "</DataArray>\n"
             << "</Cells>\n";


    // File Closure
    //-------------
    file_out << "</Piece>\n"
             << "</UnstructuredGrid>\n"
             << "</VTKFile>";
  }


  //========================================================
  /// Output in paraview format into specified file.
  ///
  /// Breaks up each element into sub-elements for plotting
  /// purposes. We assume that all elements are of the same
  /// type (fct will break (in paranoid mode) if paraview
  /// output fcts of the elements are inconsistent).
  //========================================================
  void Mesh::output_fct_paraview(
    std::ofstream& file_out,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt) const
  {
    // Change the scientific format so that E is used rather than e
    file_out.setf(std::ios_base::uppercase);

    // Decide how many elements there are to be plotted
    unsigned long number_of_elements = this->Element_pt.size();

    // Cast to finite element and return if cast fails.
    FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(0));

#ifdef PARANOID
    if (fe_pt == 0)
    {
      throw OomphLibError("Recast for FiniteElement failed for element 0!\n",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif


#ifdef PARANOID
    // Check if all elements have the same number of degrees of freedom,
    // if they don't, paraview will break
    unsigned el_zero_ndof = fe_pt->nscalar_paraview();
    for (unsigned i = 1; i < number_of_elements; i++)
    {
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));
      unsigned el_i_ndof = fe_pt->nscalar_paraview();
      if (el_zero_ndof != el_i_ndof)
      {
        std::stringstream error_stream;
        error_stream
          << "Element " << i << " has different number of degrees of freedom\n"
          << "than from previous elements, Paraview cannot handle this.\n"
          << "We suggest that the problem is broken up into submeshes instead."
          << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Make variables to hold the number of nodes and elements
    unsigned long number_of_nodes = 0;
    unsigned long total_number_of_elements = 0;

    // Loop over all the elements to find total number of plot points
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " failed" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      number_of_nodes += fe_pt->nplot_points_paraview(nplot);
      total_number_of_elements += fe_pt->nsub_elements_paraview(nplot);
    }


    // File Declaration
    //------------------

    // Insert the necessary lines plus header of file, and
    // number of nodes and elements
    file_out << "<?xml version=\"1.0\"?>\n"
             << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
             << "byte_order=\"LittleEndian\">\n"
             << "<UnstructuredGrid>\n"
             << "<Piece NumberOfPoints=\"" << number_of_nodes
             << "\" NumberOfCells=\"" << total_number_of_elements << "\">\n";


    // Point Data
    //-----------

    // Check the number of degrees of freedom
    unsigned ndof = fe_pt->nscalar_paraview();

    // Point data is going in here
    file_out << "<PointData ";

    // Insert just the first scalar name, since paraview reads everything
    // else after that as being of the same type. Get information from
    // first element.
    file_out << "Scalars=\"" << fe_pt->scalar_name_paraview(0) << "\">\n";

    // Loop over i scalar fields and j number of elements
    for (unsigned i = 0; i < ndof; i++)
    {
      file_out << "<DataArray type=\"Float32\" "
               << "Name=\"" << fe_pt->scalar_name_paraview(i) << "\" "
               << "format=\"ascii\""
               << ">\n";

      for (unsigned j = 0; j < number_of_elements; j++)
      {
        // Cast to FiniteElement and (in paranoid mode) check
        // if cast has failed.
        FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(j));

#ifdef PARANOID
        if (fe_pt == 0)
        {
          std::stringstream error_stream;
          error_stream << "Recast for element " << j << " failed" << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        fe_pt->scalar_value_fct_paraview(file_out, i, nplot, exact_soln_pt);
      }

      // Close of the DataArray
      file_out << "</DataArray>\n";
    }

    // Close off the PointData set
    file_out << "</PointData>\n";


    // Geometric Points
    //------------------

    file_out << "<Points>\n"
             << "<DataArray type=\"Float32\""
             << " NumberOfComponents=\""
             // This always has to be 3 for an unstructured grid
             << 3 << "\" "
             << "format=\"ascii\">\n";

    // Loop over all the elements to print their plot points
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " faild" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      fe_pt->output_paraview(file_out, nplot);
    }

    file_out << "</DataArray>\n"
             << "</Points>\n";


    // Cells
    //-------

    file_out
      << "<Cells>\n"
      << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    // Make counter for keeping track of all the local elements,
    // because Paraview requires global coordinates
    unsigned counter = 0;

    // Write connectivity with the local elements
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " faild" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      fe_pt->write_paraview_output_offset_information(file_out, nplot, counter);
    }

    file_out << "</DataArray>\n"
             << "<DataArray type=\"Int32\" "
             << "Name=\"offsets\" format=\"ascii\">\n";

    // Make variable that holds the current offset number
    unsigned offset_sum = 0;

    // Write the offset for the specific elements
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " failed" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      fe_pt->write_paraview_offsets(file_out, nplot, offset_sum);
    }

    file_out << "</DataArray>\n"
             << "<DataArray type=\"UInt8\" Name=\"types\">\n";

    // Loop over all elements to get the type that they have
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " failed" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      fe_pt->write_paraview_type(file_out, nplot);
    }

    file_out << "</DataArray>\n"
             << "</Cells>\n";


    // File Closure
    //-------------
    file_out << "</Piece>\n"
             << "</UnstructuredGrid>\n"
             << "</VTKFile>";
  }


  //========================================================
  /// Output in paraview format into specified file.
  ///
  /// Breaks up each element into sub-elements for plotting
  /// purposes. We assume that all elements are of the same
  /// type (fct will break (in paranoid mode) if paraview
  /// output fcts of the elements are inconsistent).
  //========================================================
  void Mesh::output_fct_paraview(
    std::ofstream& file_out,
    const unsigned& nplot,
    const double& time,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt) const
  {
    // Change the scientific format so that E is used rather than e
    file_out.setf(std::ios_base::uppercase);

    // Decide how many elements there are to be plotted
    unsigned long number_of_elements = this->Element_pt.size();

    // Cast to finite element and return if cast fails.
    FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(0));

#ifdef PARANOID
    if (fe_pt == 0)
    {
      throw OomphLibError("Recast for FiniteElement failed for element 0!\n",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif


#ifdef PARANOID
    // Check if all elements have the same number of degrees of freedom,
    // if they don't, paraview will break
    unsigned el_zero_ndof = fe_pt->nscalar_paraview();
    for (unsigned i = 1; i < number_of_elements; i++)
    {
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));
      unsigned el_i_ndof = fe_pt->nscalar_paraview();
      if (el_zero_ndof != el_i_ndof)
      {
        std::stringstream error_stream;
        error_stream
          << "Element " << i << " has different number of degrees of freedom\n"
          << "than from previous elements, Paraview cannot handle this.\n"
          << "We suggest that the problem is broken up into submeshes instead."
          << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Make variables to hold the number of nodes and elements
    unsigned long number_of_nodes = 0;
    unsigned long total_number_of_elements = 0;

    // Loop over all the elements to find total number of plot points
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " failed" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      number_of_nodes += fe_pt->nplot_points_paraview(nplot);
      total_number_of_elements += fe_pt->nsub_elements_paraview(nplot);
    }


    // File Declaration
    //------------------

    // Insert the necessary lines plus header of file, and
    // number of nodes and elements
    file_out << "<?xml version=\"1.0\"?>\n"
             << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
             << "byte_order=\"LittleEndian\">\n"
             << "<UnstructuredGrid>\n"
             << "<Piece NumberOfPoints=\"" << number_of_nodes
             << "\" NumberOfCells=\"" << total_number_of_elements << "\">\n";


    // Point Data
    //-----------

    // Check the number of degrees of freedom
    unsigned ndof = fe_pt->nscalar_paraview();

    // Point data is going in here
    file_out << "<PointData ";

    // Insert just the first scalar name, since paraview reads everything
    // else after that as being of the same type. Get information from
    // first element.
    file_out << "Scalars=\"" << fe_pt->scalar_name_paraview(0) << "\">\n";

    // Loop over i scalar fields and j number of elements
    for (unsigned i = 0; i < ndof; i++)
    {
      file_out << "<DataArray type=\"Float32\" "
               << "Name=\"" << fe_pt->scalar_name_paraview(i) << "\" "
               << "format=\"ascii\""
               << ">\n";

      for (unsigned j = 0; j < number_of_elements; j++)
      {
        // Cast to FiniteElement and (in paranoid mode) check
        // if cast has failed.
        FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(j));

#ifdef PARANOID
        if (fe_pt == 0)
        {
          std::stringstream error_stream;
          error_stream << "Recast for element " << j << " failed" << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        fe_pt->scalar_value_fct_paraview(
          file_out, i, nplot, time, exact_soln_pt);
      }

      // Close of the DataArray
      file_out << "</DataArray>\n";
    }

    // Close off the PointData set
    file_out << "</PointData>\n";


    // Geometric Points
    //------------------

    file_out << "<Points>\n"
             << "<DataArray type=\"Float32\""
             << " NumberOfComponents=\""
             // This always has to be 3 for an unstructured grid
             << 3 << "\" "
             << "format=\"ascii\">\n";

    // Loop over all the elements to print their plot points
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " faild" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      fe_pt->output_paraview(file_out, nplot);
    }

    file_out << "</DataArray>\n"
             << "</Points>\n";


    // Cells
    //-------

    file_out
      << "<Cells>\n"
      << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    // Make counter for keeping track of all the local elements,
    // because Paraview requires global coordinates
    unsigned counter = 0;

    // Write connectivity with the local elements
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " faild" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      fe_pt->write_paraview_output_offset_information(file_out, nplot, counter);
    }

    file_out << "</DataArray>\n"
             << "<DataArray type=\"Int32\" "
             << "Name=\"offsets\" format=\"ascii\">\n";

    // Make variable that holds the current offset number
    unsigned offset_sum = 0;

    // Write the offset for the specific elements
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " failed" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      fe_pt->write_paraview_offsets(file_out, nplot, offset_sum);
    }

    file_out << "</DataArray>\n"
             << "<DataArray type=\"UInt8\" Name=\"types\">\n";

    // Loop over all elements to get the type that they have
    for (unsigned i = 0; i < number_of_elements; i++)
    {
      // Cast to FiniteElement and (in paranoid mode) check
      // if cast has failed.
      FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(element_pt(i));

#ifdef PARANOID
      if (fe_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Recast for element " << i << " failed" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      fe_pt->write_paraview_type(file_out, nplot);
    }

    file_out << "</DataArray>\n"
             << "</Cells>\n";


    // File Closure
    //-------------
    file_out << "</Piece>\n"
             << "</UnstructuredGrid>\n"
             << "</VTKFile>";
  }


  //========================================================
  /// Output function for the mesh class
  ///
  /// Loop over all elements and plot (i.e. execute
  /// the element's own output() function)
  //========================================================
  void Mesh::output(std::ostream& outfile)
  {
    // Loop over the elements and call their output functions
    // Assign Element_pt_range
    unsigned long Element_pt_range = Element_pt.size();
    for (unsigned long e = 0; e < Element_pt_range; e++)
    {
      // Try to cast to FiniteElement
      FiniteElement* el_pt = dynamic_cast<FiniteElement*>(Element_pt[e]);
      if (el_pt == 0)
      {
        oomph_info << "Can't execute output(...) for non FiniteElements"
                   << std::endl;
      }
      else
      {
#ifdef OOMPH_HAS_MPI
        if (Output_halo_elements)
#endif
        {
          el_pt->output(outfile);
        }
#ifdef OOMPH_HAS_MPI
        else
        {
          if (!el_pt->is_halo())
          {
            el_pt->output(outfile);
          }
        }
#endif
      }
    }
  }

  //========================================================
  /// Output function for the mesh class
  ///
  /// Loop over all elements and plot (i.e. execute
  /// the element's own output() function). Use
  /// n_plot plot points in each coordinate direction.
  //========================================================
  void Mesh::output(std::ostream& outfile, const unsigned& n_plot)
  {
    // Loop over the elements and call their output functions
    // Assign Element_pt_range
    unsigned long Element_pt_range = Element_pt.size();

    for (unsigned long e = 0; e < Element_pt_range; e++)
    {
      // Try to cast to FiniteElement
      FiniteElement* el_pt = dynamic_cast<FiniteElement*>(Element_pt[e]);
      if (el_pt == 0)
      {
        oomph_info << "Can't execute output(...) for non FiniteElements"
                   << std::endl;
      }
      else
      {
#ifdef OOMPH_HAS_MPI
        if (Output_halo_elements)
#endif
        {
          el_pt->output(outfile, n_plot);
        }
#ifdef OOMPH_HAS_MPI
        else
        {
          if (!el_pt->is_halo())
          {
            el_pt->output(outfile, n_plot);
          }
        }
#endif
      }
    }
  }


  //========================================================
  /// Output function for the mesh class
  ///
  /// Loop over all elements and plot (i.e. execute
  /// the element's own output() function)
  /// (C style output)
  //========================================================
  void Mesh::output(FILE* file_pt)
  {
    // Loop over the elements and call their output functions
    // Assign Element_pt_range
    unsigned long Element_pt_range = Element_pt.size();
    for (unsigned long e = 0; e < Element_pt_range; e++)
    {
      // Try to cast to FiniteElement
      FiniteElement* el_pt = dynamic_cast<FiniteElement*>(Element_pt[e]);
      if (el_pt == 0)
      {
        oomph_info << "Can't execute output(...) for non FiniteElements"
                   << std::endl;
      }
      else
      {
#ifdef OOMPH_HAS_MPI
        if (Output_halo_elements)
#endif
        {
          el_pt->output(file_pt);
        }
#ifdef OOMPH_HAS_MPI
        else
        {
          if (!el_pt->is_halo())
          {
            el_pt->output(file_pt);
          }
        }
#endif
      }
    }
  }

  //========================================================
  /// Output function for the mesh class
  ///
  /// Loop over all elements and plot (i.e. execute
  /// the element's own output() function). Use
  /// n_plot plot points in each coordinate direction.
  /// (C style output)
  //========================================================
  void Mesh::output(FILE* file_pt, const unsigned& n_plot)
  {
    // Loop over the elements and call their output functions
    // Assign Element_pt_range
    unsigned long Element_pt_range = Element_pt.size();
    for (unsigned long e = 0; e < Element_pt_range; e++)
    {
      // Try to cast to FiniteElement
      FiniteElement* el_pt = dynamic_cast<FiniteElement*>(Element_pt[e]);
      if (el_pt == 0)
      {
        oomph_info << "Can't execute output(...) for non FiniteElements"
                   << std::endl;
      }
      else
      {
#ifdef OOMPH_HAS_MPI
        if (Output_halo_elements)
#endif
        {
          el_pt->output(file_pt, n_plot);
        }
#ifdef OOMPH_HAS_MPI
        else
        {
          if (!el_pt->is_halo())
          {
            el_pt->output(file_pt, n_plot);
          }
        }
#endif
      }
    }
  }


  //========================================================
  /// Output function for the mesh class
  ///
  /// Loop over all elements and plot (i.e. execute
  /// the element's own output() function). Use
  /// n_plot plot points in each coordinate direction.
  //========================================================
  void Mesh::output_fct(std::ostream& outfile,
                        const unsigned& n_plot,
                        FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Loop over the elements and call their output functions
    // Assign Element_pt_range
    unsigned long Element_pt_range = Element_pt.size();
    for (unsigned long e = 0; e < Element_pt_range; e++)
    {
      // Try to cast to FiniteElement
      FiniteElement* el_pt = dynamic_cast<FiniteElement*>(Element_pt[e]);
      if (el_pt == 0)
      {
        oomph_info << "Can't execute output_fct(...) for non FiniteElements"
                   << std::endl;
      }
      else
      {
#ifdef OOMPH_HAS_MPI
        if (Output_halo_elements)
#endif
        {
          el_pt->output_fct(outfile, n_plot, exact_soln_pt);
        }
#ifdef OOMPH_HAS_MPI
        else
        {
          if (!el_pt->is_halo())
          {
            el_pt->output_fct(outfile, n_plot, exact_soln_pt);
          }
        }
#endif
      }
    }
  }

  //========================================================
  /// Output function for the mesh class
  ///
  /// Loop over all elements and plot (i.e. execute
  /// the element's own output() function) at time t. Use
  /// n_plot plot points in each coordinate direction.
  //========================================================
  void Mesh::output_fct(std::ostream& outfile,
                        const unsigned& n_plot,
                        const double& time,
                        FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {
    // Loop over the elements and call their output functions
    // Assign Element_pt_range
    unsigned long Element_pt_range = Element_pt.size();
    for (unsigned long e = 0; e < Element_pt_range; e++)
    {
      // Try to cast to FiniteElement
      FiniteElement* el_pt = dynamic_cast<FiniteElement*>(Element_pt[e]);
      if (el_pt == 0)
      {
        oomph_info << "Can't execute output_fct(...) for non FiniteElements"
                   << std::endl;
      }
      else
      {
#ifdef OOMPH_HAS_MPI
        if (Output_halo_elements)
#endif
        {
          el_pt->output_fct(outfile, n_plot, time, exact_soln_pt);
        }
#ifdef OOMPH_HAS_MPI
        else
        {
          if (!el_pt->is_halo())
          {
            el_pt->output_fct(outfile, n_plot, time, exact_soln_pt);
          }
        }
#endif
      }
    }
  }

  //==================================================================
  /// Assign the initial values for an impulsive start, which is
  /// acheived by looping over all data in the mesh (internal element
  /// data and data stored at nodes) and setting the calling the
  /// assign_initial_values_impulsive() function for each data's
  /// timestepper
  //=================================================================
  void Mesh::assign_initial_values_impulsive()
  {
    // Loop over the elements
    unsigned long Nelement = nelement();
    for (unsigned long e = 0; e < Nelement; e++)
    {
      // Find the number of internal dofs
      unsigned Ninternal = element_pt(e)->ninternal_data();
      // Loop over internal dofs and shift the time values
      // using the internals data's timestepper
      for (unsigned j = 0; j < Ninternal; j++)
      {
        element_pt(e)
          ->internal_data_pt(j)
          ->time_stepper_pt()
          ->assign_initial_values_impulsive(element_pt(e)->internal_data_pt(j));
      }
    }

    // Loop over the nodes
    unsigned long n_node = nnode();
    for (unsigned long n = 0; n < n_node; n++)
    {
      // Assign initial values using the Node's timestepper
      Node_pt[n]->time_stepper_pt()->assign_initial_values_impulsive(
        Node_pt[n]);
      // Assign initial positions using the Node's timestepper
      Node_pt[n]
        ->position_time_stepper_pt()
        ->assign_initial_positions_impulsive(Node_pt[n]);
    }
  }

  //===============================================================
  /// Shift time-dependent data along for next timestep:
  /// Again this is achieved by looping over all data and calling
  /// the functions defined in each data object's timestepper.
  //==============================================================
  void Mesh::shift_time_values()
  {
    // Loop over the elements which shift their internal data
    // via their own timesteppers
    const unsigned long Nelement = nelement();
    for (unsigned long e = 0; e < Nelement; e++)
    {
      // Find the number of internal dofs
      const unsigned Ninternal = element_pt(e)->ninternal_data();
      // Loop over internal dofs and shift the time values
      // using the internals data's timestepper
      for (unsigned j = 0; j < Ninternal; j++)
      {
        element_pt(e)
          ->internal_data_pt(j)
          ->time_stepper_pt()
          ->shift_time_values(element_pt(e)->internal_data_pt(j));
      }
    }

    // Loop over the nodes
    const unsigned long n_node = nnode();
    for (unsigned long n = 0; n < n_node; n++)
    {
      // Shift the Data associated with the nodes with the Node's own
      // timestepper
      Node_pt[n]->time_stepper_pt()->shift_time_values(Node_pt[n]);
      // Push history of nodal positions back
      Node_pt[n]->position_time_stepper_pt()->shift_time_positions(Node_pt[n]);
    }
  }

  //=========================================================================
  /// Calculate predictions for all Data and positions associated
  /// with the mesh. This is usually only used for adaptive time-stepping
  /// when the comparison between a predicted value and the actual value
  /// is usually used to determine the change in step size. Again the
  /// loop is over all data in the mesh and individual timestepper functions
  /// for each data value are called.
  //=========================================================================
  void Mesh::calculate_predictions()
  {
    // Loop over the elements which shift their internal data
    // via their own timesteppers
    const unsigned long Nelement = nelement();
    for (unsigned long e = 0; e < Nelement; e++)
    {
      // Find the number of internal dofs
      const unsigned Ninternal = element_pt(e)->ninternal_data();
      // Loop over internal dofs and calculate predicted positions
      // using the internals data's timestepper
      for (unsigned j = 0; j < Ninternal; j++)
      {
        element_pt(e)
          ->internal_data_pt(j)
          ->time_stepper_pt()
          ->calculate_predicted_values(element_pt(e)->internal_data_pt(j));
      }
    }

    // Loop over the nodes
    const unsigned long n_node = nnode();
    for (unsigned long n = 0; n < n_node; n++)
    {
      // Calculate the predicted values at the nodes
      Node_pt[n]->time_stepper_pt()->calculate_predicted_values(Node_pt[n]);
      // Calculate the predicted positions
      Node_pt[n]->position_time_stepper_pt()->calculate_predicted_positions(
        Node_pt[n]);
    }
  }

  //===============================================================
  /// Virtual function that should be overloaded if the mesh
  /// has any mesh level storage of the timestepper
  //==================================================================
  void Mesh::set_mesh_level_time_stepper(TimeStepper* const& time_stepper_pt,
                                         const bool& preserve_existing_data)
  {
#ifdef PARANOID
    if (!Suppress_warning_about_empty_mesh_level_time_stepper_function)
    {
      std::ostringstream warning_stream;
      warning_stream
        << "Empty set_mesh_level_time_stepper() has been called.\n"
        << "This function needs to be overloaded to reset any (pointers to) \n"
        << "timesteppers for meshes that store timesteppers in locations "
           "other\n"
        << "than the Nodes or Elements;\n"
        << "e.g. SpineMeshes have SpineData with timesteppers,\n"
        << "Triangle and TetMeshes store the timestepper for use in "
           "adaptivity.\n\n\n";
      warning_stream
        << "If you are solving a continuation or bifurcation detecion\n"
        << "problem and strange things are happening, then check that\n"
        << "you don't need to overload this function for your mesh."
        << "\n This warning can be suppressed by setting:\n"
        << "Mesh::Suppress_warning_about_empty_mesh_level_time_stepper_"
           "function=true"
        << std::endl;
      OomphLibWarning(
        warning_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif
  }

  //===============================================================
  /// Set the values of auxilliary data used in continuation problems
  /// when the data is pinned.
  //==================================================================
  void Mesh::set_consistent_pinned_values_for_continuation(
    ContinuationStorageScheme* const& continuation_storage_pt)
  {
    // Loop over the nodes
    const unsigned long n_node = this->nnode();
    for (unsigned long n = 0; n < n_node; n++)
    {
      continuation_storage_pt->set_consistent_pinned_values(this->Node_pt[n]);
      continuation_storage_pt->set_consistent_pinned_positions(
        this->Node_pt[n]);
    }

    // Loop over the elements
    const unsigned long n_element = this->nelement();
    for (unsigned long e = 0; e < n_element; e++)
    {
      // Cache pointer to the elemnet
      GeneralisedElement* const elem_pt = this->element_pt(e);
      // Find the number of internal dofs
      const unsigned n_internal = elem_pt->ninternal_data();

      // Loop over internal dofs and test the data
      for (unsigned j = 0; j < n_internal; j++)
      {
        continuation_storage_pt->set_consistent_pinned_values(
          elem_pt->internal_data_pt(j));
      }
    }
  }


  //===============================================================
  /// Return true if the pointer corresponds to any data stored in
  /// the mesh and false if not
  //==================================================================
  bool Mesh::does_pointer_correspond_to_mesh_data(double* const& parameter_pt)
  {
    // Loop over the nodes
    const unsigned long n_node = this->nnode();
    for (unsigned long n = 0; n < n_node; n++)
    {
      // Check the values and positional data associated with each node
      if ((this->Node_pt[n]->does_pointer_correspond_to_value(parameter_pt)) ||
          (this->Node_pt[n]->does_pointer_correspond_to_position_data(
            parameter_pt)))
      {
        return true;
      }
    }

    // Loop over the elements
    const unsigned long n_element = this->nelement();
    for (unsigned long e = 0; e < n_element; e++)
    {
      // Cache pointer to the elemnet
      GeneralisedElement* const elem_pt = this->element_pt(e);

      // Find the number of internal dofs
      const unsigned n_internal = elem_pt->ninternal_data();

      // Loop over internal dofs and test the data
      for (unsigned j = 0; j < n_internal; j++)
      {
        if (elem_pt->internal_data_pt(j)->does_pointer_correspond_to_value(
              parameter_pt))
        {
          return true;
        }
      }
    }

    // If we get here we haven't found the data, so return false
    return false;
  }


  //===============================================================
  /// Set the time stepper associated with all the nodal data
  /// in the problem
  //==============================================================
  void Mesh::set_nodal_time_stepper(TimeStepper* const& time_stepper_pt,
                                    const bool& preserve_existing_data)
  {
    // Loop over the nodes
    const unsigned long n_node = this->nnode();
    for (unsigned long n = 0; n < n_node; n++)
    {
      // Set the timestepper associated with each node
      this->Node_pt[n]->set_time_stepper(time_stepper_pt,
                                         preserve_existing_data);
      this->Node_pt[n]->set_position_time_stepper(time_stepper_pt,
                                                  preserve_existing_data);
    }
  }

  //===============================================================
  /// Set the time stepper associated with all internal data stored
  /// in the elements in the mesh
  //===============================================================
  void Mesh::set_elemental_internal_time_stepper(
    TimeStepper* const& time_stepper_pt, const bool& preserve_existing_data)
  {
    // Loop over the elements
    const unsigned long n_element = this->nelement();
    for (unsigned long e = 0; e < n_element; e++)
    {
      // Find the number of internal dofs
      const unsigned n_internal = this->element_pt(e)->ninternal_data();

      // Loop over internal dofs and set the timestepper
      for (unsigned j = 0; j < n_internal; j++)
      {
        this->element_pt(e)->internal_data_pt(j)->set_time_stepper(
          time_stepper_pt, preserve_existing_data);
      }
    }
  }

  //========================================================================
  /// A function that upgrades an ordinary node to a boundary node.
  /// All pointers to the node from the mesh's elements are found.
  /// and replaced by pointers to the new boundary node. If the node
  /// is present in the mesh's list of nodes, that pointer is also
  /// replaced. Finally, the pointer argument node_pt addresses the new
  /// node on return from the function.
  /// We shouldn't ever really use this, but it does make life that
  /// bit easier for the lazy mesh writer.
  //=======================================================================
  void Mesh::convert_to_boundary_node(Node*& node_pt)
  {
    // Cache a list of FiniteElement pointers for use in this function.
    Vector<FiniteElement*> fe_pt(nelement(), 0);
    for (unsigned e = 0, ne = nelement(); e < ne; e++)
    {
      // Some elements may not have been build yet, just store a null pointer
      // for these cases.
      if (Element_pt[e] == 0) fe_pt[e] = 0;
      else
        fe_pt[e] = finite_element_pt(e);
    }

    // Now call the real function
    convert_to_boundary_node(node_pt, fe_pt);
  }

  // ============================================================
  /// As convert_to_boundary_node but with a vector of pre-"dynamic cast"ed
  /// pointers passed in. If this function is being called often then
  /// creating this vector and passing it in explicitly each time can give a
  /// large speed up.
  // Note: the real reason that this function is so slow in the first place
  // is because it has to loop over all elements. So if you use this function
  // O(N) times your boundary node creation complexity is O(N^2).
  // ============================================================
  void Mesh::convert_to_boundary_node(
    Node*& node_pt, const Vector<FiniteElement*>& finite_element_pt)
  {
    // If the node is already a boundary node, then return straight away,
    // we don't need to do anything
    if (dynamic_cast<BoundaryNodeBase*>(node_pt) != 0)
    {
      return;
    }

    // Loop over all the elements in the mesh and find all those in which
    // the present node is referenced and the corresponding local node number
    // in those elements.

    // Storage for elements and local node number
    std::list<std::pair<unsigned long, int>>
      list_of_elements_and_local_node_numbers;

    // Loop over all elements
    unsigned long n_element = this->nelement();
    for (unsigned long e = 0; e < n_element; e++)
    {
      // Buffer the case when we have not yet filled up the element array
      // Unfortunately, we should not assume that the array has been filled
      // in a linear order, so we can't break out early.
      if (Element_pt[e] != 0)
      {
        // Find the local node number of the passed node
        int node_number = finite_element_pt[e]->get_node_number(node_pt);
        // If the node is present in the element, add it to our list and
        // NULL out the local element entries
        if (node_number != -1)
        {
          list_of_elements_and_local_node_numbers.insert(
            list_of_elements_and_local_node_numbers.end(),
            std::make_pair(e, node_number));
          // Null it out
          finite_element_pt[e]->node_pt(node_number) = 0;
        }
      }
    } // End of loop over elements

    // If there are no entries in the list we are in real trouble
    if (list_of_elements_and_local_node_numbers.empty())
    {
      std::ostringstream error_stream;
      error_stream << "Node " << node_pt
                   << " is not contained in any elements in the Mesh."
                   << std::endl
                   << "How was it created then?" << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


    // Create temporary storage for a pointer to the old node.
    // This is required because if we have passed a reference to the
    // first element that we find, constructing the new node
    // will over-write our pointer and we'll get segmentation faults.
    Node* old_node_pt = node_pt;

    // We now create the new node by using the first element in the list
    std::list<std::pair<unsigned long, int>>::iterator list_it =
      list_of_elements_and_local_node_numbers.begin();

    // Create a new boundary node, using the timestepper from the
    // original node
    Node* new_node_pt =
      finite_element_pt[list_it->first]->construct_boundary_node(
        list_it->second, node_pt->time_stepper_pt());

    // Now copy all the information accross from the old node

    // Can we cast the node to a solid node
    SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(new_node_pt);
    // If it's a solid node, do the casting
    if (solid_node_pt != 0)
    {
      solid_node_pt->copy(dynamic_cast<SolidNode*>(old_node_pt));
    }
    else
    {
      new_node_pt->copy(old_node_pt);
    }

    // Loop over all other elements in the list and set their pointers
    // to the new node
    for (++list_it; // Increment the iterator
         list_it != list_of_elements_and_local_node_numbers.end();
         ++list_it)
    {
      finite_element_pt[list_it->first]->node_pt(list_it->second) = new_node_pt;
    }

    // Finally, find the position of the node in the global mesh
    Vector<Node*>::iterator it =
      std::find(Node_pt.begin(), Node_pt.end(), old_node_pt);

    // If it is in the mesh, update the pointer
    if (it != Node_pt.end())
    {
      *it = new_node_pt;
    }

    // Can now delete the old node
    delete old_node_pt;

    // Replace the passed pointer by a pointer to the new node
    // Note that in most cases, this will be wasted work because the node
    // pointer will either the pointer in the mesh's or an element's
    // node_pt vector. Still assignment is quicker than an if to check this.
    node_pt = new_node_pt;
  }

#ifdef OOMPH_HAS_MPI

  //========================================================================
  /// Setup shared node scheme: Shared node lookup scheme contains
  /// a unique correspondence between all nodes on the halo/haloed
  /// elements between two processors.
  //========================================================================
  void Mesh::setup_shared_node_scheme()
  {
    // Determine the shared nodes lookup scheme - all nodes located on the
    // halo(ed) elements between two domains.  This scheme is necessary in order
    // to identify master nodes that may not be present in the halo-haloed
    // element lookup scheme between two processors (for example, if the node
    // is on an element which is in a lookup scheme between two higher-numbered
    // processors)

    double t_start = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start = TimingHelpers::timer();
    }

    // Store number of processors and current process
    int n_proc = Comm_pt->nproc();
    int my_rank = Comm_pt->my_rank();

    // Need to clear the shared node scheme first
    Shared_node_pt.clear();

    for (int d = 0; d < n_proc; d++)
    {
      // map of bools for whether the node has been shared,
      // initialised to 0 (false) for each domain d
      std::map<Node*, bool> node_shared;

      // For all domains lower than the current domain: Do halos first
      // then haloed, to ensure correct order in lookup scheme from
      // the other side
      if (d < my_rank)
      {
        // Get the nodes from the halo elements first
        Vector<GeneralisedElement*> halo_elem_pt(this->halo_element_pt(d));
        unsigned nhalo_elem = halo_elem_pt.size();

        for (unsigned e = 0; e < nhalo_elem; e++)
        {
          // Get finite element
          FiniteElement* el_pt = dynamic_cast<FiniteElement*>(halo_elem_pt[e]);
          if (el_pt != 0)
          {
            // Loop over nodes
            unsigned nnod = el_pt->nnode();
            for (unsigned j = 0; j < nnod; j++)
            {
              Node* nod_pt = el_pt->node_pt(j);

              // Add it as a shared node from current domain
              if (!node_shared[nod_pt])
              {
                this->add_shared_node_pt(d, nod_pt);
                node_shared[nod_pt] = true;
              }

            } // end loop over nodes
          }

        } // end loop over elements

        // Now get any left over nodes on the haloed elements
        Vector<GeneralisedElement*> haloed_elem_pt(this->haloed_element_pt(d));
        unsigned nhaloed_elem = haloed_elem_pt.size();

        for (unsigned e = 0; e < nhaloed_elem; e++)
        {
          // Get element
          FiniteElement* el_pt =
            dynamic_cast<FiniteElement*>(haloed_elem_pt[e]);
          if (el_pt != 0)
          {
            // Loop over the nodes
            unsigned nnod = el_pt->nnode();
            for (unsigned j = 0; j < nnod; j++)
            {
              Node* nod_pt = el_pt->node_pt(j);

              // Add it as a shared node from current domain
              if (!node_shared[nod_pt])
              {
                this->add_shared_node_pt(d, nod_pt);
                node_shared[nod_pt] = true;
              }

            } // end loop over nodes
          }
        } // end loop over elements
      }

      // If the domain is bigger than the current rank: Do haloed first
      // then halo, to ensure correct order in lookup scheme from
      // the other side
      if (d > my_rank)
      {
        // Get the nodes from the haloed elements first
        Vector<GeneralisedElement*> haloed_elem_pt(this->haloed_element_pt(d));
        unsigned nhaloed_elem = haloed_elem_pt.size();

        for (unsigned e = 0; e < nhaloed_elem; e++)
        {
          // Get element
          FiniteElement* el_pt =
            dynamic_cast<FiniteElement*>(haloed_elem_pt[e]);
          if (el_pt != 0)
          {
            // Loop over nodes
            unsigned nnod = el_pt->nnode();
            for (unsigned j = 0; j < nnod; j++)
            {
              Node* nod_pt = el_pt->node_pt(j);

              // Add it as a shared node from current domain
              if (!node_shared[nod_pt])
              {
                this->add_shared_node_pt(d, nod_pt);
                node_shared[nod_pt] = true;
              }

            } // end loop over nodes
          }
        } // end loop over elements

        // Now get the nodes from any halo elements left over
        Vector<GeneralisedElement*> halo_elem_pt(this->halo_element_pt(d));
        unsigned nhalo_elem = halo_elem_pt.size();

        for (unsigned e = 0; e < nhalo_elem; e++)
        {
          // Get element
          FiniteElement* el_pt = dynamic_cast<FiniteElement*>(halo_elem_pt[e]);
          if (el_pt != 0)
          {
            // Loop over nodes
            unsigned nnod = el_pt->nnode();
            for (unsigned j = 0; j < nnod; j++)
            {
              Node* nod_pt = el_pt->node_pt(j);

              // Add it as a shared node from current domain
              if (!node_shared[nod_pt])
              {
                this->add_shared_node_pt(d, nod_pt);
                node_shared[nod_pt] = true;
              }

            } // end loop over nodes
          }
        } // end loop over elements

      } // end if (d ...)

    } // end loop over processes


    double t_end = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for identification of shared nodes: "
                 << t_end - t_start << std::endl;
      oomph_info.stream_pt()->flush();
    }
  }

  //========================================================================
  /// Synchronise shared node lookup schemes to cater for the
  /// the case where:
  /// (1) a certain node on the current processor is halo with proc p
  ///     (i.e. its non-halo counterpart lives on processor p)
  /// (2) that node also exists (also as a halo) on another processor
  ///     (q, say) where its non-halo counter part is also known to be
  ///     on processor p.
  /// However, without calling this function the current processor does not
  /// necessarily know that it shares a node with processor q. This
  /// information can be required, e.g. when synchronising hanging node
  /// schemes over all processors.
  //========================================================================
  void Mesh::synchronise_shared_nodes(const bool& report_stats)
  {
    double t_start = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start = TimingHelpers::timer();
    }

    double tt_start = 0.0;
    double tt_end = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      tt_start = TimingHelpers::timer();
    }

    // Storage for current processor and number of processors
    int n_proc = Comm_pt->nproc();
    int my_rank = Comm_pt->my_rank();


#ifdef PARANOID
    // Has some twit filled in shared nodes with own process?!
    // Check at start of function
    if (Shared_node_pt[my_rank].size() != 0)
    {
      throw OomphLibError(
        "Processor has shared nodes with itself! Something's gone wrong!",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // Stage 1: Populate the set of of processor IDs that
    // each haloed node on current processor is haloed by.
    std::map<Node*, std::set<int>> shared_domain_set;

    // Associate unique number with any haloed nodes on this processor
    std::map<Node*, unsigned> global_haloed_node_number_plus_one;
    unsigned global_haloed_count = 0;

    // Loop over domains
    for (int d = 0; d < n_proc; d++)
    {
      // Don't talk to yourself
      if (d != my_rank)
      {
        // Loop over haloed nodes
        unsigned nnod_haloed = this->nhaloed_node(d);
        for (unsigned j = 0; j < nnod_haloed; j++)
        {
          Node* nod_pt = this->haloed_node_pt(d, j);
          shared_domain_set[nod_pt].insert(d);
          if (global_haloed_node_number_plus_one[nod_pt] == 0)
          {
            global_haloed_node_number_plus_one[nod_pt] =
              global_haloed_count + 1;
            global_haloed_count++;
          }
        }
      }
    }

    if (Global_timings::Doc_comprehensive_timings)
    {
      tt_end = TimingHelpers::timer();
      oomph_info << "Time for initial classification in "
                    "Mesh::synchronise_shared_nodes(): "
                 << tt_end - tt_start << std::endl;
      tt_start = TimingHelpers::timer();
    }


    // Stage 2: All-to-all communication to inform all processors that
    // hold halo nodes with current processor about all domains that the current
    // processor shares nodes with [This will allow the other processors to add
    // these nodes to their shared node lookup schemes with processors
    // that they currently "can't see"].

    // Data to be sent to each processor
    Vector<int> send_n(n_proc, 0);

    // Storage for all values to be sent to all processors
    Vector<unsigned> send_data;

    // Start location within send_data for data to be sent to each processor
    Vector<int> send_displacement(n_proc, 0);


    // Loop over haloed nodes with other domains
    for (int domain = 0; domain < n_proc; domain++)
    {
      // Set the offset for the current processor
      send_displacement[domain] = send_data.size();

      // Every processor works on haloed nodes with proc "domain" and
      // sends associations across to that domain. No need to talk to
      // yourself...
      if (domain != my_rank)
      {
        // Send total number of global haloed nodes
        // send_data.push_back(global_haloed_count);

        // Loop over haloed nodes
        unsigned nnod_haloed = this->nhaloed_node(domain);
        for (unsigned j = 0; j < nnod_haloed; j++)
        {
          Node* nod_pt = this->haloed_node_pt(domain, j);

          // Send global ID of haloed node
          send_data.push_back(global_haloed_node_number_plus_one[nod_pt] - 1);

          // Get the set of domains that halo this node
          std::set<int> tmp_shared_domain_set = shared_domain_set[nod_pt];

          // Insert number of shared domains into send data
          unsigned n_shared_domains = tmp_shared_domain_set.size();
          send_data.push_back(n_shared_domains);

          // Add shared domains
          for (std::set<int>::iterator it = tmp_shared_domain_set.begin();
               it != tmp_shared_domain_set.end();
               it++)
          {
            send_data.push_back(*it);
          }
        }
      }

      // Find the number of data added to the vector
      send_n[domain] = send_data.size() - send_displacement[domain];
    }


    // Storage for the number of data to be received from each processor
    Vector<int> receive_n(n_proc, 0);

    // Now send numbers of data to be sent between all processors
    MPI_Alltoall(
      &send_n[0], 1, MPI_INT, &receive_n[0], 1, MPI_INT, Comm_pt->mpi_comm());


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
                  Comm_pt->mpi_comm());


    if (Global_timings::Doc_comprehensive_timings)
    {
      tt_end = TimingHelpers::timer();
      oomph_info << "Time for alltoall in Mesh::synchronise_shared_nodes(): "
                 << tt_end - tt_start << std::endl;
      tt_start = TimingHelpers::timer();
    }


    if (Global_timings::Doc_comprehensive_timings)
    {
      oomph_info << "Starting vector to set conversion in "
                 << "Mesh::synchronise_shared_nodes() for a total of "
                 << nshared_node() << " nodes\n";
      tt_start = TimingHelpers::timer();
    }


    // Copy vector-based representation of shared nodes into
    // sets for faster search
    Vector<std::set<Node*>> shared_node_set(n_proc);
    for (int d = 0; d < n_proc; d++)
    {
      unsigned n_vector = Shared_node_pt[d].size();
      for (unsigned i = 0; i < n_vector; i++)
      {
        shared_node_set[d].insert(Shared_node_pt[d][i]);
      }
    }


    if (Global_timings::Doc_comprehensive_timings)
    {
      tt_end = TimingHelpers::timer();
      oomph_info
        << "Time for vector to set in Mesh::synchronise_shared_nodes(): "
        << tt_end - tt_start << std::endl;
      tt_start = TimingHelpers::timer();
    }


    // Now use the received data
    for (int send_rank = 0; send_rank < n_proc; send_rank++)
    {
      // Don't bother to do anything for the processor corresponding to the
      // current processor or if no data were received from this processor
      if ((send_rank != my_rank) && (receive_n[send_rank] != 0))
      {
        // Counter for the data within the large array
        unsigned count = receive_displacement[send_rank];

        // Read total number of global haloed nodes
        // unsigned n_global_haloed_nodes_on_send_proc=receive_data[count++];

        // Storage for nodes and associated domains:
        // domains_map[global_haloed_node_number].first = node
        // domains_map[global_haloed_node_number].second = set of domains
        //                                                 this node is
        //                                                 associated with.
        std::map<unsigned, std::pair<Node*, std::set<unsigned>>> domains_map;

        // Loop over halo nodes with sending processor
        unsigned nnod_halo = this->nhalo_node(send_rank);
        for (unsigned j = 0; j < nnod_halo; j++)
        {
          Node* nod_pt = this->halo_node_pt(send_rank, j);

          // Read unique ID of haloed node on send proc
          unsigned haloed_node_id_on_send_proc = receive_data[count++];

          // Read out number of shared domains into send data
          unsigned n_shared_domains = receive_data[count++];

          // Prepare set of domains
          std::set<unsigned> domain_set;

          // Read 'em
          for (unsigned i = 0; i < n_shared_domains; i++)
          {
            int shared_domain = receive_data[count++];

            // Record in set
            domain_set.insert(shared_domain);
          }

          // Add entry:
          domains_map[haloed_node_id_on_send_proc] =
            std::make_pair(nod_pt, domain_set);

        } // end of loop over halo nodes


        // Now add new shared nodes in order
#ifdef PARANOID
        int previous_one = -1;
#endif
        for (std::map<unsigned, std::pair<Node*, std::set<unsigned>>>::iterator
               it = domains_map.begin();
             it != domains_map.end();
             it++)
        {
          // Super-paranoid: Check that the map actually sorted entries
          // by key (as it should)
#ifdef PARANOID
          if (int((*it).first) < previous_one)
          {
            std::ostringstream error_stream;
            error_stream << "Map did not store entries in order of key\n "
                         << "Current key:  " << (*it).first
                         << "; previous one: " << previous_one << "\n"
                         << "Need to rewrite this...\n";
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
          previous_one = (*it).first;
#endif

          // Extract node
          Node* nod_pt = (*it).second.first;

          // Extract set of domains
          std::set<unsigned> domain_set((*it).second.second);

          // Read 'em
          for (std::set<unsigned>::iterator itt = domain_set.begin();
               itt != domain_set.end();
               itt++)
          {
            int shared_domain = (*itt);

            // No need to add shared nodes with oneself!
            if (shared_domain != my_rank)
            {
              // Is node already listed in shared node scheme? Find it
              // and get iterator to entry
              std::set<Node*>::iterator ittt =
                shared_node_set[shared_domain].find(nod_pt);

              // If it's not in there already iterator points to end of
              // set
              if (ittt == shared_node_set[shared_domain].end())
              {
                // Now add it
                add_shared_node_pt(shared_domain, nod_pt);

                // Update set
                shared_node_set[shared_domain].insert(nod_pt);
              }
            }
          }
        }

      } // end of any data is received and ignore own domain

    } // end of loop over send ranks


#ifdef PARANOID
    // Has some twit filled in shared nodes with own process?!
    // Check at end pf function.
    if (Shared_node_pt[my_rank].size() != 0)
    {
      throw OomphLibError(
        "Processor has shared nodes with itself! Something's gone wrong!",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif


    if (Global_timings::Doc_comprehensive_timings)
    {
      tt_end = TimingHelpers::timer();
      oomph_info
        << "Time for final processing in Mesh::synchronise_shared_nodes(): "
        << tt_end - tt_start << std::endl;
      tt_start = TimingHelpers::timer();
    }

    if (Global_timings::Doc_comprehensive_timings)
    {
      double t_end = TimingHelpers::timer();
      oomph_info << "Total time for Mesh::synchronise_shared_nodes(): "
                 << t_end - t_start << std::endl;
    }
  }


  //========================================================================
  /// Classify all halo and haloed information in the mesh
  //========================================================================
  void Mesh::classify_halo_and_haloed_nodes(DocInfo& doc_info,
                                            const bool& report_stats)
  {
    //  MemoryUsage::doc_memory_usage(
    //   "at beginning of Mesh::classify_halo_and_haloed_nodes");

    double t_start = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start = TimingHelpers::timer();
    }

    double tt_start = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      tt_start = TimingHelpers::timer();
    }

    // Set up shared nodes scheme
    setup_shared_node_scheme();

    // MemoryUsage::doc_memory_usage("after setup shared node scheme");

    double tt_end = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      tt_end = TimingHelpers::timer();
      oomph_info << "Time for Mesh::setup_shared_node_scheme() "
                 << " Mesh::classify_halo_and_haloed_nodes(): "
                 << tt_end - tt_start << std::endl;
      oomph_info.stream_pt()->flush();
      tt_start = TimingHelpers::timer();
    }


    // Wipe existing storage schemes for halo(ed) nodes
    Halo_node_pt.clear();
    Haloed_node_pt.clear();

    // Storage for number of processors and current processor
    int n_proc = Comm_pt->nproc();
    int my_rank = Comm_pt->my_rank();
    MPI_Status status;

    // Determine which processors the nodes are associated with
    // and hence who's in charge
    std::map<Data*, std::set<unsigned>> processors_associated_with_data;
    std::map<Data*, unsigned> processor_in_charge;

    // Loop over all processors and associate any nodes on the halo
    // elements involved with that processor
    for (int domain = 0; domain < n_proc; domain++)
    {
      // Get vector of halo elements by copy operation
      Vector<GeneralisedElement*> halo_elem_pt(this->halo_element_pt(domain));

      // Loop over halo elements associated with this adjacent domain
      unsigned nelem = halo_elem_pt.size();
      for (unsigned e = 0; e < nelem; e++)
      {
        // Get element only have nodes if a finite element
        FiniteElement* finite_el_pt =
          dynamic_cast<FiniteElement*>(halo_elem_pt[e]);
        if (finite_el_pt != 0)
        {
          // Loop over nodes
          unsigned nnod = finite_el_pt->nnode();
          for (unsigned j = 0; j < nnod; j++)
          {
            Node* nod_pt = finite_el_pt->node_pt(j);
            // Associate node with this domain
            processors_associated_with_data[nod_pt].insert(domain);

            // Do the same if the node is solid
            SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
            if (solid_nod_pt != 0)
            {
              processors_associated_with_data[solid_nod_pt
                                                ->variable_position_pt()]
                .insert(domain);
            }
          }
        }
      }
    }


    // Loop over all [non-halo] elements and associate their nodes
    // with current procesor
    unsigned nelem = this->nelement();
    for (unsigned e = 0; e < nelem; e++)
    {
      FiniteElement* finite_el_pt =
        dynamic_cast<FiniteElement*>(this->element_pt(e));

      // Only visit non-halos and finite elements
      if ((finite_el_pt != 0) && (!finite_el_pt->is_halo()))
      {
        // Loop over nodes
        unsigned nnod = finite_el_pt->nnode();
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = finite_el_pt->node_pt(j);

          // Associate this node with current processor
          processors_associated_with_data[nod_pt].insert(my_rank);

          // do the same if we have a SolidNode
          SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
          if (solid_nod_pt != 0)
          {
            processors_associated_with_data[solid_nod_pt
                                              ->variable_position_pt()]
              .insert(my_rank);
          }
        }
      }
    }


    if (Global_timings::Doc_comprehensive_timings)
    {
      tt_end = TimingHelpers::timer();
      oomph_info
        << "Time for setup loops in Mesh::classify_halo_and_haloed_nodes: "
        << tt_end - tt_start << std::endl;
      oomph_info.stream_pt()->flush();
      tt_start = TimingHelpers::timer();
    }

    // MemoryUsage::doc_memory_usage("after setup loops");

    // At this point we need to "synchronise" the nodes on halo(ed) elements
    // so that the processors_associated_with_data agrees for the same node
    // on all processors. Strategy: All local nodes have just had their
    // association recorded. Now loop over all haloed elements
    // and send the association of their nodes to the corresponding
    // halo processors where they update/augment the association of the
    // nodes of the corresponding halo elements.

    // Loop over all domains
    for (int d = 0; d < n_proc; d++)
    {
      // Prepare vector to send/receive
      Vector<unsigned> processors_associated_with_data_on_other_proc;

      if (d != my_rank)
      {
        // Communicate the processors associated with data on haloed elements

        // Get haloed elements
        Vector<GeneralisedElement*> haloed_elem_pt(this->haloed_element_pt(d));

        // Initialise counter for this haloed layer
        unsigned count_data = 0;

        // Loop over haloed elements
        unsigned n_haloed_elem = haloed_elem_pt.size();
        for (unsigned e = 0; e < n_haloed_elem; e++)
        {
          // Only nodes in finite elements
          FiniteElement* haloed_el_pt =
            dynamic_cast<FiniteElement*>(haloed_elem_pt[e]);
          if (haloed_el_pt != 0)
          {
            // Loop over nodes
            unsigned n_node = haloed_el_pt->nnode();
            for (unsigned j = 0; j < n_node; j++)
            {
              Node* nod_pt = haloed_el_pt->node_pt(j);

              // Number of processors associated with this node
              unsigned n_assoc = processors_associated_with_data[nod_pt].size();

              // This number needs to be sent
              processors_associated_with_data_on_other_proc.push_back(n_assoc);
              count_data++;

              // Now add the process IDs associated to the vector to be sent
              std::set<unsigned> procs_set =
                processors_associated_with_data[nod_pt];
              for (std::set<unsigned>::iterator it = procs_set.begin();
                   it != procs_set.end();
                   it++)
              {
                processors_associated_with_data_on_other_proc.push_back(*it);
                count_data++;
              }
            }
          }
        }


        // Send the information
        MPI_Send(&count_data, 1, MPI_UNSIGNED, d, 0, Comm_pt->mpi_comm());
        if (count_data != 0)
        {
          MPI_Send(&processors_associated_with_data_on_other_proc[0],
                   count_data,
                   MPI_UNSIGNED,
                   d,
                   1,
                   Comm_pt->mpi_comm());
        }
      }
      else
      {
        // Receive the processors associated with data onto halo elements
        for (int dd = 0; dd < n_proc; dd++)
        {
          if (dd != my_rank) // (my_rank=d)
          {
            // We will be looping over the halo elements with process dd
            Vector<GeneralisedElement*> halo_elem_pt(this->halo_element_pt(dd));
            unsigned n_halo_elem = halo_elem_pt.size();
            unsigned count_data = 0;
            MPI_Recv(&count_data,
                     1,
                     MPI_UNSIGNED,
                     dd,
                     0,
                     Comm_pt->mpi_comm(),
                     &status);
            if (count_data != 0)
            {
              processors_associated_with_data_on_other_proc.resize(count_data);
              MPI_Recv(&processors_associated_with_data_on_other_proc[0],
                       count_data,
                       MPI_UNSIGNED,
                       dd,
                       1,
                       Comm_pt->mpi_comm(),
                       &status);

              // Reset counter and loop through nodes on halo elements
              count_data = 0;
              for (unsigned e = 0; e < n_halo_elem; e++)
              {
                FiniteElement* halo_el_pt =
                  dynamic_cast<FiniteElement*>(halo_elem_pt[e]);
                if (halo_el_pt != 0)
                {
                  unsigned n_node = halo_el_pt->nnode();
                  for (unsigned j = 0; j < n_node; j++)
                  {
                    Node* nod_pt = halo_el_pt->node_pt(j);

                    // Get number of processors associated with data that was
                    // sent
                    unsigned n_assoc =
                      processors_associated_with_data_on_other_proc[count_data];
                    count_data++;

                    for (unsigned i_assoc = 0; i_assoc < n_assoc; i_assoc++)
                    {
                      // Get the process ID
                      unsigned sent_domain =
                        processors_associated_with_data_on_other_proc
                          [count_data];
                      count_data++;

                      // Add it to this processor's list of IDs
                      processors_associated_with_data[nod_pt].insert(
                        sent_domain);

                      // If the node is solid then add the ID to the solid data
                      SolidNode* solid_nod_pt =
                        dynamic_cast<SolidNode*>(nod_pt);
                      if (solid_nod_pt != 0)
                      {
                        processors_associated_with_data
                          [solid_nod_pt->variable_position_pt()]
                            .insert(sent_domain);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (Global_timings::Doc_comprehensive_timings)
    {
      tt_end = TimingHelpers::timer();
      oomph_info
        << "Time for pt2pt send/recv in Mesh::classify_halo_and_haloed_nodes: "
        << tt_end - tt_start << std::endl;
      oomph_info.stream_pt()->flush();
      tt_start = TimingHelpers::timer();
    }


    // MemoryUsage::doc_memory_usage("after pt2pt send/recv");

    // Loop over all nodes on the present processor and put the highest-numbered
    // processor associated with each node "in charge" of the node
    unsigned nnod = this->nnode();
    for (unsigned j = 0; j < nnod; j++)
    {
      Node* nod_pt = this->node_pt(j);

      // Reset halo status of node to false
      nod_pt->set_nonhalo();

      // If it's a SolidNode then the halo status of the data
      // associated with that must also be reset to false
      SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
      if (solid_nod_pt != 0)
      {
        solid_nod_pt->variable_position_pt()->set_nonhalo();
      }

      // Now put the highest-numbered one in charge
      unsigned proc_max = 0;
      std::set<unsigned> procs_set = processors_associated_with_data[nod_pt];
      for (std::set<unsigned>::iterator it = procs_set.begin();
           it != procs_set.end();
           it++)
      {
        if (*it > proc_max) proc_max = *it;
      }
      processor_in_charge[nod_pt] = proc_max;

      // Do the same if we have a SolidNode
      if (solid_nod_pt != 0)
      {
        // Now put the highest-numbered one in charge
        unsigned proc_max_solid = 0;
        std::set<unsigned> procs_set_solid =
          processors_associated_with_data[solid_nod_pt->variable_position_pt()];
        for (std::set<unsigned>::iterator it = procs_set_solid.begin();
             it != procs_set_solid.end();
             it++)
        {
          if (*it > proc_max_solid) proc_max_solid = *it;
        }
        processor_in_charge[solid_nod_pt->variable_position_pt()] =
          proc_max_solid;
      }
    }


    // First stab at determining halo nodes. They are located on the halo
    // elements and the processor in charge differs from the
    // current processor

    // Only count nodes once (map is initialised to 0 = false)
    std::map<Node*, bool> done;

    // Loop over all processors
    for (int domain = 0; domain < n_proc; domain++)
    {
      // Get vector of halo elements by copy operation
      Vector<GeneralisedElement*> halo_elem_pt(this->halo_element_pt(domain));

      // Loop over halo elements associated with this adjacent domain
      unsigned nelem = halo_elem_pt.size();

      for (unsigned e = 0; e < nelem; e++)
      {
        // Get element
        GeneralisedElement* el_pt = halo_elem_pt[e];

        // Can if be cast to a finite element
        FiniteElement* finite_el_pt = dynamic_cast<FiniteElement*>(el_pt);
        if (finite_el_pt != 0)
        {
          // Loop over nodes
          unsigned nnod = finite_el_pt->nnode();
          for (unsigned j = 0; j < nnod; j++)
          {
            Node* nod_pt = finite_el_pt->node_pt(j);

            // Have we done this node already?
            if (!done[nod_pt])
            {
              // Is the other processor/domain in charge of this node?
              int proc_in_charge = processor_in_charge[nod_pt];

              if (proc_in_charge != my_rank)
              {
                // To keep the order of the nodes consistent with that
                // in the haloed node lookup scheme, only
                // allow it to be added when the current domain is in charge
                if (proc_in_charge == int(domain))
                {
                  // Add it as being halo node whose non-halo counterpart
                  // is located on processor proc_in_charge
                  this->add_halo_node_pt(proc_in_charge, nod_pt);

                  // The node itself needs to know it is a halo
                  nod_pt->set_halo(proc_in_charge);

                  // If it's a SolidNode then the data associated with that
                  // must also be halo
                  SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
                  if (solid_nod_pt != 0)
                  {
                    solid_nod_pt->variable_position_pt()->set_halo(
                      proc_in_charge);
                  }

                  // We're done with this node
                  done[nod_pt] = true;
                }
              }
            }
          }

        } // End of finite element case

        // Now make sure internal data on halo elements is also halo
        unsigned nintern_data = el_pt->ninternal_data();
        for (unsigned iintern = 0; iintern < nintern_data; iintern++)
        {
          el_pt->internal_data_pt(iintern)->set_halo(domain);
        }
      }
    }


    // First stab at determining haloed nodes. They are located on the haloed
    // elements and the processor in charge is the current processor

    // Loop over processors that share haloes with the current one
    for (int domain = 0; domain < n_proc; domain++)
    {
      // Only count nodes once (map is initialised to 0 = false)
      std::map<Node*, bool> node_done;

      // Get vector of haloed elements by copy operation
      Vector<GeneralisedElement*> haloed_elem_pt(
        this->haloed_element_pt(domain));

      // Loop over haloed elements associated with this adjacent domain
      unsigned nelem = haloed_elem_pt.size();

      for (unsigned e = 0; e < nelem; e++)
      {
        // Get element
        GeneralisedElement* el_pt = haloed_elem_pt[e];

        // Can it be cast to a finite element
        FiniteElement* finite_el_pt = dynamic_cast<FiniteElement*>(el_pt);
        if (finite_el_pt != 0)
        {
          // Loop over nodes
          unsigned nnod = finite_el_pt->nnode();
          for (unsigned j = 0; j < nnod; j++)
          {
            Node* nod_pt = finite_el_pt->node_pt(j);

            // Have we done this node already?
            if (!node_done[nod_pt])
            {
              // Is the current processor/domain in charge of this node?
              int proc_in_charge = processor_in_charge[nod_pt];

              if (proc_in_charge == my_rank)
              {
                // Add it as being haloed from specified domain
                this->add_haloed_node_pt(domain, nod_pt);

                // We're done with this node
                node_done[nod_pt] = true;
              }
            }
          }
        }
      }
    }


    if (Global_timings::Doc_comprehensive_timings)
    {
      tt_end = TimingHelpers::timer();
      oomph_info
        << "Time for first classific in Mesh::classify_halo_and_haloed_nodes: "
        << tt_end - tt_start << std::endl;
      oomph_info.stream_pt()->flush();
      tt_start = TimingHelpers::timer();
    }

    // MemoryUsage::doc_memory_usage("after first classific");


    // Find any overlooked halo nodes: These are any nodes on the halo/haloed
    // elements (i.e. precisely the nodes currently contained in the shared
    // node scheme) that have not been classified as haloes (yet) though they
    // should have been because another processor is in charge of them.
    // This arises when the "overlooked halo node" is not part of the
    // halo/haloed element lookup scheme between the current processor
    // and the processor that holds the non-halo counterpart. E.g. we're
    // on proc 3. A node at the very edge of its halo layer also exists
    // on procs 0 and 1 with 1 being "in charge". However, the node in
    // question is not part of the halo/haloed element lookup scheme between
    // processor 1 and 3 so in the classification performed above, we never
    // visit it so it's overlooked. The code below rectifies this by going
    // through the intermediate processor (here proc 0) that contains the node
    // in lookup schemes with the halo processor (here proc 3, this one) and the
    // one that contains the non-halo counterpart (here proc 1).


    // Counter for number of overlooked halos (if there aren't any we don't
    // need any comms below)
    unsigned n_overlooked_halo = 0;

    // Record previously overlooked halo nodes so they can be
    // added to the shared node lookup scheme (in a consistent order) below
    Vector<Vector<Node*>> over_looked_halo_node_pt(n_proc);

    // Record previously overlooked haloed nodes so they can be
    // added to the shared node lookup scheme (in a consistent order) below
    Vector<Vector<Node*>> over_looked_haloed_node_pt(n_proc);

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
        unsigned nnod = nshared_node(domain);
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = shared_node_pt(domain, j);

          // Is a different-numbered processor/domain in charge of this node?
          int proc_in_charge = processor_in_charge[nod_pt];
          if ((proc_in_charge != my_rank) && !(nod_pt->is_halo()))
          {
            // Add it as being halo node whose non-halo counterpart
            // is located on processor proc_in_charge
            this->add_halo_node_pt(proc_in_charge, nod_pt);

            // We have another one...
            n_overlooked_halo++;
            over_looked_halo_node_pt[proc_in_charge].push_back(nod_pt);

            // The node itself needs to know it is a halo
            nod_pt->set_halo(proc_in_charge);

            // If it's a SolidNode then the data associated with that
            // must also be halo
            SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
            if (solid_nod_pt != 0)
            {
              solid_nod_pt->variable_position_pt()->set_halo(proc_in_charge);
            }

            // Send shared node ID and processor in charge info to
            // "intermediate" processor (i.e. the processor that has
            // the lookup schemes to talk to both
            send_data.push_back(j);
            send_data.push_back(proc_in_charge);
          }
        }
      }

      // End of data
      send_data.push_back(-1);

      // Find the number of data added to the vector
      send_n[domain] = send_data.size() - send_displacement[domain];
    }


    // Check if any processor has stumbled across overlooked halos
    // (if not we can omit the comms below)
    unsigned global_max_n_overlooked_halo = 0;
    MPI_Allreduce(&n_overlooked_halo,
                  &global_max_n_overlooked_halo,
                  1,
                  MPI_UNSIGNED,
                  MPI_MAX,
                  Comm_pt->mpi_comm());


    oomph_info << "Global max number of overlooked haloes: "
               << global_max_n_overlooked_halo << std::endl;

    if (Global_timings::Doc_comprehensive_timings)
    {
      tt_end = TimingHelpers::timer();
      oomph_info << "Time for setup 1st alltoalls in "
                    "Mesh::classify_halo_and_haloed_nodes: "
                 << tt_end - tt_start << std::endl;
      oomph_info.stream_pt()->flush();
      tt_start = TimingHelpers::timer();
    }

    // MemoryUsage::doc_memory_usage("after setup 1st alltoalls");

    // Any comms needed?
    if (global_max_n_overlooked_halo > 0)
    {
      // Storage for the number of data to be received from each processor
      Vector<int> receive_n(n_proc, 0);

      // Now send numbers of data to be sent between all processors
      MPI_Alltoall(
        &send_n[0], 1, MPI_INT, &receive_n[0], 1, MPI_INT, Comm_pt->mpi_comm());


      if (Global_timings::Doc_comprehensive_timings)
      {
        tt_end = TimingHelpers::timer();
        oomph_info
          << "Time for 1st alltoall in Mesh::classify_halo_and_haloed_nodes: "
          << tt_end - tt_start << std::endl;
        oomph_info.stream_pt()->flush();
        tt_start = TimingHelpers::timer();
      }

      // MemoryUsage::doc_memory_usage("after 1st alltoall");


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
                    Comm_pt->mpi_comm());


      if (Global_timings::Doc_comprehensive_timings)
      {
        tt_end = TimingHelpers::timer();
        oomph_info
          << "Time for 2nd alltoall in Mesh::classify_halo_and_haloed_nodes: "
          << tt_end - tt_start << std::endl;
        oomph_info.stream_pt()->flush();
        tt_start = TimingHelpers::timer();
      }

      // MemoryUsage::doc_memory_usage("after 2nd alltoall");

      // Provide storage for data to be sent to processor that used to be
      // in charge
      Vector<Vector<int>> send_data_for_proc_in_charge(n_proc);

      // Now use the received data
      for (int send_rank = 0; send_rank < n_proc; send_rank++)
      {
        // Don't bother to do anything for the processor corresponding to the
        // current processor or if no data were received from this processor
        if ((send_rank != my_rank) && (receive_n[send_rank] != 0))
        {
          // Counter for the data within the large array
          unsigned count = receive_displacement[send_rank];

          // Unpack until we reach "end of data" indicator (-1)
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
              // Shared halo node number in lookup scheme between intermediate
              // (i.e. this) processor and the one that has the overlooked halo
              unsigned j = unsigned(next_one);

              // Processor in charge:
              unsigned proc_in_charge = unsigned(receive_data[count++]);

              // Find actual node from shared node lookup scheme
              Node* nod_pt = shared_node_pt(send_rank, j);


              // Note: This search seems relatively cheap
              //       and in the tests done, did not benefit
              //       from conversion to map-based search
              //       as in
              //       TreeBasedRefineableMeshBase::synchronise_hanging_nodes()


              // Locate its index in lookup scheme with proc in charge
              bool found = false;
              unsigned nnod = nshared_node(proc_in_charge);
              for (unsigned jj = 0; jj < nnod; jj++)
              {
                if (nod_pt == shared_node_pt(proc_in_charge, jj))
                {
                  found = true;

                  // Shared node ID in lookup scheme with intermediate (i.e.
                  // this) processor
                  send_data_for_proc_in_charge[proc_in_charge].push_back(jj);

                  // Processor that holds the overlooked halo node
                  send_data_for_proc_in_charge[proc_in_charge].push_back(
                    send_rank);

                  break;
                }
              }
              if (!found)
              {
                std::ostringstream error_stream;
                error_stream
                  << "Failed to find node that is shared node " << j
                  << " (with processor " << send_rank
                  << ") \n in shared node lookup scheme with processor "
                  << proc_in_charge << " which is in charge.\n";
                throw OomphLibError(error_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
            }
          }
        }
      } // End of data is received


      if (Global_timings::Doc_comprehensive_timings)
      {
        tt_end = TimingHelpers::timer();
        oomph_info << "Time for 1st setup 3rd alltoall in "
                      "Mesh::classify_halo_and_haloed_nodes: "
                   << tt_end - tt_start << std::endl;
        oomph_info.stream_pt()->flush();
        tt_start = TimingHelpers::timer();
      }

      // MemoryUsage::doc_memory_usage("after 1st setup for 3rd alltoall");

      // Data to be sent to each processor
      Vector<int> all_send_n(n_proc, 0);

      // Storage for all values to be sent to all processors
      Vector<int> all_send_data;

      // Start location within send_data for data to be sent to each processor
      Vector<int> all_send_displacement(n_proc, 0);

      // Collate data
      for (int domain = 0; domain < n_proc; domain++)
      {
        // Set the offset for the current processor
        all_send_displacement[domain] = all_send_data.size();

        // Don't bother to do anything if the processor in the loop is the
        // current processor
        if (domain != my_rank)
        {
          unsigned n = send_data_for_proc_in_charge[domain].size();
          for (unsigned j = 0; j < n; j++)
          {
            all_send_data.push_back(send_data_for_proc_in_charge[domain][j]);
          }
        }

        // End of data
        all_send_data.push_back(-1);

        // Find the number of data added to the vector
        all_send_n[domain] =
          all_send_data.size() - all_send_displacement[domain];
      }


      if (Global_timings::Doc_comprehensive_timings)
      {
        tt_end = TimingHelpers::timer();
        oomph_info << "Time for 2nd setup 3rd alltoall in "
                      "Mesh::classify_halo_and_haloed_nodes: "
                   << tt_end - tt_start << std::endl;
        oomph_info.stream_pt()->flush();
        tt_start = TimingHelpers::timer();
      }

      // MemoryUsage::doc_memory_usage("after 2nd setup 3rd alltoall");

      // Storage for the number of data to be received from each processor
      Vector<int> all_receive_n(n_proc, 0);

      // Now send numbers of data to be sent between all processors
      MPI_Alltoall(&all_send_n[0],
                   1,
                   MPI_INT,
                   &all_receive_n[0],
                   1,
                   MPI_INT,
                   Comm_pt->mpi_comm());


      if (Global_timings::Doc_comprehensive_timings)
      {
        tt_end = TimingHelpers::timer();
        oomph_info
          << "Time for 3rd alltoall in Mesh::classify_halo_and_haloed_nodes: "
          << tt_end - tt_start << std::endl;
        oomph_info.stream_pt()->flush();
        tt_start = TimingHelpers::timer();
      }

      // MemoryUsage::doc_memory_usage("after 3rd alltoall");

      // We now prepare the data to be received
      // by working out the displacements from the received data
      Vector<int> all_receive_displacement(n_proc, 0);
      int all_receive_data_count = 0;

      for (int rank = 0; rank < n_proc; ++rank)
      {
        // Displacement is number of data received so far
        all_receive_displacement[rank] = all_receive_data_count;
        all_receive_data_count += all_receive_n[rank];
      }

      // Now resize the receive buffer for all data from all processors
      // Make sure that it has a size of at least one
      if (all_receive_data_count == 0)
      {
        ++all_receive_data_count;
      }
      Vector<int> all_receive_data(all_receive_data_count);

      // Make sure that the send buffer has size at least one
      // so that we don't get a segmentation fault
      if (all_send_data.size() == 0)
      {
        all_send_data.resize(1);
      }

      // Now send the data between all the processors
      MPI_Alltoallv(&all_send_data[0],
                    &all_send_n[0],
                    &all_send_displacement[0],
                    MPI_INT,
                    &all_receive_data[0],
                    &all_receive_n[0],
                    &all_receive_displacement[0],
                    MPI_INT,
                    Comm_pt->mpi_comm());


      if (Global_timings::Doc_comprehensive_timings)
      {
        tt_end = TimingHelpers::timer();
        oomph_info
          << "Time for 4th alltoall in Mesh::classify_halo_and_haloed_nodes: "
          << tt_end - tt_start << std::endl;
        oomph_info.stream_pt()->flush();
        tt_start = TimingHelpers::timer();
      }

      // MemoryUsage::doc_memory_usage("after 4th alltoall");


      // Now use the received data
      for (int send_rank = 0; send_rank < n_proc; send_rank++)
      {
        // Don't bother to do anything for the processor corresponding to the
        // current processor or if no data were received from this processor
        if ((send_rank != my_rank) && (all_receive_n[send_rank] != 0))
        {
          // Counter for the data within the large array
          unsigned count = all_receive_displacement[send_rank];

          // Unpack until we reach "end of data" indicator (-1)
          while (true)
          {
            // Read next entry
            int next_one = all_receive_data[count++];

            if (next_one == -1)
            {
              break;
            }
            else
            {
              // Shared node ID in lookup scheme with intermediate (sending)
              // processor
              unsigned j = unsigned(next_one);

              // Get pointer to previously overlooked halo
              Node* nod_pt = shared_node_pt(send_rank, j);

              // Proc where overlooked halo is
              unsigned proc_with_overlooked_halo = all_receive_data[count++];

              // Add it as being haloed from specified domain
              this->add_haloed_node_pt(proc_with_overlooked_halo, nod_pt);

              // Record new haloed node so it an be added to the shared
              // node lookup scheme (in a consistent order) below.
              over_looked_haloed_node_pt[proc_with_overlooked_halo].push_back(
                nod_pt);
            }
          }
        }
      } // End of data is received

      if (Global_timings::Doc_comprehensive_timings)
      {
        tt_end = TimingHelpers::timer();
        oomph_info << "Time for postproc 4th alltoall in "
                      "Mesh::classify_halo_and_haloed_nodes: "
                   << tt_end - tt_start << std::endl;
        oomph_info.stream_pt()->flush();
        tt_start = TimingHelpers::timer();
      }

      // MemoryUsage::doc_memory_usage("after postprocess 4th alltoall");

      // Now add previously overlooked halo/haloed nodes to shared node
      // lookup scheme in consistent order
      for (int d = 0; d < n_proc; d++)
      {
        // For all domains lower than the current domain: Do halos first
        // then haloed, to ensure correct order in lookup scheme from
        // the other side
        if (d < my_rank)
        {
          unsigned nnod = over_looked_halo_node_pt[d].size();
          for (unsigned j = 0; j < nnod; j++)
          {
            this->add_shared_node_pt(d, over_looked_halo_node_pt[d][j]);
          }
          nnod = over_looked_haloed_node_pt[d].size();
          for (unsigned j = 0; j < nnod; j++)
          {
            this->add_shared_node_pt(d, over_looked_haloed_node_pt[d][j]);
          }
        }
        else if (d > my_rank)
        {
          unsigned nnod = over_looked_haloed_node_pt[d].size();
          for (unsigned j = 0; j < nnod; j++)
          {
            this->add_shared_node_pt(d, over_looked_haloed_node_pt[d][j]);
          }
          nnod = over_looked_halo_node_pt[d].size();
          for (unsigned j = 0; j < nnod; j++)
          {
            this->add_shared_node_pt(d, over_looked_halo_node_pt[d][j]);
          }
        }
      }


      // Doc stats
      if (report_stats)
      {
        // Report total number of halo(ed) and shared nodes for this process
        oomph_info << "BEFORE SYNCHRONISE SHARED NODES Processor " << my_rank
                   << " holds " << this->nnode() << " nodes of which "
                   << this->nhalo_node() << " are halo nodes \n while "
                   << this->nhaloed_node() << " are haloed nodes, and "
                   << this->nshared_node() << " are shared nodes." << std::endl;

        // Report number of halo(ed) and shared nodes with each domain
        // from the current process
        for (int iproc = 0; iproc < n_proc; iproc++)
        {
          // Get vector of halo elements by copy operation
          Vector<GeneralisedElement*> halo_elem_pt(
            this->halo_element_pt(iproc));
          Vector<GeneralisedElement*> haloed_elem_pt(
            this->haloed_element_pt(iproc));
          oomph_info << "With process " << iproc << ", there are "
                     << this->nhalo_node(iproc) << " halo nodes, and "
                     << std::endl
                     << this->nhaloed_node(iproc) << " haloed nodes, and "
                     << this->nshared_node(iproc) << " shared nodes"
                     << std::endl
                     << halo_elem_pt.size() << " halo elements and "
                     << haloed_elem_pt.size() << " haloed elements\n";
        }
      }

      //  // Doc stats
      //  if (report_stats)
      //   {
      //    // Report total number of halo(ed) and shared nodes for this process
      //    oomph_info << "BEFORE SYNCHRONISE SHARED NODES Processor " <<
      //    my_rank
      //               << " holds " << this->nnode()
      //               << " nodes of which " << this->nhalo_node()
      //               << " are halo nodes \n while " << this->nhaloed_node()
      //               << " are haloed nodes, and " << this->nshared_node()
      //               << " are shared nodes." << std::endl;

      //    // Report number of halo(ed) and shared nodes with each domain
      //    // from the current process
      //    for (int iproc=0;iproc<n_proc;iproc++)
      //     {
      //      // Get vector of halo elements by copy operation
      //      Vector<GeneralisedElement*>
      //      halo_elem_pt(this->halo_element_pt(iproc));
      //      Vector<GeneralisedElement*> haloed_elem_pt(
      //       this->haloed_element_pt(iproc));
      //      oomph_info << "With process " << iproc << ", there are "
      //                 << this->nhalo_node(iproc) << " halo nodes, and " <<
      //                 std::endl
      //                 << this->nhaloed_node(iproc) << " haloed nodes, and "
      //                 << this->nshared_node(iproc) << " shared nodes" <<
      //                 std::endl
      //                 <<  halo_elem_pt.size() << " halo elements and "
      //                 <<  haloed_elem_pt.size() << " haloed elements\n";
      //     }
      //   }

    } // end if comms reqd because we encountered overlooked halo elements


    // MemoryUsage::doc_memory_usage("before sync halo nodes");

    // Synchronise shared nodes
    synchronise_shared_nodes(report_stats);

    // MemoryUsage::doc_memory_usage("after sync halo nodes");

#ifdef PARANOID
    // Has some twit filled in haloed nodes with own process?!
    if (Haloed_node_pt[my_rank].size() != 0)
    {
      throw OomphLibError(
        "Processor has haloed nodes with itself! Something's gone wrong!",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
    // Has some twit filled in halo nodes with own process?!
    if (Halo_node_pt[my_rank].size() != 0)
    {
      throw OomphLibError(
        "Processor has halo nodes with itself! Something's gone wrong!",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
    // Has some twit filled in root haloed elements with own process?!
    if (Root_haloed_element_pt[my_rank].size() != 0)
    {
      throw OomphLibError("Processor has root haloed elements with itself! "
                          "Something's gone wrong!",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // Has some twit filled in root halo elements with own process?!
    if (Root_halo_element_pt[my_rank].size() != 0)
    {
      throw OomphLibError(
        "Processor has root halo elements with itself! Something's gone wrong!",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Doc stats
    if (report_stats)
    {
      // Report total number of halo(ed) and shared nodes for this process
      oomph_info << "Processor " << my_rank << " holds " << this->nnode()
                 << " nodes of which " << this->nhalo_node()
                 << " are halo nodes \n while " << this->nhaloed_node()
                 << " are haloed nodes, and " << this->nshared_node()
                 << " are shared nodes." << std::endl;

      // Report number of halo(ed) and shared nodes with each domain
      // from the current process
      for (int iproc = 0; iproc < n_proc; iproc++)
      {
        // Get vector of halo elements by copy operation
        Vector<GeneralisedElement*> halo_elem_pt(this->halo_element_pt(iproc));
        Vector<GeneralisedElement*> haloed_elem_pt(
          this->haloed_element_pt(iproc));
        oomph_info << "With process " << iproc << ", there are "
                   << this->nhalo_node(iproc) << " halo nodes, and "
                   << std::endl
                   << this->nhaloed_node(iproc) << " haloed nodes, and "
                   << this->nshared_node(iproc) << " shared nodes" << std::endl
                   << halo_elem_pt.size() << " halo elements and "
                   << haloed_elem_pt.size() << " haloed elements\n";
      }
    }


    // MemoryUsage::doc_memory_usage("before resize halo nodes");

    // Now resize halo nodes if required (can be over-ruled from the outside
    // either by user (for (risky!) efficienty saving) or from overloaded
    // version of classify_... in refineable version of that function
    // where resize_halo_nodes() is called after synchronising hanging nodes.
    if (!Resize_halo_nodes_not_required)
    {
      resize_halo_nodes();
    }

    // MemoryUsage::doc_memory_usage("after resize halo nodes");

    if (Global_timings::Doc_comprehensive_timings)
    {
      double t_end = TimingHelpers::timer();
      oomph_info << "Total time for Mesh::classify_halo_and_halo_nodes(): "
                 << t_end - t_start << std::endl;
      oomph_info.stream_pt()->flush();
    }

    //  MemoryUsage::doc_memory_usage(
    //   "at end of Mesh::classify_halo_and_halo_nodes()");
  }


  //========================================================================
  /// Helper function that resizes halo nodes to the same
  /// size as their non-halo counterparts if required. (A discrepancy
  /// can arise if a FaceElement that introduces additional unknowns
  /// are attached to a bulk element that shares a node with a haloed element.
  /// In that case the joint node between haloed and non-haloed element
  /// is resized on that processor but not on the one that holds the
  /// halo counterpart (because no FaceElement is attached to the halo
  /// element)
  //=========================================================================
  void Mesh::resize_halo_nodes()
  {
    double t_start = 0.0;
    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start = TimingHelpers::timer();
    }

    MPI_Status status;

    // Nuffink needs to be done if mesh isn't distributed
    if (is_mesh_distributed())
    {
      // Storage for current processor and number of processors
      int n_proc = Comm_pt->nproc();
      int my_rank = Comm_pt->my_rank();

      // Loop over domains on which non-halo counter parts of my halo nodes live
      for (int d = 0; d < n_proc; d++)
      {
        // On current processor: Receive data for my halo nodes with proc d.
        // Elsewhere: Send haloed data with proc d.
        if (d == my_rank)
        {
          // Loop over domains that hold non-halo counterparts of my halo nodes
          for (int dd = 0; dd < n_proc; dd++)
          {
            // Don't talk to yourself
            if (dd != d)
            {
              // How many of my nodes are halo nodes whose non-halo
              // counterpart is located on processor dd?
              int nnod_halo = this->nhalo_node(dd);
              int nnod_ext_halo = this->nexternal_halo_node(dd);
              if ((nnod_halo + nnod_ext_halo) != 0)
              {
                // Receive from processor dd number of haloed nodes (entry 0),
                // external haloed nodes (entry 1) and total number of
                // unsigneds to be sent below (entry 2)
                Vector<int> tmp(3);
                MPI_Recv(
                  &tmp[0], 3, MPI_INT, dd, 0, Comm_pt->mpi_comm(), &status);

#ifdef PARANOID
                // Check that number of halo/haloed nodes match
                int nnod_haloed = tmp[0];
                if (nnod_haloed != nnod_halo)
                {
                  std::ostringstream error_message;
                  error_message << "Clash in numbers of halo and haloed nodes "
                                << std::endl;
                  error_message << " between procs " << dd << " and " << d
                                << ": " << nnod_haloed << " " << nnod_halo
                                << std::endl;
                  throw OomphLibError(error_message.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }

                // Check that number of external halo/haloed nodes match
                int nnod_ext_haloed = tmp[1];
                if (nnod_ext_haloed != nnod_ext_halo)
                {
                  std::ostringstream error_message;
                  error_message
                    << "Clash in numbers of external halo and haloed nodes "
                    << std::endl;
                  error_message << " between procs " << dd << " and " << d
                                << ": " << nnod_ext_haloed << " "
                                << nnod_ext_halo << std::endl;
                  throw OomphLibError(error_message.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
#endif

                // How many unsigneds are we about to receive
                unsigned n_rec = tmp[2];

                // Get strung-together data from other proc
                Vector<unsigned> unsigned_rec_data(n_rec);
                MPI_Recv(&unsigned_rec_data[0],
                         n_rec,
                         MPI_UNSIGNED,
                         dd,
                         0,
                         Comm_pt->mpi_comm(),
                         &status);

                // Step through the flat-packed unsigneds
                unsigned count = 0;

                // Normal and external halo nodes
                for (unsigned loop = 0; loop < 2; loop++)
                {
                  unsigned hi_nod = nnod_halo;
                  if (loop == 1)
                  {
                    hi_nod = nnod_ext_halo;
                  }
                  for (int j = 0; j < int(hi_nod); j++)
                  {
                    // Which node are we dealing with
                    Node* nod_pt = 0;
                    if (loop == 0)
                    {
                      nod_pt = this->halo_node_pt(dd, j);
                    }
                    else
                    {
                      nod_pt = this->external_halo_node_pt(dd, j);
                    }

                    // How many values do we have locally?
                    unsigned nval_local = nod_pt->nvalue();

                    // Read number of values on other side
                    unsigned nval_other = unsigned_rec_data[count++];

                    if (nval_local != nval_other)
                    {
                      nod_pt->resize(nval_other);
                    }

                    // Read number of entries in resize map
                    unsigned nentry = unsigned_rec_data[count++];
                    if (nentry != 0)
                    {
                      // Is current node a boundary node?
                      BoundaryNodeBase* bnod_pt =
                        dynamic_cast<BoundaryNodeBase*>(nod_pt);
#ifdef PARANOID
                      if (bnod_pt == 0)
                      {
                        throw OomphLibError(
                          "Failed to cast node to boundary node even though "
                          "we've received data for boundary node",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
                      }
#endif

                      // Create storage for map if it doesn't already exist
                      bool already_existed = true;
                      if (
                        bnod_pt
                          ->index_of_first_value_assigned_by_face_element_pt() ==
                        0)
                      {
                        bnod_pt
                          ->index_of_first_value_assigned_by_face_element_pt() =
                          new std::map<unsigned, unsigned>;
                        already_existed = false;
                      }

                      // Get pointer to the map of indices associated with
                      // additional values created by face elements
                      std::map<unsigned, unsigned>* map_pt =
                        bnod_pt
                          ->index_of_first_value_assigned_by_face_element_pt();

                      // Loop over number of entries in map (as received)
                      for (unsigned i = 0; i < nentry; i++)
                      {
                        // Read out pairs...
                        unsigned first_received = unsigned_rec_data[count++];
                        unsigned second_received = unsigned_rec_data[count++];

                        // If it exists check that values are consistent:
                        if (already_existed)
                        {
#ifdef PARANOID
                          if ((*map_pt)[first_received] != second_received)
                          {
                            std::ostringstream error_message;
                            error_message << "Existing map entry for map entry "
                                          << i << " for node located at ";
                            unsigned n = nod_pt->ndim();
                            for (unsigned ii = 0; ii < n; ii++)
                            {
                              error_message << nod_pt->position(ii) << " ";
                            }
                            error_message
                              << "Key: " << first_received << " "
                              << "Local value: " << (*map_pt)[first_received]
                              << " "
                              << "Received value: " << second_received
                              << std::endl;
                            throw OomphLibError(error_message.str(),
                                                OOMPH_CURRENT_FUNCTION,
                                                OOMPH_EXCEPTION_LOCATION);
                          }
#endif
                        }
                        // Else assign
                        else
                        {
                          (*map_pt)[first_received] = second_received;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        // Send my haloed nodes whose halo counterparts are located on processor
        // d
        else
        {
          // Storage for number of haloed nodes (entry 0), number of external
          // haloed nodes (entry 1) and total number of unsigneds to be sent
          // below (entry 2)
          Vector<int> tmp(3);
          int nnod_haloed = this->nhaloed_node(d);
          tmp[0] = nnod_haloed;
          int nnod_ext_haloed = this->nexternal_haloed_node(d);
          tmp[1] = nnod_ext_haloed;
          if ((nnod_haloed + nnod_ext_haloed) != 0)
          {
            // Now string together the data
            Vector<unsigned> unsigned_send_data;

            // Normal and external haloed nodes
            for (unsigned loop = 0; loop < 2; loop++)
            {
              unsigned hi_nod = nnod_haloed;
              if (loop == 1)
              {
                hi_nod = nnod_ext_haloed;
              }
              for (int j = 0; j < int(hi_nod); j++)
              {
                // Which node are we dealing with?
                Node* nod_pt = 0;
                if (loop == 0)
                {
                  nod_pt = this->haloed_node_pt(d, j);
                }
                else
                {
                  nod_pt = this->external_haloed_node_pt(d, j);
                }

                // Add number of values of node
                unsigned_send_data.push_back(nod_pt->nvalue());

                // Get pointer to the map of indices associated with
                // additional values created by face elements

                // Is it a boundary node?
                BoundaryNodeBase* bnod_pt =
                  dynamic_cast<BoundaryNodeBase*>(nod_pt);
                if (bnod_pt == 0)
                {
                  // Not a boundary node -- there are zero map entries to follow
                  unsigned_send_data.push_back(0);
                }
                else
                {
                  // It's a boundary node: Check if it's been resized
                  std::map<unsigned, unsigned>* map_pt =
                    bnod_pt->index_of_first_value_assigned_by_face_element_pt();

                  // No additional values created -- there are zero map entries
                  // to follow
                  if (map_pt == 0)
                  {
                    unsigned_send_data.push_back(0);
                  }
                  // Created additional values
                  else
                  {
                    // How many map entries were there
                    unsigned_send_data.push_back(map_pt->size());

                    // Loop over entries in map and add to send data
                    for (std::map<unsigned, unsigned>::iterator p =
                           map_pt->begin();
                         p != map_pt->end();
                         p++)
                    {
                      unsigned_send_data.push_back((*p).first);
                      unsigned_send_data.push_back((*p).second);
                    }
                  }
                }
              }
            }

            // How many values are there in total?
            int n_send = unsigned_send_data.size();
            tmp[2] = n_send;

            // Send the counts across to the other processor
            MPI_Send(&tmp[0], 3, MPI_INT, d, 0, Comm_pt->mpi_comm());

            // Send it across to the processor
            MPI_Send(&unsigned_send_data[0],
                     n_send,
                     MPI_UNSIGNED,
                     d,
                     0,
                     Comm_pt->mpi_comm());
          }
        }
      }
    }

    if (Global_timings::Doc_comprehensive_timings)
    {
      double t_end = TimingHelpers::timer();
      oomph_info << "Total time for Mesh::resize_halo_nodes(): "
                 << t_end - t_start << std::endl;
    }
  }


  //========================================================================
  /// Get all the halo data stored in the mesh and add pointers to
  /// the data to the map, indexed by global equation number
  //=========================================================================
  void Mesh::get_all_halo_data(std::map<unsigned, double*>& map_of_halo_data)
  {
    // Loop over the map of Halo_nodes
    for (std::map<unsigned, Vector<Node*>>::iterator it = Halo_node_pt.begin();
         it != Halo_node_pt.end();
         ++it)
    {
      // Find the number of nodes
      unsigned n_node = (it->second).size();
      // Loop over them all
      for (unsigned n = 0; n < n_node; n++)
      {
        // Add the Node's values (including any solid values) to
        // the map
        (it->second)[n]->add_value_pt_to_map(map_of_halo_data);
      }
    } // End of loop over halo nodes

    // Now loop over all the halo elements and add their internal data
    // Loop over the map of Halo_nodes
    for (std::map<unsigned, Vector<GeneralisedElement*>>::iterator it =
           Root_halo_element_pt.begin();
         it != Root_halo_element_pt.end();
         ++it)
    {
      // Find the number of root elements
      unsigned n_element = (it->second).size();
      for (unsigned e = 0; e < n_element; e++)
      {
        GeneralisedElement* el_pt = (it->second)[e];

        // Is it a refineable element?
        RefineableElement* ref_el_pt = dynamic_cast<RefineableElement*>(el_pt);
        if (ref_el_pt != 0)
        {
          // Vector of pointers to leaves in tree emanating from
          // current root halo element
          Vector<Tree*> leaf_pt;
          ref_el_pt->tree_pt()->stick_leaves_into_vector(leaf_pt);

          // Loop over leaves and add their objects (the finite elements)
          // to vector
          unsigned nleaf = leaf_pt.size();
          for (unsigned l = 0; l < nleaf; l++)
          {
            leaf_pt[l]->object_pt()->add_internal_value_pt_to_map(
              map_of_halo_data);
          }
        }
        else
        {
          el_pt->add_internal_value_pt_to_map(map_of_halo_data);
        }
      }
    }

    /// Repeat for the external data
    for (std::map<unsigned, Vector<Node*>>::iterator it =
           External_halo_node_pt.begin();
         it != External_halo_node_pt.end();
         ++it)
    {
      // Find the number of nodes
      unsigned n_node = (it->second).size();
      // Loop over them all
      for (unsigned n = 0; n < n_node; n++)
      {
        // Add the Node's values (including any solid values) to
        // the map
        (it->second)[n]->add_value_pt_to_map(map_of_halo_data);
      }
    } // End of loop over halo nodes

    // Now loop over all the halo elements and add their internal data
    // Loop over the map of Halo_nodes
    for (std::map<unsigned, Vector<GeneralisedElement*>>::iterator it =
           External_halo_element_pt.begin();
         it != External_halo_element_pt.end();
         ++it)
    {
      // Find the number of root elements
      unsigned n_element = (it->second).size();
      for (unsigned e = 0; e < n_element; e++)
      {
        (it->second)[e]->add_internal_value_pt_to_map(map_of_halo_data);
      }
    }
  }


  //========================================================================
  /// Get halo node stats for this distributed mesh:
  /// Average/max/min number of halo nodes over all processors.
  /// \b Careful: Involves MPI Broadcasts and must therefore
  /// be called on all processors!
  //========================================================================
  void Mesh::get_halo_node_stats(double& av_number,
                                 unsigned& max_number,
                                 unsigned& min_number)
  {
    // Storage for number of processors and current processor
    int n_proc = Comm_pt->nproc();
    int my_rank = Comm_pt->my_rank();

    // Create vector to hold number of halo nodes
    Vector<int> nhalo_nodes(n_proc);

    // Stick own number of halo nodes into appropriate entry
    int nhalo_node_local = nhalo_node();

    // Gather information on root processor: First argument group
    // specifies what is to be sent (one int from each procssor, indicating
    // the number of dofs on it), the second group indicates where
    // the results are to be gathered (in rank order) on root processor.
    MPI_Gather(&nhalo_node_local,
               1,
               MPI_INT,
               &nhalo_nodes[0],
               1,
               MPI_INT,
               0,
               Comm_pt->mpi_comm());

    // Initialise stats
    av_number = 0.0;
    int max = -1;
    int min = 1000000000;

    if (my_rank == 0)
    {
      for (int i = 0; i < n_proc; i++)
      {
        av_number += double(nhalo_nodes[i]);
        if (int(nhalo_nodes[i]) > max) max = nhalo_nodes[i];
        if (int(nhalo_nodes[i]) < min) min = nhalo_nodes[i];
      }
      av_number /= double(n_proc);
    }

    // Now broadcast the result back out
    MPI_Bcast(&max, 1, MPI_INT, 0, Comm_pt->mpi_comm());
    MPI_Bcast(&min, 1, MPI_INT, 0, Comm_pt->mpi_comm());
    MPI_Bcast(&av_number, 1, MPI_DOUBLE, 0, Comm_pt->mpi_comm());

    max_number = max;
    min_number = min;
  }


  //========================================================================
  /// Get haloed node stats for this distributed mesh:
  /// Average/max/min number of haloed nodes over all processors.
  /// \b Careful: Involves MPI Broadcasts and must therefore
  /// be called on all processors!
  //========================================================================
  void Mesh::get_haloed_node_stats(double& av_number,
                                   unsigned& max_number,
                                   unsigned& min_number)
  {
    // Storage for number of processors and current processor
    int n_proc = Comm_pt->nproc();
    int my_rank = Comm_pt->my_rank();

    // Create vector to hold number of haloed nodes
    Vector<int> nhaloed_nodes(n_proc);

    // Stick own number of haloed nodes into appropriate entry
    int nhaloed_node_local = nhaloed_node();

    // Gather information on root processor: First argument group
    // specifies what is to be sent (one int from each procssor, indicating
    // the number of dofs on it), the second group indicates where
    // the results are to be gathered (in rank order) on root processor.
    MPI_Gather(&nhaloed_node_local,
               1,
               MPI_INT,
               &nhaloed_nodes[0],
               1,
               MPI_INT,
               0,
               Comm_pt->mpi_comm());

    // Initialise stats
    av_number = 0.0;
    int max = -1;
    int min = 1000000000;

    if (my_rank == 0)
    {
      for (int i = 0; i < n_proc; i++)
      {
        av_number += double(nhaloed_nodes[i]);
        if (int(nhaloed_nodes[i]) > max) max = nhaloed_nodes[i];
        if (int(nhaloed_nodes[i]) < min) min = nhaloed_nodes[i];
      }
      av_number /= double(n_proc);
    }

    // Now broadcast the result back out
    MPI_Bcast(&max, 1, MPI_INT, 0, Comm_pt->mpi_comm());
    MPI_Bcast(&min, 1, MPI_INT, 0, Comm_pt->mpi_comm());
    MPI_Bcast(&av_number, 1, MPI_DOUBLE, 0, Comm_pt->mpi_comm());

    max_number = max;
    min_number = min;
  }

  //========================================================================
  /// Distribute the mesh. Add to vector of deleted elements.
  //========================================================================
  void Mesh::distribute(OomphCommunicator* comm_pt,
                        const Vector<unsigned>& element_domain,
                        Vector<GeneralisedElement*>& deleted_element_pt,
                        DocInfo& doc_info,
                        const bool& report_stats,
                        const bool& overrule_keep_as_halo_element_status)
  {
    // Store communicator
    Comm_pt = comm_pt;

    // Storage for number of processors and current processor
    int n_proc = comm_pt->nproc();
    int my_rank = comm_pt->my_rank();

    // Storage for number of elements and number of nodes on this mesh
    unsigned nelem = this->nelement();
    unsigned nnod = this->nnode();

    std::ostringstream filename;

    // Doc the partitioning (only on processor 0)
    //-------------------------------------------
    if (doc_info.is_doc_enabled())
    {
      if (my_rank == 0)
      {
        // Open files for doc of element partitioning
        Vector<std::ofstream*> domain_file(n_proc);
        for (int d = 0; d < n_proc; d++)
        {
          // Note: doc_info.number() was set in Problem::distribute(...) to
          // reflect the submesh number
          // Clear the filename
          filename.str("");
          filename << doc_info.directory() << "/domain" << d << "-"
                   << doc_info.number() << ".dat";
          domain_file[d] = new std::ofstream(filename.str().c_str());
        }

        // Doc
        for (unsigned e = 0; e < nelem; e++)
        {
          // If we can't cast to a finite element, we can't output because
          // there is no output function
          FiniteElement* f_el_pt =
            dynamic_cast<FiniteElement*>(this->element_pt(e));
          if (f_el_pt != 0)
          {
            f_el_pt->output(*domain_file[element_domain[e]], 5);
          }
        }

        for (int d = 0; d < n_proc; d++)
        {
          domain_file[d]->close();
          delete domain_file[d];
          domain_file[d] = 0;
        }
      }
    }

    // Loop over all elements, associate all
    //--------------------------------------
    // nodes with the highest-numbered processor and record all
    //---------------------------------------------------------
    // processors the node is associated with
    //---------------------------------------

    // Storage for processors in charge and processors associated with data
    std::map<Data*, std::set<unsigned>> processors_associated_with_data;
    std::map<Data*, unsigned> processor_in_charge;

    // For all nodes set the processor in charge to zero
    for (unsigned j = 0; j < nnod; j++)
    {
      Node* nod_pt = this->node_pt(j);
      processor_in_charge[nod_pt] = 0;
    }

    // Loop over elements
    for (unsigned e = 0; e < nelem; e++)
    {
      // Get an element and its domain
      FiniteElement* el_pt = dynamic_cast<FiniteElement*>(this->element_pt(e));
      // Element will only have nodes if it is a finite element
      if (el_pt != 0)
      {
        unsigned el_domain = element_domain[e];

        // Associate nodes with highest numbered processor
        unsigned nnod = el_pt->nnode();
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = el_pt->node_pt(j);

          // processor in charge was initialised to 0 above
          if (el_domain > processor_in_charge[nod_pt])
          {
            processor_in_charge[nod_pt] = el_domain;
          }
          processors_associated_with_data[nod_pt].insert(el_domain);
        }
      }
    }

    // Doc the partitioning (only on processor 0)
    //-------------------------------------------
    if (doc_info.is_doc_enabled())
    {
      if (my_rank == 0)
      {
        // Open files for doc of node partitioning
        Vector<std::ofstream*> node_file(n_proc);
        for (int d = 0; d < n_proc; d++)
        {
          // Note: doc_info.number() was set in Problem::distribute(...) to
          // reflect the submesh number...
          // Clear the filename
          filename.str("");
          filename << doc_info.directory() << "/node" << d << "-"
                   << doc_info.number() << ".dat";
          node_file[d] = new std::ofstream(filename.str().c_str());
        }

        // Doc
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = this->node_pt(j);
          const unsigned n_dim = nod_pt->ndim();
          for (unsigned i = 0; i < n_dim; i++)
          {
            *node_file[processor_in_charge[nod_pt]] << nod_pt->x(i) << " ";
          }
          *node_file[processor_in_charge[nod_pt]] << "\n";
        }
        for (int d = 0; d < n_proc; d++)
        {
          node_file[d]->close();
          delete node_file[d];
          node_file[d] = 0;
        }
      }
    }

    // Declare all nodes as obsolete. We'll
    // change this setting for all nodes that must be retained
    // further down
    for (unsigned j = 0; j < nnod; j++)
    {
      this->node_pt(j)->set_obsolete();
    }


    // Backup old mesh data and flush mesh
    //-------------------------------------

    // Backup pointers to elements in this mesh
    Vector<GeneralisedElement*> backed_up_el_pt(nelem);
    for (unsigned e = 0; e < nelem; e++)
    {
      backed_up_el_pt[e] = this->element_pt(e);
    }

    // Info. used to re-generate the boundary element scheme after the
    // deletion of elements not belonging to the current processor)

    // Get the number of boundary elements before the deletion of non
    // retained elements
    const unsigned tmp_nboundary = this->nboundary();
    Vector<unsigned> ntmp_boundary_elements(tmp_nboundary);
    // In case that there are regions, also get the number of boundary
    // elements in each region
    Vector<Vector<unsigned>> ntmp_boundary_elements_in_region(tmp_nboundary);
    // Also get the finite element version of the elements and back them
    // up
    Vector<FiniteElement*> backed_up_f_el_pt(nelem);

    // Only do this if the mesh is a TriangleMeshBase
    TriangleMeshBase* triangle_mesh_pt = dynamic_cast<TriangleMeshBase*>(this);
    bool is_a_triangle_mesh_base_mesh = false;
    if (triangle_mesh_pt != 0)
    {
      // Set the flag to indicate we are working with a TriangleMeshBase
      // mesh
      is_a_triangle_mesh_base_mesh = true;

      // Are there regions?
      const unsigned n_regions = triangle_mesh_pt->nregion();

      // Loop over the boundaries
      for (unsigned ib = 0; ib < tmp_nboundary; ib++)
      {
        // Get the number of boundary elements
        ntmp_boundary_elements[ib] = this->nboundary_element(ib);

        // Resize the container
        ntmp_boundary_elements_in_region[ib].resize(n_regions);

        // Loop over the regions
        for (unsigned rr = 0; rr < n_regions; rr++)
        {
          // Get the region id
          const unsigned region_id =
            static_cast<unsigned>(triangle_mesh_pt->region_attribute(rr));

          // Store the number of element in the region (notice we are
          // using the region index not the region id to refer to the
          // region)
          ntmp_boundary_elements_in_region[ib][rr] =
            triangle_mesh_pt->nboundary_element_in_region(ib, region_id);

        } // for (rr < n_regions)

      } // for (ib < tmp_nboundary)

      for (unsigned e = 0; e < nelem; e++)
      {
        // Get the finite element version of the elements and back them
        // up
        backed_up_f_el_pt[e] = this->finite_element_pt(e);
      }

    } // if (triangle_mesh_pt!=0)

    // Flush the mesh storage
    this->flush_element_storage();

    // Delete any storage of external elements and nodes
    this->delete_all_external_storage();

    // Boolean to indicate which element is to be retained
    std::vector<bool> element_retained(nelem, false);

    // Storage for element numbers of root halo elements that will be
    // retained on current processor: root_halo_element[p][j]
    // stores the element number (in the order in which the elements are stored
    // in backed_up_el_pt) of the j-th root halo element with processor
    // p.
    Vector<Vector<int>> root_halo_element(n_proc);

    // Dummy entry to make sure we always have something to send
    for (int p = 0; p < n_proc; p++)
    {
      root_halo_element[p].push_back(-1);
    }

    // Determine which elements are going to end up on which processor
    //----------------------------------------------------------------
    unsigned number_of_retained_elements = 0;

    // Loop over all backed up elements
    nelem = backed_up_el_pt.size();
    for (unsigned e = 0; e < nelem; e++)
    {
      // Get element and its domain
      GeneralisedElement* el_pt = backed_up_el_pt[e];
      unsigned el_domain = element_domain[e];

      // If element is located on current processor add it back to the mesh
      if (el_domain == unsigned(my_rank))
      {
        // Add element to current processor
        element_retained[e] = true;
        number_of_retained_elements++;
      }
      // Otherwise we may still need it if it's a halo element:
      else
      {
        // If this current mesh has been told to keep all elements as halos,
        // OR the element itself knows that it must be kept then
        // keep it
        if ((this->Keep_all_elements_as_halos) ||
            (el_pt->must_be_kept_as_halo()))
        {
          if (!overrule_keep_as_halo_element_status)
          {
            // Add as root halo element whose non-halo counterpart is
            // located on processor el_domain
            if (!element_retained[e])
            {
              root_halo_element[el_domain].push_back(e);
              element_retained[e] = true;
              number_of_retained_elements++;
            }
          }
        }
        // Otherwise: Is one of the nodes associated with the current processor?
        else
        {
          // Can only have nodes if this is a finite element
          FiniteElement* finite_el_pt = dynamic_cast<FiniteElement*>(el_pt);
          if (finite_el_pt != 0)
          {
            unsigned n_node = finite_el_pt->nnode();
            for (unsigned n = 0; n < n_node; n++)
            {
              Node* nod_pt = finite_el_pt->node_pt(n);

              // Keep element? (use stl find?)
              unsigned keep_it = false;
              for (std::set<unsigned>::iterator it =
                     processors_associated_with_data[nod_pt].begin();
                   it != processors_associated_with_data[nod_pt].end();
                   it++)
              {
                if (*it == unsigned(my_rank))
                {
                  keep_it = true;
                  // Break out of the loop over processors
                  break;
                }
              }

              // Add a root halo element either if keep_it=true
              if (keep_it)
              {
                // Add as root halo element whose non-halo counterpart is
                // located on processor el_domain
                if (!element_retained[e])
                {
                  root_halo_element[el_domain].push_back(e);
                  element_retained[e] = true;
                  number_of_retained_elements++;
                }
                // Now break out of loop over nodes
                break;
              }
            }
          }
        } // End of testing for halo by virtue of shared nodes
      } // End of halo element conditions
    } // end of loop over elements

    // First check that the number of elements is greater than zero, when
    // working with submeshes it may be possible that some of them have
    // no elements (face element meshes) since those may have been
    // deleted in "Problem::actions_before_distribute()"
    if (nelem > 0)
    {
      // Check that we are working with a TriangleMeshBase mesh, if
      // that is the case then we need to create shared boundaries
      if (is_a_triangle_mesh_base_mesh)
      {
        // Creation of shared boundaries
        // ------------------------------
        // All processors share the same boundary id for the created
        // shared boundary. We need all the elements on all processors,
        // that is why this step is performed before the deletion of the
        // elements not associated to the current processor.
        // N.B.: This applies only to unstructured meshes
        this->create_shared_boundaries(comm_pt,
                                       element_domain,
                                       backed_up_el_pt,
                                       backed_up_f_el_pt,
                                       processors_associated_with_data,
                                       overrule_keep_as_halo_element_status);
      } // if (is_a_triangle_mesh_base_mesh)
    } // if (nelem > 0)

    // NOTE: No need to add additional layer of halo elements.
    //       Procedure for removing "overlooked" halo nodes in
    //       deals classify_halo_and_haloed_nodes() deals
    //       with the problem addressed here. [code that dealt with this
    //       problem at distribution stage has been removed]

    // Store the finite element pointer version of the elements that are
    // about to be deleted, used to reset the boundary elements info
    Vector<FiniteElement*> deleted_f_element_pt;

    // Copy the elements associated with the actual
    // current processor into its own permanent storage.
    // Do it in the order in which the elements appeared originally
    nelem = backed_up_el_pt.size();
    for (unsigned e = 0; e < nelem; e++)
    {
      GeneralisedElement* el_pt = backed_up_el_pt[e];
      if (element_retained[e])
      {
        this->add_element_pt(el_pt);
      }
      else
      {
        // Flush the object attached to the tree for this element?
        RefineableElement* ref_el_pt = dynamic_cast<RefineableElement*>(el_pt);
        if (ref_el_pt != 0)
        {
          ref_el_pt->tree_pt()->flush_object();
        }


        // Store pointer to the element that's about to be deleted.

        // Only for structured meshes since this "deleted_element_pt"
        // vector is used in the "problem" class to set null pointer to
        // the deleted elements in the Base_mesh_element_pt structure
        if (!is_a_triangle_mesh_base_mesh)
        {
          deleted_element_pt.push_back(el_pt);
        } // if (!is_a_triangle_mesh_base_mesh)

        if (is_a_triangle_mesh_base_mesh)
        {
          // Store pointer to the finite element that's about to be deleted
          deleted_f_element_pt.push_back(backed_up_f_el_pt[e]);
        }

        // Delete the element
        delete el_pt;
      }
    }

    // Copy the root halo elements associated with the actual
    // current processor into its own permanent storage; the order
    // here is somewhat random but we compensate for that by
    // ensuring that the corresponding haloed elements are
    // added in the same order below
#ifdef PARANOID
    std::map<unsigned, bool> done;
#endif
    for (int d = 0; d < n_proc; d++)
    {
      nelem = root_halo_element[d].size();
      for (unsigned e = 0; e < nelem; e++)
      {
        int number = root_halo_element[d][e];
        if (number >= 0)
        {
#ifdef PARANOID
          if (done[number])
          {
            std::ostringstream error_message;
            error_message << "Have already added element " << number
                          << " as root halo element\n"
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
          done[number] = true;
#endif
          this->add_root_halo_element_pt(d, backed_up_el_pt[number]);
        }
      }
    }


    // Now get root haloed elements: root_haloed_element[p][j] stores
    // the element number (in the order in which the elements are stored
    // in backedup_el_pt) of the j-th rooted halo element with processor
    // p. On proc my_proc this the same as root_haloed_element[my_proc][j]
    // on processor p, so get the information by gatherv operations.
    Vector<Vector<unsigned>> root_haloed_element(n_proc);

    // Find out number of root halo elements with each other proc
    Vector<int> nhalo(n_proc, 0);
    Vector<int> nhaloed(n_proc, 0);
    for (int p = 0; p < n_proc; p++)
    {
      nhalo[p] = root_halo_element[p].size();
    }

    // Each processor sends number of halo elements it has with processor
    // p to processor p where this information is stored in nhaloed[...]
    for (int p = 0; p < n_proc; p++)
    {
      // Gather the p-th entries in nhalo from every processor on
      // processor p and store them in nhaloed consecutively
      // starting at beginning
      MPI_Gather(
        &nhalo[p], 1, MPI_INT, &nhaloed[0], 1, MPI_INT, p, comm_pt->mpi_comm());
    }

    // In the big sequence of concatenated root halo elements (enumerated
    // individually on the various processors) where do the root halo
    // elements from a given processor start? Also figure out how many
    // root haloed elements there are in total by summing up their numbers
    Vector<int> start_index(n_proc, 0);
    unsigned total_number_of_root_haloed_elements = 0;
    for (int i_proc = 0; i_proc < n_proc; i_proc++)
    {
      total_number_of_root_haloed_elements += nhaloed[i_proc];
      if (i_proc != 0)
      {
        start_index[i_proc] =
          total_number_of_root_haloed_elements - nhaloed[i_proc];
      }
      else
      {
        start_index[0] = 0;
      }
    }

    // Storage for all root haloed elements from the various processors, one
    // after the other, with some padding from negative entries to avoid
    // zero length vectors
    Vector<int> all_root_haloed_element(total_number_of_root_haloed_elements,
                                        0);

    // Now send the ids of the relevant elements via gatherv
    for (int p = 0; p < n_proc; p++)
    {
      // Gather the p-th entries in nhalo from every processor on
      // processor p and store them in nhaloed consecutively
      // starting at beginning
      MPI_Gatherv(&root_halo_element[p][0], // pointer to first entry in vector
                                            // to be gathered on processor p
                  nhalo[p], // Number of entries to be sent
                  MPI_INT,
                  &all_root_haloed_element[0], // Target -- this will store
                                               // the element numbers of
                                               // all root haloed elements
                                               // received from other processors
                                               // in order
                  &nhaloed[0], // Pointer to the vector containing the lengths
                               // of the vectors received from elsewhere
                  &start_index[0], // "offset" for storage of vector received
                                   // from various processors in the global
                                   // concatenated vector
                  MPI_INT,
                  p, // processor that gathers the information
                  comm_pt->mpi_comm());
    }


    // Determine root haloed elements
    //-------------------------------

    // Loop over all other processors
    unsigned count = 0;
    for (int d = 0; d < n_proc; d++)
    {
#ifdef PARANOID
      std::map<unsigned, bool> done;
#endif

      // Loop over root haloed elements with specified processor
      unsigned n = nhaloed[d];
      for (unsigned e = 0; e < n; e++)
      {
        int number = all_root_haloed_element[count];
        count++;

        // Ignore padded -1s that were only inserted to avoid
        // zero sized vectors
        if (number >= 0)
        {
          // Get pointer to element
          GeneralisedElement* el_pt = backed_up_el_pt[number];

          // Halo elements can't be haloed themselves
          if (!el_pt->is_halo())
          {
#ifdef PARANOID
            if (done[number])
            {
              std::ostringstream error_message;
              error_message << "Have already added element " << number
                            << " as root haloed element\n"
                            << std::endl;
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
            done[number] = true;
#endif

            // Current element is haloed by other processor
            this->add_root_haloed_element_pt(d, el_pt);
          }
        }
      }
    }


    // Doc stats
    if (report_stats)
    {
      oomph_info << "Processor " << my_rank << " holds " << this->nelement()
                 << " elements of which " << this->nroot_halo_element()
                 << " are root halo elements \n while "
                 << this->nroot_haloed_element() << " are root haloed elements"
                 << std::endl;
    }


    // Loop over all retained elements and mark their nodes
    //-----------------------------------------------------
    // as to be retained too (some double counting going on here)
    //-----------------------------------------------------------
    nelem = this->nelement();
    for (unsigned e = 0; e < nelem; e++)
    {
      FiniteElement* f_el_pt =
        dynamic_cast<FiniteElement*>(this->element_pt(e));

      // If we have a finite element
      if (f_el_pt != 0)
      {
        // Loop over nodes
        unsigned nnod = f_el_pt->nnode();
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = f_el_pt->node_pt(j);
          nod_pt->set_non_obsolete();
        }
      }
    }


    // Now remove the pruned nodes
    this->prune_dead_nodes();


#ifdef OOMPH_HAS_TRIANGLE_LIB
    if (is_a_triangle_mesh_base_mesh)
    {
      triangle_mesh_pt->reset_boundary_element_info(
        ntmp_boundary_elements,
        ntmp_boundary_elements_in_region,
        deleted_f_element_pt);
    } // if (tri_mesh_pt!=0)
    else
    {
#endif // #ifdef OOMPH_HAS_TRIANGLE_LIB

      // Temporarly set the mesh as distributed
      this->setup_boundary_element_info();

#ifdef OOMPH_HAS_TRIANGLE_LIB
    }
#endif

    // Re-setup tree forest. (Call this every time even if
    // a (distributed) mesh has no elements on this processor.
    // We still need to participate in communication.)
    TreeBasedRefineableMeshBase* ref_mesh_pt =
      dynamic_cast<TreeBasedRefineableMeshBase*>(this);
    if (ref_mesh_pt != 0)
    {
      ref_mesh_pt->setup_tree_forest();
    }

    // Classify nodes
    classify_halo_and_haloed_nodes(doc_info, report_stats);

    // Doc?
    //-----
    if (doc_info.is_doc_enabled())
    {
      doc_mesh_distribution(doc_info);
    }
  }


  //========================================================================
  /// (Irreversibly) prune halo(ed) elements and nodes, usually
  /// after another round of refinement, to get rid of
  /// excessively wide halo layers. Note that the current
  /// mesh will be now regarded as the base mesh and no unrefinement
  /// relative to it will be possible once this function
  /// has been called.
  //========================================================================
  void Mesh::prune_halo_elements_and_nodes(
    Vector<GeneralisedElement*>& deleted_element_pt,
    DocInfo& doc_info,
    const bool& report_stats)
  {
    // MemoryUsage::doc_memory_usage("at start of mesh-level prunes");


#ifdef OOMPH_HAS_MPI
    // Delete any external element storage before performing the redistribution
    // (in particular, external halo nodes that are on mesh boundaries)
    this->delete_all_external_storage();
#endif

    TreeBasedRefineableMeshBase* ref_mesh_pt =
      dynamic_cast<TreeBasedRefineableMeshBase*>(this);
    if (ref_mesh_pt != 0)
    {
      // Storage for number of processors and current processor
      int n_proc = Comm_pt->nproc();
      int my_rank = Comm_pt->my_rank();

      // Doc stats
      if (report_stats)
      {
        oomph_info
          << "\n----------------------------------------------------\n";
        oomph_info << "Before pruning: Processor " << my_rank << " holds "
                   << this->nelement() << " elements of which "
                   << this->nroot_halo_element()
                   << " are root halo elements \n while "
                   << this->nroot_haloed_element()
                   << " are root haloed elements" << std::endl;

        // Report total number of halo(ed) and shared nodes for this process
        oomph_info << "Before pruning: Processor " << my_rank << " holds "
                   << this->nnode() << " nodes of which " << this->nhalo_node()
                   << " are halo nodes \n while " << this->nhaloed_node()
                   << " are haloed nodes, and " << this->nshared_node()
                   << " are shared nodes." << std::endl;

        // Report number of halo(ed) and shared nodes with each domain
        // from the current process
        for (int iproc = 0; iproc < n_proc; iproc++)
        {
          oomph_info << "Before pruning: With process " << iproc
                     << ", there are " << this->nhalo_node(iproc)
                     << " halo nodes, and " << std::endl
                     << this->nhaloed_node(iproc) << " haloed nodes, and "
                     << this->nshared_node(iproc) << " shared nodes"
                     << std::endl;
        }
        oomph_info
          << "----------------------------------------------------\n\n";
      }


      double t_start = 0.0;
      if (Global_timings::Doc_comprehensive_timings)
      {
        t_start = TimingHelpers::timer();
      }

      // Declare all nodes as obsolete. We'll
      // change this setting for all nodes that must be retained
      // further down
      unsigned nnod = this->nnode();
      for (unsigned j = 0; j < nnod; j++)
      {
        this->node_pt(j)->set_obsolete();
      }

      // Backup pointers to elements in this mesh (they must be finite elements
      // beacuse it's a refineable mesh)
      unsigned nelem = this->nelement();
      Vector<FiniteElement*> backed_up_el_pt(nelem);
      std::map<RefineableElement*, bool> keep_element;
      for (unsigned e = 0; e < nelem; e++)
      {
        backed_up_el_pt[e] = this->finite_element_pt(e);
      }

      // MemoryUsage::doc_memory_usage("after backed up elements");

      // Get the min and max refinement level, and current refinement pattern
      unsigned min_ref = 0;
      unsigned max_ref = 0;

      // Get min and max refinement level
      unsigned local_min_ref = 0;
      unsigned local_max_ref = 0;
      ref_mesh_pt->get_refinement_levels(local_min_ref, local_max_ref);

      // Reconcile between processors: If (e.g. following distribution/pruning)
      // the mesh has no elements on this processor) then ignore its
      // contribution to the poll of max/min refinement levels
      int int_local_min_ref = local_min_ref;
      int int_local_max_ref = local_max_ref;

      if (nelem == 0)
      {
        int_local_min_ref = INT_MAX;
        int_local_max_ref = INT_MIN;
      }

      int int_min_ref = 0;
      int int_max_ref = 0;

      MPI_Allreduce(&int_local_min_ref,
                    &int_min_ref,
                    1,
                    MPI_INT,
                    MPI_MIN,
                    Comm_pt->mpi_comm());
      MPI_Allreduce(&int_local_max_ref,
                    &int_max_ref,
                    1,
                    MPI_INT,
                    MPI_MAX,
                    Comm_pt->mpi_comm());

      min_ref = unsigned(int_min_ref);
      max_ref = unsigned(int_max_ref);

      // Get refinement pattern
      Vector<Vector<unsigned>> current_refined;
      ref_mesh_pt->get_refinement_pattern(current_refined);

      // get_refinement_pattern refers to the elements at each level
      // that were refined when proceeding to the next level
      int local_n_ref = current_refined.size();

      // Bypass if no elements
      if (nelem == 0)
      {
        local_n_ref = INT_MIN;
      }

      // Reconcile between processors
      int int_n_ref = 0;
      MPI_Allreduce(
        &local_n_ref, &int_n_ref, 1, MPI_INT, MPI_MAX, Comm_pt->mpi_comm());
      unsigned n_ref = int(int_n_ref);


      double t_end = 0.0;
      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info
          << "Time for establishing refinement levels in "
          << " Mesh::prune_halo_elements_and_nodes() [includes comms]: "
          << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

      // MemoryUsage::doc_memory_usage("after establishing refinement levels");

      // Bypass everything until next comms if no elements
      if (nelem > 0)
      {
        // Loop over all elements; keep those on the min refinement level
        // Need to go back to the level indicated by min_ref
        unsigned base_level = n_ref - (max_ref - min_ref);

        // This is the new minimum uniform refinement level
        // relative to total unrefined base mesh
        ref_mesh_pt->uniform_refinement_level_when_pruned() = min_ref;

        // Get the elements at the specified "base" refinement level
        Vector<RefineableElement*> base_level_elements_pt;
        ref_mesh_pt->get_elements_at_refinement_level(base_level,
                                                      base_level_elements_pt);

        // Loop over the elements at this level
        unsigned n_base_el = base_level_elements_pt.size();
        for (unsigned e = 0; e < n_base_el; e++)
        {
          // Extract correct element...
          RefineableElement* ref_el_pt = base_level_elements_pt[e];


          // Check it exists
          if (ref_el_pt != 0)
          {
            // Keep all non-halo elements, remove excess halos
            if ((!ref_el_pt->is_halo()) || (ref_el_pt->must_be_kept_as_halo()))
            {
              keep_element[ref_el_pt] = true;

              // Loop over this non-halo element's nodes and retain them too
              unsigned nnod = ref_el_pt->nnode();
              for (unsigned j = 0; j < nnod; j++)
              {
                ref_el_pt->node_pt(j)->set_non_obsolete();
              }
            }
          }

        } // end loop over base level elements
      }

      // Synchronise refinement level when pruned over all meshes even if they
      // were empty (in which case the uniform refinement level is still zero
      // so go for max
      unsigned n_unif = 0;
      unsigned n_unif_local =
        ref_mesh_pt->uniform_refinement_level_when_pruned();
      MPI_Allreduce(
        &n_unif_local, &n_unif, 1, MPI_UNSIGNED, MPI_MAX, Comm_pt->mpi_comm());
      ref_mesh_pt->uniform_refinement_level_when_pruned() = n_unif;


      t_end = 0.0;
      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info
          << "Time for synchronising refinement levels in "
          << " Mesh::prune_halo_elements_and_nodes() [includes comms]: "
          << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }


      // MemoryUsage::doc_memory_usage("after synchronising refinement levels");


      // Now work on which "root" halo elements to keep at this level
      // Can't use the current set directly; however,
      // we know the refinement level of the current halo, so
      // it is possible to go from that backwards to find the "father
      // halo element" necessary to complete this step

      // Temp map of vectors holding the pointers to the root halo elements
      std::map<unsigned, Vector<FiniteElement*>> tmp_root_halo_element_pt;

      // Temp map of vectors holding the pointers to the root haloed elements
      std::map<unsigned, Vector<FiniteElement*>> tmp_root_haloed_element_pt;

      // Make sure we only push back each element once
      std::map<unsigned, std::map<FiniteElement*, bool>>
        tmp_root_halo_element_is_retained;

      // Map to store if a halo element survives
      std::map<FiniteElement*, bool> halo_element_is_retained;

      for (int domain = 0; domain < n_proc; domain++)
      {
        // Get vector of halo elements with processor domain by copy operation
        Vector<GeneralisedElement*> halo_elem_pt(this->halo_element_pt(domain));

        // Loop over halo elements associated with this adjacent domain
        unsigned nelem = halo_elem_pt.size();
        for (unsigned e = 0; e < nelem; e++)
        {
          // Get element
          RefineableElement* ref_el_pt =
            dynamic_cast<RefineableElement*>(halo_elem_pt[e]);

          // An element should only be kept if its refinement
          // level is the same as the minimum refinement level
          unsigned halo_el_level = ref_el_pt->refinement_level();

          RefineableElement* el_pt = 0;
          if (halo_el_level == min_ref)
          {
            // Already at the correct level
            el_pt = ref_el_pt;
          }
          else
          {
            // Need to go up the tree to the father element at min_ref
            RefineableElement* father_el_pt;
            ref_el_pt->get_father_at_refinement_level(min_ref, father_el_pt);
            el_pt = father_el_pt;
          }

          // Now loop over nodes in the halo element and retain it as
          // halo element if any of it's nodes are shared with one of the
          // non-halo-elements that we've already identified earlier -- these
          // were set to non-obsolete above.
          unsigned nnod = el_pt->nnode();
          for (unsigned j = 0; j < nnod; j++)
          {
            // Node has been reclaimed by one of the non-halo elements above
            Node* nod_pt = el_pt->node_pt(j);
            if (!nod_pt->is_obsolete())
            {
              // Keep element and add it to preliminary storage for
              // halo elements associated with current neighbouring domain
              if (!tmp_root_halo_element_is_retained[domain][el_pt])
              {
                tmp_root_halo_element_pt[domain].push_back(el_pt);
                tmp_root_halo_element_is_retained[domain][el_pt] = true;
              }
              keep_element[el_pt] = true;
              halo_element_is_retained[el_pt] = true;
              break;
            }
          }
        }
      }

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for setup of retention pattern in "
                   << " Mesh::prune_halo_elements_and_nodes(): "
                   << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

      // MemoryUsage::doc_memory_usage("after setup of retention pattern");

      // Make sure everybody finishes this part
      MPI_Barrier(Comm_pt->mpi_comm());

      // Now all processors have decided (independently) which of their
      // (to-be root) halo elements they wish to retain. Now we need to figure
      // out which of their elements are haloed and add them in the appropriate
      // order into the haloed element scheme. For this we exploit that
      // the halo and haloed elements are accessed in the same order on
      // all processors!

      // Identify haloed elements on domain d
      for (int d = 0; d < n_proc; d++)
      {
        // Loop over domains that halo this domain
        for (int dd = 0; dd < n_proc; dd++)
        {
          // Dont't talk to yourself
          if (d != dd)
          {
            // If we're identifying my haloed elements:
            if (d == my_rank)
            {
              // Get vector all elements that are currently haloed by domain dd
              Vector<GeneralisedElement*> haloed_elem_pt(
                this->haloed_element_pt(dd));

              // Create a vector of ints to indicate if the halo element
              // on processor dd processor was kept
              unsigned nelem = haloed_elem_pt.size();
              Vector<int> halo_kept(nelem);

              // Receive this vector from processor dd
              if (nelem != 0)
              {
                MPI_Status status;
                MPI_Recv(&halo_kept[0],
                         nelem,
                         MPI_INT,
                         dd,
                         0,
                         Comm_pt->mpi_comm(),
                         &status);

                // Classify haloed element accordingly
                for (unsigned e = 0; e < nelem; e++)
                {
                  RefineableElement* ref_el_pt =
                    dynamic_cast<RefineableElement*>(haloed_elem_pt[e]);

                  // An element should only be kept if its refinement
                  // level is the same as the minimum refinement level
                  unsigned haloed_el_level = ref_el_pt->refinement_level();

                  // Go up the tree to the correct level
                  RefineableElement* el_pt;

                  if (haloed_el_level == min_ref)
                  {
                    // Already at the correct level
                    el_pt = ref_el_pt;
                  }
                  else
                  {
                    // Need to go up the tree to the father element at min_ref
                    RefineableElement* father_el_pt;
                    ref_el_pt->get_father_at_refinement_level(min_ref,
                                                              father_el_pt);
                    el_pt = father_el_pt;
                  }

                  if (halo_kept[e] == 1)
                  {
                    // I am being haloed by processor dd
                    // Only keep it if it's not already in the storage
                    bool already_root_haloed = false;
                    unsigned n_root_haloed =
                      tmp_root_haloed_element_pt[dd].size();
                    for (unsigned e_root = 0; e_root < n_root_haloed; e_root++)
                    {
                      if (el_pt == tmp_root_haloed_element_pt[dd][e_root])
                      {
                        already_root_haloed = true;
                        break;
                      }
                    }
                    if (!already_root_haloed)
                    {
                      tmp_root_haloed_element_pt[dd].push_back(el_pt);
                    }
                  }
                }
              }
            }
            else
            {
              // If we're dealing with my halo elements:
              if (dd == my_rank)
              {
                // Find (current) halo elements on processor dd whose non-halo
                // is on processor d
                Vector<GeneralisedElement*> halo_elem_pt(
                  this->halo_element_pt(d));

                // Create a vector of ints to indicate if the halo
                // element was kept
                unsigned nelem = halo_elem_pt.size();
                Vector<int> halo_kept(nelem, 0);
                for (unsigned e = 0; e < nelem; e++)
                {
                  RefineableElement* ref_el_pt =
                    dynamic_cast<RefineableElement*>(halo_elem_pt[e]);

                  // An element should only be kept if its refinement
                  // level is the same as the minimum refinement level
                  unsigned halo_el_level = ref_el_pt->refinement_level();

                  // Go up the tree to the correct level
                  RefineableElement* el_pt;
                  if (halo_el_level == min_ref)
                  {
                    // Already at the correct level
                    el_pt = ref_el_pt;
                  }
                  else
                  {
                    // Need to go up the tree to the father element at min_ref
                    RefineableElement* father_el_pt;
                    ref_el_pt->get_father_at_refinement_level(min_ref,
                                                              father_el_pt);
                    el_pt = father_el_pt;
                  }

                  if (halo_element_is_retained[el_pt])
                  {
                    halo_kept[e] = 1;
                  }
                }

                // Now send this vector to processor d to tell it which of
                // the haloed elements (which are listed in the same order)
                // are to be retained as haloed elements.
                if (nelem != 0)
                {
                  MPI_Send(
                    &halo_kept[0], nelem, MPI_INT, d, 0, Comm_pt->mpi_comm());
                }
              }
            }
          }
        }
      }

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for pt2pt comms of retention pattern in "
                   << " Mesh::prune_halo_elements_and_nodes(): "
                   << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

      // MemoryUsage::doc_memory_usage("after pt2pt comms of retention
      // pattern");


      // Backup pointers to nodes in this mesh
      nnod = this->nnode();
      Vector<Node*> backed_up_nod_pt(nnod);
      for (unsigned j = 0; j < nnod; j++)
      {
        backed_up_nod_pt[j] = this->node_pt(j);
      }

      // Flush the mesh storage
      this->flush_element_and_node_storage();

      // Storage for non-active trees/elements that are to be deleted
      std::set<Tree*> trees_to_be_deleted_pt;

      // Loop over all backed-up elements
      nelem = backed_up_el_pt.size();
      for (unsigned e = 0; e < nelem; e++)
      {
        RefineableElement* ref_el_pt =
          dynamic_cast<RefineableElement*>(backed_up_el_pt[e]);

        // Get refinement level
        unsigned level = ref_el_pt->refinement_level();

        // Go up the tree to the correct level
        RefineableElement* el_pt;

        if (level == min_ref)
        {
          // Already at the correct level
          el_pt = ref_el_pt;
        }
        else
        {
          // Need to go up the tree to the father element at min_ref
          RefineableElement* father_el_pt;
          ref_el_pt->get_father_at_refinement_level(min_ref, father_el_pt);
          el_pt = father_el_pt;
        }

        // If the base element is going to be kept, then add the current element
        // to the "new" mesh
        if (keep_element[el_pt])
        {
          this->add_element_pt(ref_el_pt);
        }
        else
        {
          // Flush the object attached to the tree for this element?
          RefineableElement* my_el_pt =
            dynamic_cast<RefineableElement*>(ref_el_pt);
          if (my_el_pt != 0)
          {
            my_el_pt->tree_pt()->flush_object();
          }

          // Get associated tree root
          Tree* tmp_tree_root_pt = my_el_pt->tree_pt()->root_pt();

          // Get all the tree nodes
          Vector<Tree*> tmp_all_tree_nodes_pt;
          tmp_tree_root_pt->stick_all_tree_nodes_into_vector(
            tmp_all_tree_nodes_pt);

          // Loop over all of them and delete associated object/
          // and flush any record of it in tree
          unsigned n_tree = tmp_all_tree_nodes_pt.size();
          for (unsigned j = 0; j < n_tree; j++)
          {
            if (tmp_all_tree_nodes_pt[j]->object_pt() != 0)
            {
              unsigned lev =
                tmp_all_tree_nodes_pt[j]->object_pt()->refinement_level();
              if (lev <= min_ref)
              {
                if (!keep_element[tmp_all_tree_nodes_pt[j]->object_pt()])
                {
                  trees_to_be_deleted_pt.insert(tmp_all_tree_nodes_pt[j]);
                }
              }
            }
          }

          // Delete the element
          deleted_element_pt.push_back(ref_el_pt);
          delete ref_el_pt;
        }
      }

      // MemoryUsage::doc_memory_usage("before deleting superfluous elements");

      // Delete superfluous elements
      for (std::set<Tree*>::iterator it = trees_to_be_deleted_pt.begin();
           it != trees_to_be_deleted_pt.end();
           it++)
      {
        Tree* tree_pt = (*it);
        if (tree_pt->object_pt() != 0)
        {
          deleted_element_pt.push_back(tree_pt->object_pt());
          delete tree_pt->object_pt();
          tree_pt->flush_object();
        }
      }

      // MemoryUsage::doc_memory_usage("after deleting superfluous elements");


      // Wipe the storage scheme for (root) halo(ed) elements and then re-assign
      Root_haloed_element_pt.clear();
      Root_halo_element_pt.clear();
      for (int domain = 0; domain < n_proc; domain++)
      {
        unsigned nelem = tmp_root_halo_element_pt[domain].size();
        for (unsigned e = 0; e < nelem; e++)
        {
          Root_halo_element_pt[domain].push_back(
            tmp_root_halo_element_pt[domain][e]);
        }

        nelem = tmp_root_haloed_element_pt[domain].size();
        for (unsigned e = 0; e < nelem; e++)
        {
          Root_haloed_element_pt[domain].push_back(
            tmp_root_haloed_element_pt[domain][e]);
        }
      }

      //    MemoryUsage::doc_memory_usage(
      //     "after wiping storage scheme for root halo/ed els");

      // Loop over all retained elements at this level and mark their nodes
      //-------------------------------------------------------------------
      // as to be retained too (some double counting going on here)
      //-----------------------------------------------------------
      nelem = this->nelement();
      for (unsigned e = 0; e < nelem; e++)
      {
        FiniteElement* el_pt = this->finite_element_pt(e);

        // Loop over nodes
        unsigned nnod = el_pt->nnode();
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = el_pt->node_pt(j);
          nod_pt->set_non_obsolete();
        }
      }

      // Complete rebuild of mesh by adding retained nodes
      // Note that they are added in the order in which they
      // occured in the original mesh as this guarantees the
      // synchronisity between the serialised access to halo
      // and haloed nodes from different processors.
      nnod = backed_up_nod_pt.size();
      for (unsigned j = 0; j < nnod; j++)
      {
        Node* nod_pt = backed_up_nod_pt[j];
        if (!nod_pt->is_obsolete())
        {
          // Not obsolete so add it back to the mesh
          this->add_node_pt(nod_pt);
        }
      }

      // MemoryUsage::doc_memory_usage("after adding nodes back in");

      // Prune and rebuild mesh
      //-----------------------

      // Now remove the pruned nodes from the boundary lookup scheme
      this->prune_dead_nodes();

      // MemoryUsage::doc_memory_usage("after prune dead nodes");

      // And finally re-setup the boundary lookup scheme for elements
      this->setup_boundary_element_info();

      // MemoryUsage::doc_memory_usage("after setup boundary info");

      // Re-setup tree forest if needed. (Call this every time even if
      // a (distributed) mesh has no elements on this processor.
      // We still need to participate in communication.)
      TreeBasedRefineableMeshBase* ref_mesh_pt =
        dynamic_cast<TreeBasedRefineableMeshBase*>(this);
      if (ref_mesh_pt != 0)
      {
        ref_mesh_pt->setup_tree_forest();
      }

      // MemoryUsage::doc_memory_usage("after setup tree forest");

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info
          << "Time for local rebuild of mesh from retained elements in "
          << " Mesh::prune_halo_elements_and_nodes(): " << t_end - t_start
          << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Classify nodes
      classify_halo_and_haloed_nodes(doc_info, report_stats);

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for Mesh::classify_halo_and_haloed_nodes() "
                   << " Mesh::prune_halo_elements_and_nodes(): "
                   << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

      // MemoryUsage::doc_memory_usage("after
      // classify_halo_and_haloed_nodes()");

      // Doc?
      //-----
      if (doc_info.is_doc_enabled())
      {
        oomph_info << "Outputting distribution in " << doc_info.directory()
                   << " " << doc_info.number() << std::endl;
        doc_mesh_distribution(doc_info);
      }

      // Finally: Reorder the nodes within the mesh's node vector
      // to establish a standard ordering regardless of the sequence
      // of mesh refinements -- this is required to allow dump/restart
      // on refined meshes
      this->reorder_nodes();


      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for Mesh::reorder_nodes() "
                   << " Mesh::prune_halo_elements_and_nodes(): "
                   << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

      // MemoryUsage::doc_memory_usage("after reorder nodes");

      // Doc stats
      if (report_stats)
      {
        oomph_info
          << "\n----------------------------------------------------\n";
        oomph_info << "After pruning:  Processor " << my_rank << " holds "
                   << this->nelement() << " elements of which "
                   << this->nroot_halo_element()
                   << " are root halo elements \n while "
                   << this->nroot_haloed_element()
                   << " are root haloed elements" << std::endl;

        // Report total number of halo(ed) and shared nodes for this process
        oomph_info << "After pruning: Processor " << my_rank << " holds "
                   << this->nnode() << " nodes of which " << this->nhalo_node()
                   << " are halo nodes \n while " << this->nhaloed_node()
                   << " are haloed nodes, and " << this->nshared_node()
                   << " are shared nodes." << std::endl;

        // Report number of halo(ed) and shared nodes with each domain
        // from the current process
        for (int iproc = 0; iproc < n_proc; iproc++)
        {
          oomph_info << "After pruning: With process " << iproc
                     << ", there are " << this->nhalo_node(iproc)
                     << " halo nodes, and " << std::endl
                     << this->nhaloed_node(iproc) << " haloed nodes, and "
                     << this->nshared_node(iproc) << " shared nodes"
                     << std::endl;
        }
        oomph_info
          << "----------------------------------------------------\n\n";
      }
    }

    // MemoryUsage::doc_memory_usage("end of mesh level prune");
  }


  //========================================================================
  ///  Get efficiency of mesh distribution: In an ideal distribution
  /// without halo overhead, each processor would only hold its own
  /// elements. Efficieny per processor =  (number of non-halo elements)/
  /// (total number of elements).
  //========================================================================
  void Mesh::get_efficiency_of_mesh_distribution(double& av_efficiency,
                                                 double& max_efficiency,
                                                 double& min_efficiency)
  {
    // Storage for number of processors and current processor
    int n_proc = Comm_pt->nproc();
    int my_rank = Comm_pt->my_rank();

    // Create vector to hold number of elements and halo elements
    Vector<int> nhalo_elements(n_proc);
    Vector<int> n_elements(n_proc);

    // Count total number of halo elements
    unsigned count = 0;
    for (int d = 0; d < n_proc; d++)
    {
      Vector<GeneralisedElement*> halo_elem_pt(halo_element_pt(d));
      count += halo_elem_pt.size();
    }

    // Stick own number into appropriate entry
    int nhalo_element_local = count;
    int n_element_local = nelement();

    // Gather information on root processor: First argument group
    // specifies what is to be sent (one int from each procssor, indicating
    // the number of elements on it), the second group indicates where
    // the results are to be gathered (in rank order) on root processor.
    MPI_Gather(&nhalo_element_local,
               1,
               MPI_INT,
               &nhalo_elements[0],
               1,
               MPI_INT,
               0,
               Comm_pt->mpi_comm());
    MPI_Gather(&n_element_local,
               1,
               MPI_INT,
               &n_elements[0],
               1,
               MPI_INT,
               0,
               Comm_pt->mpi_comm());

    // Initialise stats
    av_efficiency = 0.0;
    double max = -1.0;
    double min = 1000000000.0;

    if (my_rank == 0)
    {
      for (int i = 0; i < n_proc; i++)
      {
        unsigned nel = n_elements[i];
        double eff = 0.0;
        if (nel != 0)
        {
          eff = double(n_elements[i] - nhalo_elements[i]) / double(nel);
        }
        av_efficiency += eff;
        if (eff > max) max = eff;
        if (eff < min) min = eff;
      }
      av_efficiency /= double(n_proc);
    }

    // Now broadcast the result back out
    MPI_Bcast(&max, 1, MPI_DOUBLE, 0, Comm_pt->mpi_comm());
    MPI_Bcast(&min, 1, MPI_DOUBLE, 0, Comm_pt->mpi_comm());
    MPI_Bcast(&av_efficiency, 1, MPI_DOUBLE, 0, Comm_pt->mpi_comm());

    max_efficiency = max;
    min_efficiency = min;
  }


  //========================================================================
  /// Doc the mesh distribution
  //========================================================================
  void Mesh::doc_mesh_distribution(DocInfo& doc_info)
  {
    // Storage for current processor and number of processors
    int my_rank = Comm_pt->my_rank();
    int n_proc = Comm_pt->nproc();

    std::ostringstream filename;
    std::ostringstream filename2;
    std::ofstream some_file;
    std::ofstream some_file2;

    // Doc elements on this processor
    filename << doc_info.directory() << "/" << doc_info.label()
             << "elements_on_proc" << my_rank << "_" << doc_info.number()
             << ".dat";
    some_file.open(filename.str().c_str());
    this->output(some_file, 5);
    some_file.close();

    // Doc non-halo elements on this processor
    filename.str("");
    filename << doc_info.directory() << "/" << doc_info.label()
             << "non_halo_elements_on_proc" << my_rank << "_"
             << doc_info.number() << ".dat";
    some_file.open(filename.str().c_str());

    // Get to elements on processor
    unsigned nelem = this->nelement();
    for (unsigned e = 0; e < nelem; e++)
    {
      // Get the element
      GeneralisedElement* el_pt = this->element_pt(e);
      // Only output if not  a halo element
      if (!el_pt->is_halo())
      {
        FiniteElement* f_el_pt = dynamic_cast<FiniteElement*>(el_pt);
        // Indicate a generalised element
        if (f_el_pt == 0)
        {
          some_file << "Generalised Element " << e << "\n";
        }
        else
        // output if non-halo and a finite element
        {
          f_el_pt->output(some_file, 5);
        }
      }
    }

    some_file.close();


    // Doc halo elements on this processor
    filename.str("");
    filename << doc_info.directory() << "/" << doc_info.label()
             << "halo_elements_on_proc" << my_rank << "_" << doc_info.number()
             << ".dat";
    some_file.open(filename.str().c_str());
    for (int domain = 0; domain < n_proc; domain++)
    {
      filename2.str("");
      filename2 << doc_info.directory() << "/" << doc_info.label()
                << "halo_elements_with_proc" << domain << "_on_proc" << my_rank
                << "_" << doc_info.number() << ".dat";
      some_file2.open(filename2.str().c_str());

      // Get vector of halo elements by copy operation
      Vector<GeneralisedElement*> halo_elem_pt(this->halo_element_pt(domain));
      unsigned nelem = halo_elem_pt.size();
      for (unsigned e = 0; e < nelem; e++)
      {
        FiniteElement* f_el_pt = dynamic_cast<FiniteElement*>(halo_elem_pt[e]);
        if (f_el_pt != 0)
        {
#ifdef PARANOID
          // Check if it's refineable and if so if it's a leaf.
          // If not something must have gone wrong.
          RefineableElement* ref_el_pt =
            dynamic_cast<RefineableElement*>(f_el_pt);
          if (ref_el_pt != 0)
          {
            if (!(ref_el_pt->tree_pt()->is_leaf()))
            {
              std::ostringstream error_message;
              error_message
                << "Haloed element is not a leaf element. This shouldn't happen"
                << std::endl;
              error_message << "Here are the nodal positions: " << std::endl;
              unsigned nnod = ref_el_pt->nnode();
              for (unsigned j = 0; j < nnod; j++)
              {
                Node* nod_pt = ref_el_pt->node_pt(j);
                unsigned n_dim = nod_pt->ndim();
                for (unsigned i = 0; i < n_dim; i++)
                {
                  error_message << nod_pt->x(i) << " ";
                }
                error_message << "\n";
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
            }
          }
#endif

          f_el_pt->output(some_file, 5);
          f_el_pt->output(some_file2, 5);
        }
        // Indicate a generalised element
        else
        {
          some_file << "Generalised Element " << e << "\n";
          some_file2 << "Generalised Element " << e << "\n";
        }
      }
      some_file2.close();
    }
    some_file.close();

    // Doc haloed elements on this processor
    filename.str("");
    filename << doc_info.directory() << "/" << doc_info.label()
             << "haloed_elements_on_proc" << my_rank << "_" << doc_info.number()
             << ".dat";
    some_file.open(filename.str().c_str());
    for (int domain = 0; domain < n_proc; domain++)
    {
      filename2.str("");
      filename2 << doc_info.directory() << "/" << doc_info.label()
                << "haloed_elements_with_proc" << domain << "_on_proc"
                << my_rank << "_" << doc_info.number() << ".dat";
      some_file2.open(filename2.str().c_str());

      // Get vector of haloed elements by copy operation
      Vector<GeneralisedElement*> haloed_elem_pt(
        this->haloed_element_pt(domain));
      unsigned nelem = haloed_elem_pt.size();
      for (unsigned e = 0; e < nelem; e++)
      {
        // Is it a finite element
        FiniteElement* finite_el_pt =
          dynamic_cast<FiniteElement*>(haloed_elem_pt[e]);
        if (finite_el_pt != 0)
        {
#ifdef PARANOID
          // Check if it's refineable and if so if it's a leaf.
          // If not something must have gone wrong.
          RefineableElement* ref_el_pt =
            dynamic_cast<RefineableElement*>(finite_el_pt);
          if (ref_el_pt != 0)
          {
            if (!(ref_el_pt->tree_pt()->is_leaf()))
            {
              std::ostringstream error_message;
              error_message
                << "Haloed element is not a leaf element. This shouldn't happen"
                << std::endl;
              error_message << "Here are the nodal positions: " << std::endl;
              unsigned nnod = ref_el_pt->nnode();
              for (unsigned j = 0; j < nnod; j++)
              {
                Node* nod_pt = ref_el_pt->node_pt(j);
                unsigned n_dim = nod_pt->ndim();
                for (unsigned i = 0; i < n_dim; i++)
                {
                  error_message << nod_pt->x(i) << " ";
                }
                error_message << "\n";
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
            }
          }
#endif

          finite_el_pt->output(some_file, 5);
          finite_el_pt->output(some_file2, 5);
        }
        // Indicate a generalised element
        else
        {
          some_file << "Generalised  Element " << e << "\n";
          some_file2 << "Generalised  Element " << e << "\n";
        }
      }
      some_file2.close();
    }
    some_file.close();


    // Doc non-halo nodes on this processor
    filename.str("");
    filename << doc_info.directory() << "/" << doc_info.label()
             << "non_halo_nodes_on_proc" << my_rank << "_" << doc_info.number()
             << ".dat";
    some_file.open(filename.str().c_str());
    unsigned nnod = this->nnode();
    for (unsigned j = 0; j < nnod; j++)
    {
      Node* nod_pt = this->node_pt(j);
      if (!nod_pt->is_halo())
      {
        unsigned ndim = nod_pt->ndim();
        for (unsigned i = 0; i < ndim; i++)
        {
          some_file << nod_pt->x(i) << " ";
        }
        some_file << nod_pt->nvalue() << " " << nod_pt->non_halo_proc_ID()
                  << " " << nod_pt->is_halo() << " " << nod_pt->hang_code()
                  << std::endl;
      }
    }
    some_file.close();


    // Doc nodes on this processor
    filename.str("");
    filename << doc_info.directory() << "/" << doc_info.label()
             << "nodes_on_proc" << my_rank << "_" << doc_info.number()
             << ".dat";
    some_file.open(filename.str().c_str());
    for (unsigned j = 0; j < nnod; j++)
    {
      Node* nod_pt = this->node_pt(j);
      unsigned ndim = nod_pt->ndim();
      for (unsigned i = 0; i < ndim; i++)
      {
        some_file << nod_pt->x(i) << " ";
      }
      some_file << nod_pt->nvalue() << " " << nod_pt->non_halo_proc_ID() << " "
                << nod_pt->is_halo() << " " << nod_pt->hang_code() << std::endl;
    }
    some_file.close();

    // Doc solid nodes on this processor
    filename.str("");
    filename << doc_info.directory() << "/" << doc_info.label()
             << "solid_nodes_on_proc" << my_rank << "_" << doc_info.number()
             << ".dat";
    some_file.open(filename.str().c_str());
    unsigned nsnod = this->nnode();
    for (unsigned j = 0; j < nsnod; j++)
    {
      Node* nod_pt = this->node_pt(j);
      SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
      if (solid_nod_pt != 0)
      {
        unsigned ndim = solid_nod_pt->ndim();
        for (unsigned i = 0; i < ndim; i++)
        {
          some_file << nod_pt->x(i) << " ";
        }
        some_file << nod_pt->nvalue() << " " << nod_pt->non_halo_proc_ID()
                  << " " << nod_pt->is_halo() << " " << nod_pt->hang_code()
                  << "\n";
      }
    }
    some_file.close();


    // Doc halo nodes on this processor
    filename.str("");
    filename << doc_info.directory() << "/" << doc_info.label()
             << "halo_nodes_on_proc" << my_rank << "_" << doc_info.number()
             << ".dat";
    some_file.open(filename.str().c_str());
    for (int domain = 0; domain < n_proc; domain++)
    {
      filename2.str("");
      filename2 << doc_info.directory() << "/" << doc_info.label()
                << "halo_nodes_with_proc" << domain << "_on_proc" << my_rank
                << "_" << doc_info.number() << ".dat";
      some_file2.open(filename2.str().c_str());
      unsigned nnod = this->nhalo_node(domain);
      for (unsigned j = 0; j < nnod; j++)
      {
        Node* nod_pt = this->halo_node_pt(domain, j);
        unsigned ndim = nod_pt->ndim();
        for (unsigned i = 0; i < ndim; i++)
        {
          some_file << nod_pt->x(i) << " ";
          some_file2 << nod_pt->x(i) << " ";
        }
        some_file << nod_pt->nvalue() << " " << nod_pt->non_halo_proc_ID()
                  << " " << nod_pt->is_halo() << " " << nod_pt->hang_code()
                  << "\n";
        some_file2 << nod_pt->nvalue() << " " << nod_pt->non_halo_proc_ID()
                   << " " << nod_pt->is_halo() << " " << nod_pt->hang_code()
                   << "\n";
      }
      some_file2.close();
    }
    some_file.close();


    // Doc haloed nodes on this processor
    filename.str("");
    filename << doc_info.directory() << "/" << doc_info.label()
             << "haloed_nodes_on_proc" << my_rank << "_" << doc_info.number()
             << ".dat";
    some_file.open(filename.str().c_str());
    for (int domain = 0; domain < n_proc; domain++)
    {
      filename2.str("");
      filename2 << doc_info.directory() << "/" << doc_info.label()
                << "haloed_nodes_with_proc" << domain << "_on_proc" << my_rank
                << "_" << doc_info.number() << ".dat";
      some_file2.open(filename2.str().c_str());
      unsigned nnod = this->nhaloed_node(domain);
      for (unsigned j = 0; j < nnod; j++)
      {
        Node* nod_pt = this->haloed_node_pt(domain, j);
        unsigned ndim = nod_pt->ndim();
        for (unsigned i = 0; i < ndim; i++)
        {
          some_file << nod_pt->x(i) << " ";
          some_file2 << nod_pt->x(i) << " ";
        }
        some_file << nod_pt->nvalue() << " " << nod_pt->non_halo_proc_ID()
                  << " " << nod_pt->is_halo() << " " << nod_pt->hang_code()
                  << "\n";
        some_file2 << nod_pt->nvalue() << " " << nod_pt->non_halo_proc_ID()
                   << " " << nod_pt->is_halo() << " " << nod_pt->hang_code()
                   << "\n";
      }
      some_file2.close();
    }
    some_file.close();


    // Doc shared nodes on this processor
    filename.str("");
    filename << doc_info.directory() << "/" << doc_info.label()
             << "shared_nodes_on_proc" << my_rank << "_" << doc_info.number()
             << ".dat";
    some_file.open(filename.str().c_str());
    for (int domain = 0; domain < n_proc; domain++)
    {
      filename2.str("");
      filename2 << doc_info.directory() << "/" << doc_info.label()
                << "shared_nodes_with_proc" << domain << "_on_proc" << my_rank
                << "_" << doc_info.number() << ".dat";
      some_file2.open(filename2.str().c_str());
      unsigned nnod = this->nshared_node(domain);
      for (unsigned j = 0; j < nnod; j++)
      {
        Node* nod_pt = this->shared_node_pt(domain, j);
        unsigned ndim = nod_pt->ndim();
        for (unsigned i = 0; i < ndim; i++)
        {
          some_file << nod_pt->x(i) << " ";
          some_file2 << nod_pt->x(i) << " ";
        }
        some_file << nod_pt->nvalue() << " " << nod_pt->non_halo_proc_ID()
                  << " " << nod_pt->is_halo() << "\n";
        some_file2 << nod_pt->nvalue() << " " << nod_pt->non_halo_proc_ID()
                   << " " << nod_pt->is_halo() << " " << nod_pt->hang_code()
                   << "\n";
      }
      some_file2.close();
    }
    some_file.close();


    // Doc mesh
    filename.str("");
    filename << doc_info.directory() << "/" << doc_info.label() << "mesh"
             << my_rank << "_" << doc_info.number() << ".dat";
    some_file.open(filename.str().c_str());
    this->output(some_file, 5);
    some_file.close();


    // Doc boundary scheme
    filename.str("");
    filename << doc_info.directory() << "/" << doc_info.label() << "boundaries"
             << my_rank << "_" << doc_info.number() << ".dat";
    some_file.open(filename.str().c_str());
    this->output_boundaries(some_file);
    some_file.close();


    // Doc elements next to boundaries scheme
    // if set up
    if (Lookup_for_elements_next_boundary_is_setup)
    {
      // How many finite elements are adjacent to boundary b?
      unsigned nbound = this->nboundary();
      for (unsigned b = 0; b < nbound; b++)
      {
        filename.str("");
        filename << doc_info.directory() << "/" << doc_info.label()
                 << "boundary_elements" << my_rank << "_" << b << "_"
                 << doc_info.number() << ".dat";
        some_file.open(filename.str().c_str());
        unsigned nelem = this->nboundary_element(b);
        for (unsigned e = 0; e < nelem; e++)
        {
          this->boundary_element_pt(b, e)->output(some_file, 5);
        }
        some_file.close();
      }
    }
  }


  //========================================================================
  /// Check the halo/haloed/shared node/element schemes on the Mesh
  //========================================================================
  void Mesh::check_halo_schemes(DocInfo& doc_info,
                                double& max_permitted_error_for_halo_check)
  {
    oomph_info << "Doing check_halo_schemes for mesh: " << typeid(*this).name()
               << std::endl;

    // Moved this from the Problem class so that it would work better
    // in multiple mesh problems; there remains a simple "wrapper"
    // function in the Problem class that calls this for each (sub)mesh.

    MPI_Status status;
    std::ostringstream filename;
    std::ofstream shared_file;
    std::ofstream halo_file;
    std::ofstream haloed_file;
    std::ofstream ext_halo_file;
    std::ofstream ext_haloed_file;

    // Storage for current processor and number of processors
    int n_proc = Comm_pt->nproc();
    int my_rank = Comm_pt->my_rank();


    // Check the shared node scheme first: if this is incorrect then
    // the halo(ed) node scheme is likely to be wrong too

    // Doc shared nodes lookup schemes
    //-------------------------------------
    if (doc_info.is_doc_enabled())
    {
      // Loop over domains for shared nodes
      for (int dd = 0; dd < n_proc; dd++)
      {
        filename.str("");
        filename << doc_info.directory() << "/shared_node_check"
                 << doc_info.label() << "_on_proc" << my_rank << "_with_proc"
                 << dd << "_" << doc_info.number() << ".dat";
        shared_file.open(filename.str().c_str());
        shared_file << "ZONE " << std::endl;

        unsigned nnod = nshared_node(dd);
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = shared_node_pt(dd, j);
          unsigned ndim = nod_pt->ndim();
          for (unsigned i = 0; i < ndim; i++)
          {
            shared_file << nod_pt->position(i) << " ";
          }
          shared_file << j << " " << nod_pt << std::endl;
        }
        // Dummy output for processor that doesn't share nodes
        // (needed for tecplot)
        if ((nnod == 0) && (nelement() != 0))
        {
          FiniteElement* f_el_pt = dynamic_cast<FiniteElement*>(element_pt(0));
          // If it's a generalised element mesh dummy output
          if (f_el_pt == 0)
          {
            shared_file << "0.0" << std::endl;
          }
          else
          {
            unsigned ndim = f_el_pt->node_pt(0)->ndim();
            if (ndim == 2)
            {
              shared_file << " 1.0 1.1 " << std::endl;
            }
            else
            {
              shared_file << " 1.0 1.1 1.1" << std::endl;
            }
          }
        }
        shared_file.close();
      }
    }


    // Check for duplicates in shared node scheme
    for (int d = 0; d < n_proc; d++)
    {
      unsigned n_vector = Shared_node_pt[d].size();
      std::set<Node*> shared_node_set;
      for (unsigned i = 0; i < n_vector; i++)
      {
        unsigned old_size = shared_node_set.size();
        shared_node_set.insert(Shared_node_pt[d][i]);
        unsigned new_size = shared_node_set.size();
        if (old_size == new_size)
        {
          Node* nod_pt = Shared_node_pt[d][i];
          oomph_info << "Repeated node in shared node lookup scheme: " << i
                     << "-th node with proc " << d << " : " << nod_pt
                     << " at: ";
          unsigned n = nod_pt->ndim();
          for (unsigned ii = 0; ii < n; ii++)
          {
            oomph_info << nod_pt->x(ii) << " ";
          }
          oomph_info << std::endl;
        }
      }
      unsigned n_set = shared_node_set.size();
      if (n_vector != n_set)
      {
        std::ostringstream warning_stream;
        warning_stream
          << "WARNING: " << std::endl
          << "There seem to be repeated entries in shared node scheme "
          << "with proc " << d << ".\n"
          << "n_vector=" << n_vector << "; n_set=" << n_set << std::endl;
        if (doc_info.is_doc_enabled())
        {
          warning_stream << "Shared node scheme has been output in files like\n"
                         << filename.str() << std::endl;
        }
        else
        {
          warning_stream << "Re-run with doc_info enabled to see where the "
                            "shared nodes are.\n";
        }
        OomphLibWarning(warning_stream.str(),
                        "Mesh::check_halo_schemes",
                        OOMPH_EXCEPTION_LOCATION);
      }
    }


    // Check shared nodes lookup schemes
    //----------------------------------
    double max_error = 0.0;

    // Loop over domains for shared nodes
    for (int d = 0; d < n_proc; d++)
    {
      // Are my shared nodes being checked?
      if (d == my_rank)
      {
        // Loop over domains for shared nodes
        for (int dd = 0; dd < n_proc; dd++)
        {
          // Don't talk to yourself
          if (dd != d)
          {
            // How many of my nodes are shared nodes with processor dd?
            int nnod_shared = nshared_node(dd);

            if (nnod_shared != 0)
            {
              // Receive from processor dd how many of his nodes are shared
              // with this processor
              int nnod_share = 0;
              MPI_Recv(
                &nnod_share, 1, MPI_INT, dd, 0, Comm_pt->mpi_comm(), &status);

              if (nnod_shared != nnod_share)
              {
                std::ostringstream error_message;

                error_message << "Clash in numbers of shared nodes! "
                              << std::endl;
                error_message << "# of shared nodes on proc " << dd << ": "
                              << nnod_shared << std::endl;
                error_message << "# of shared nodes on proc " << d << ": "
                              << nnod_share << std::endl;
                error_message << "(Re-)run Problem::check_halo_schemes() with "
                                 "DocInfo object"
                              << std::endl;
                error_message << "to identify the problem" << std::endl;
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }


              unsigned nod_dim = finite_element_pt(0)->node_pt(0)->ndim();

              // Get strung-together nodal positions from other processor
              Vector<double> other_nodal_positions(nod_dim * nnod_share);
              MPI_Recv(&other_nodal_positions[0],
                       nod_dim * nnod_share,
                       MPI_DOUBLE,
                       dd,
                       0,
                       Comm_pt->mpi_comm(),
                       &status);

              // Check
              unsigned count = 0;
              for (int j = 0; j < nnod_share; j++)
              {
                Vector<double> x_shared(nod_dim);
                for (unsigned i = 0; i < nod_dim; i++)
                {
                  x_shared[i] = shared_node_pt(dd, j)->position(i);
                }
                Vector<double> x_share(nod_dim);
                for (unsigned i = 0; i < nod_dim; i++)
                {
                  x_share[i] = other_nodal_positions[count];
                  count++;
                }
                double error = 0.0;
                for (unsigned i = 0; i < nod_dim; i++)
                {
                  error +=
                    (x_shared[i] - x_share[i]) * (x_shared[i] - x_share[i]);
                }
                error = sqrt(error);

                // Doc if relevant
                if (error > max_permitted_error_for_halo_check)
                {
                  oomph_info << "Error in shared nodes: ";
                  for (unsigned i = 0; i < nod_dim; i++)
                  {
                    oomph_info << x_shared[i] << " ";
                  }
                  for (unsigned i = 0; i < nod_dim; i++)
                  {
                    oomph_info << x_share[i] << " ";
                  }
                  oomph_info << error << std::endl;
                  oomph_info << "shared node: " << shared_node_pt(dd, j)
                             << std::endl;
                  if (shared_node_pt(dd, j)->is_hanging())
                  {
                    oomph_info << "shared node: " << shared_node_pt(dd, j)
                               << " is hanging with masters" << std::endl;
                    for (unsigned m = 0;
                         m < shared_node_pt(dd, j)->hanging_pt()->nmaster();
                         m++)
                    {
                      oomph_info
                        << "master: "
                        << shared_node_pt(dd, j)->hanging_pt()->master_node_pt(
                             m)
                        << std::endl;
                    }
                  }
                }

                // Keep tracking max
                if (error > max_error)
                {
                  max_error = error;
                }
              }
            }
          }
        }
      }
      // My shared nodes are not being checked: Send my shared nodes
      // to the other processor
      else
      {
        int nnod_share = nshared_node(d);

        if (nnod_share != 0)
        {
          // Send it across to the processor whose shared nodes are being
          // checked
          MPI_Send(&nnod_share, 1, MPI_INT, d, 0, Comm_pt->mpi_comm());

          unsigned nod_dim = finite_element_pt(0)->node_pt(0)->ndim();

          // Now string together the nodal positions of all shared nodes
          Vector<double> nodal_positions(nod_dim * nnod_share);
          unsigned count = 0;
          for (int j = 0; j < nnod_share; j++)
          {
            for (unsigned i = 0; i < nod_dim; i++)
            {
              nodal_positions[count] = shared_node_pt(d, j)->position(i);
              count++;
            }
          }
          // Send it across to the processor whose shared nodes are being
          // checked
          MPI_Send(&nodal_positions[0],
                   nod_dim * nnod_share,
                   MPI_DOUBLE,
                   d,
                   0,
                   Comm_pt->mpi_comm());
        }
      }
    }

    oomph_info << "Max. error for shared nodes " << max_error << std::endl;
    if (max_error > max_permitted_error_for_halo_check)
    {
      std::ostringstream error_message;
      error_message << "This is bigger than the permitted threshold "
                    << max_permitted_error_for_halo_check << std::endl;
      error_message
        << "If you believe this to be acceptable for your problem\n"
        << "increase Problem::Max_permitted_error_for_halo_check and re-run \n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Now check the halo/haloed element lookup scheme

    // Doc halo/haoloed element lookup schemes
    //-----------------------------------------
    if (doc_info.is_doc_enabled())
    {
      // Loop over domains for halo elements
      for (int dd = 0; dd < n_proc; dd++)
      {
        filename.str("");
        filename << doc_info.directory() << "/halo_element_check"
                 << doc_info.label() << "_on_proc" << my_rank << "_with_proc"
                 << dd << "_" << doc_info.number() << ".dat";
        halo_file.open(filename.str().c_str());

        // Get vectors of halo/haloed elements by copy operation
        Vector<GeneralisedElement*> halo_elem_pt(halo_element_pt(dd));

        unsigned nelem = halo_elem_pt.size();

        for (unsigned e = 0; e < nelem; e++)
        {
          halo_file << "ZONE " << std::endl;
          // Can I cast to a finite element
          FiniteElement* finite_el_pt =
            dynamic_cast<FiniteElement*>(halo_elem_pt[e]);
          if (finite_el_pt != 0)
          {
            unsigned nnod = finite_el_pt->nnode();
            for (unsigned j = 0; j < nnod; j++)
            {
              Node* nod_pt = finite_el_pt->node_pt(j);
              unsigned ndim = nod_pt->ndim();
              for (unsigned i = 0; i < ndim; i++)
              {
                halo_file << nod_pt->position(i) << " ";
              }
              halo_file << std::endl;
            }
          }
        }
        halo_file.close();
      }

      // Loop over domains for halo elements
      for (int d = 0; d < n_proc; d++)
      {
        filename.str("");
        filename << doc_info.directory() << "/haloed_element_check"
                 << doc_info.label() << "_on_proc" << my_rank << "_with_proc"
                 << d << "_" << doc_info.number() << ".dat";
        haloed_file.open(filename.str().c_str());

        // Get vectors of halo/haloed elements by copy operation
        Vector<GeneralisedElement*> haloed_elem_pt(haloed_element_pt(d));

        unsigned nelem2 = haloed_elem_pt.size();
        for (unsigned e = 0; e < nelem2; e++)
        {
          haloed_file << "ZONE " << std::endl;
          // Is it a finite element
          FiniteElement* finite_el_pt =
            dynamic_cast<FiniteElement*>(haloed_elem_pt[e]);
          if (finite_el_pt != 0)
          {
            // Check if it's refineable and if so if it's a leaf.
            // If not something must have gone wrong.
            RefineableElement* ref_el_pt =
              dynamic_cast<RefineableElement*>(finite_el_pt);
            if (ref_el_pt != 0)
            {
              if (!(ref_el_pt->tree_pt()->is_leaf()))
              {
                std::ostringstream error_message;
                error_message << "Haloed element is not a leaf element. This "
                                 "shouldn't happen"
                              << std::endl;
                error_message << "Here are the nodal positions: " << std::endl;
                unsigned nnod = ref_el_pt->nnode();
                for (unsigned j = 0; j < nnod; j++)
                {
                  Node* nod_pt = ref_el_pt->node_pt(j);
                  unsigned n_dim = nod_pt->ndim();
                  for (unsigned i = 0; i < n_dim; i++)
                  {
                    error_message << nod_pt->x(i) << " ";
                  }
                  error_message << "\n";
                  throw OomphLibError(error_message.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
              }
            }
            unsigned nnod2 = finite_el_pt->nnode();
            for (unsigned j = 0; j < nnod2; j++)
            {
              Node* nod_pt = finite_el_pt->node_pt(j);
              unsigned ndim = nod_pt->ndim();
              for (unsigned i = 0; i < ndim; i++)
              {
                haloed_file << nod_pt->position(i) << " ";
              }
              haloed_file << std::endl;
            }
          }
        }
        haloed_file.close();
      }
    } // end of if doc flag

    // Check halo/haloed element lookup schemes
    //-----------------------------------------
    max_error = 0.0;
    bool shout = false;
    bool shout_and_terminate = false;

    // Loop over domains for haloed elements
    for (int d = 0; d < n_proc; d++)
    {
      // Are my haloed elements being checked?
      if (d == my_rank)
      {
        // Loop over domains for halo elements
        for (int dd = 0; dd < n_proc; dd++)
        {
          // Don't talk to yourself
          if (dd != d)
          {
            // Get vectors of haloed elements by copy operation
            Vector<GeneralisedElement*> haloed_elem_pt(haloed_element_pt(dd));

            // How many of my elements are haloed elements whose halo
            // counterpart is located on processor dd?
            int nelem_haloed = haloed_elem_pt.size();

            if (nelem_haloed != 0)
            {
              // Receive from processor dd how many of his elements are halo
              // nodes whose non-halo counterparts are located here
              int nelem_halo = 0;
              MPI_Recv(
                &nelem_halo, 1, MPI_INT, dd, 0, Comm_pt->mpi_comm(), &status);
              if (nelem_halo != nelem_haloed)
              {
                std::ostringstream error_message;
                error_message
                  << "Clash in numbers of halo and haloed elements! "
                  << std::endl;
                error_message << "# of haloed elements whose halo counterpart "
                                 "lives on proc "
                              << dd << ": " << nelem_haloed << std::endl;
                error_message << "# of halo elements whose non-halo "
                                 "counterpart lives on proc "
                              << d << ": " << nelem_halo << std::endl;
                error_message << "(Re-)run Problem::check_halo_schemes() with "
                                 "DocInfo object"
                              << std::endl;
                error_message << "to identify the problem" << std::endl;
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }


              // We can only check nodal stuff for meshes of finite elements
              if (dynamic_cast<FiniteElement*>(this->element_pt(0)))
              {
                // Get strung-together elemental nodal positions
                // from other processor
                // unsigned nnod_per_el=finite_element_pt(0)->nnode();
                unsigned nod_dim = finite_element_pt(0)->node_pt(0)->ndim();
                // Vector<double> other_nodal_positions
                // (nod_dim*nnod_per_el*nelem_halo);
                // MPI_Recv(&other_nodal_positions[0],nod_dim*nnod_per_el*nelem_halo,
                //         MPI_DOUBLE,dd,0,comm_pt->mpi_comm(),&status);
                unsigned n_nodal_positions = 0;
                MPI_Recv(&n_nodal_positions,
                         1,
                         MPI_UNSIGNED,
                         dd,
                         0,
                         Comm_pt->mpi_comm(),
                         &status);
                Vector<double> other_nodal_positions(n_nodal_positions);
                if (n_nodal_positions > 0)
                {
                  MPI_Recv(&other_nodal_positions[0],
                           n_nodal_positions,
                           MPI_DOUBLE,
                           dd,
                           0,
                           Comm_pt->mpi_comm(),
                           &status);
                }


                // Receive hanging info to be checked
                Vector<int> other_nodal_hangings;
                unsigned n_other_nodal_hangings = 0;
                MPI_Recv(&n_other_nodal_hangings,
                         1,
                         MPI_UNSIGNED,
                         dd,
                         1,
                         Comm_pt->mpi_comm(),
                         &status);
                if (n_other_nodal_hangings > 0)
                {
                  other_nodal_hangings.resize(n_other_nodal_hangings);
                  MPI_Recv(&other_nodal_hangings[0],
                           n_other_nodal_hangings,
                           MPI_INT,
                           dd,
                           1,
                           Comm_pt->mpi_comm(),
                           &status);
                }

                // If documenting, open the output files
                if (doc_info.is_doc_enabled())
                {
                  filename.str("");
                  filename << doc_info.directory() << "/error_haloed_check"
                           << doc_info.label() << "_on_proc" << my_rank
                           << "_with_proc" << dd << "_" << doc_info.number()
                           << ".dat";
                  haloed_file.open(filename.str().c_str());
                  filename.str("");
                  filename << doc_info.directory() << "/error_halo_check"
                           << doc_info.label() << "_on_proc" << my_rank
                           << "_with_proc" << dd << "_" << doc_info.number()
                           << ".dat";
                  halo_file.open(filename.str().c_str());
                }

                unsigned count = 0;
                unsigned count_hanging = 0;
                for (int e = 0; e < nelem_haloed; e++)
                {
                  FiniteElement* finite_el_pt =
                    dynamic_cast<FiniteElement*>(haloed_elem_pt[e]);

                  if (finite_el_pt != 0)
                  {
                    unsigned nnod_this_el = finite_el_pt->nnode();
                    for (unsigned j = 0; j < nnod_this_el; j++)
                    {
                      Node* nod_pt = finite_el_pt->node_pt(j);
                      // unsigned nod_dim = nod_pt->ndim();

                      // Testing POSITIONS, not x location
                      // (cf hanging nodes, nodes.h)
                      Vector<double> x_haloed(nod_dim);
                      for (unsigned i = 0; i < nod_dim; i++)
                      {
                        x_haloed[i] = nod_pt->position(i);
                      }

                      Vector<double> x_halo(nod_dim);
                      for (unsigned i = 0; i < nod_dim; i++)
                      {
                        x_halo[i] = other_nodal_positions[count];
                        ++count;
                      }

                      double error = 0.0;
                      shout = false;
                      for (unsigned i = 0; i < nod_dim; i++)
                      {
                        error +=
                          (x_haloed[i] - x_halo[i]) * (x_haloed[i] - x_halo[i]);
                      }
                      error = sqrt(error);

                      if (error > max_error)
                      {
                        max_error = error;
                      }
                      double tol = 1.0e-12;
                      if (error > tol)
                      {
                        oomph_info
                          << "Discrepancy between nodal coordinates of halo(ed)"
                          << "element larger than tolerance (" << tol
                          << ")\n  Error: " << error << "\n";
                        shout = true;
                      }

                      unsigned nval = nod_pt->nvalue();
                      int nval_other = other_nodal_hangings[count_hanging];
                      count_hanging++;
                      if (int(nval) != nval_other)
                      {
                        oomph_info
                          << "Number of values of node, " << nval
                          << ", does not match number of values on other proc, "
                          << nval_other << std::endl;
                        shout = true;
                        shout_and_terminate = true;
                      }

                      // Is other node geometrically hanging?
                      int other_geom_hanging = 0;

                      // Check hangingness/number of master nodes
                      for (int i = -1; i < int(nval); i++)
                      {
                        int nmaster_other = other_nodal_hangings[count_hanging];
                        count_hanging++;

                        // Record geom hang status of other node
                        if (i == -1) other_geom_hanging = nmaster_other;

                        // Value is hanging on local proc: Does it have the same
                        // number of masters as its counterpart on other proc?
                        if (nod_pt->is_hanging(i))
                        {
                          unsigned nmaster = nod_pt->hanging_pt(i)->nmaster();
                          if (int(nmaster) != nmaster_other)
                          {
                            oomph_info
                              << "Number of master nodes for hanging value "
                              << i << " of node, " << nmaster
                              << ", does not match number of master "
                              << "nodes on other proc, " << nmaster_other
                              << std::endl;
                            shout = true;
                            shout_and_terminate = true;
                          }
                        }
                        // Value is not hanging on local proc: It had better
                        // not have any masters (i.e. be hanging) on the other
                        // proc
                        else
                        {
                          if (nmaster_other != 0)
                          {
                            oomph_info
                              << "Value " << i
                              << " of node is not hanging whereas "
                              << " node on other proc has " << nmaster_other
                              << " masters and therefore is hanging. \n";
                            shout = true;
                            shout_and_terminate = true;
                          }
                        }
                      }

                      if (shout)
                      {
                        // Report error. NOTE: ERROR IS THROWN BELOW ONCE
                        // ALL THIS HAS BEEN PROCESSED.

                        oomph_info
                          << "Error(s) displayed above are for "
                          << "domain with non-halo (i.e. haloed) elem: " << dd
                          << "\n";
                        oomph_info
                          << "Domain with    halo                elem: " << d
                          << "\n";
                        switch (nod_dim)
                        {
                          case 1:
                            oomph_info
                              << "Current processor is " << my_rank << "\n"
                              << "Nodal positions: " << x_halo[0] << "\n"
                              << "and haloed:      " << x_haloed[0]
                              << "\n"
                              //<< "Node pointer: " << finite_el_pt->node_pt(j)
                              << "\n";
                            break;
                          case 2:
                            oomph_info
                              << "Current processor is " << my_rank << "\n"
                              << "Nodal positions: " << x_halo[0] << " "
                              << x_halo[1] << "\n"
                              << "and haloed:      " << x_haloed[0] << " "
                              << x_haloed[1]
                              << std::endl
                              //<< "Node pointer: " << finite_el_pt->node_pt(j)
                              << "\n";
                            break;
                          case 3:
                            oomph_info
                              << "Current processor is " << my_rank << "\n"
                              << "Nodal positions: " << x_halo[0] << " "
                              << x_halo[1] << " " << x_halo[2] << "\n"
                              << "and haloed:      " << x_haloed[0] << " "
                              << x_haloed[1] << " " << x_haloed[2]
                              << std::endl
                              //<< "Node pointer: " << finite_el_pt->node_pt(j)
                              << "\n";
                            break;
                          default:
                            throw OomphLibError(
                              "Nodal dimension not equal to 1, 2 or 3\n",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
                        }


                        // If documenting, write to output files
                        if (doc_info.is_doc_enabled())
                        {
                          for (unsigned i = 0; i < nod_dim; i++)
                          {
                            haloed_file << x_haloed[i] << " ";
                            halo_file << x_halo[i] << "  ";
                          }
                          haloed_file << error << " " << my_rank << " " << dd
                                      << " "
                                      << finite_el_pt->node_pt(j)->is_hanging()
                                      << std::endl;
                          halo_file << error << " " << my_rank << " " << dd
                                    << " " << other_geom_hanging << std::endl;
                        }
                      }
                    } // j<nnod_per_el
                  }
                } // e<nelem_haloed

                // If documenting, close output files
                if (doc_info.is_doc_enabled())
                {
                  haloed_file.close();
                  halo_file.close();
                }
              }
            }
          }
        }
      }
      // My haloed elements are not being checked: Send my halo elements
      // whose non-halo counterparts are located on processor d
      else
      {
        // Get vectors of halo elements by copy operation
        Vector<GeneralisedElement*> halo_elem_pt(halo_element_pt(d));

        // How many of my elements are halo elements whose non-halo
        // counterpart is located on processor d?
        unsigned nelem_halo = halo_elem_pt.size();

        if (nelem_halo != 0)
        {
          // Send it across to the processor whose haloed nodes are being
          // checked
          MPI_Send(&nelem_halo, 1, MPI_UNSIGNED, d, 0, Comm_pt->mpi_comm());

          // Only bother if the mesh consists of finite elements
          if (dynamic_cast<FiniteElement*>(this->element_pt(0)))
          {
            // Now string together the nodal positions of all halo nodes
            // Use this only to work out roughly how much space to
            // reserve for the vector. Then we can push data cheaply
            // while not assuming all elements have the same number
            // of nodes.
            unsigned nnod_first_el = finite_element_pt(0)->nnode();
            unsigned nod_dim = finite_element_pt(0)->node_pt(0)->ndim();
            Vector<double> nodal_positions;
            nodal_positions.reserve(nod_dim * nnod_first_el * nelem_halo);

            // Storage for hang information
            Vector<int> nodal_hangings;

            // unsigned count=0;
            for (unsigned e = 0; e < nelem_halo; e++)
            {
              FiniteElement* finite_el_pt =
                dynamic_cast<FiniteElement*>(halo_elem_pt[e]);
              if (finite_el_pt != 0)
              {
                unsigned nnod_this_el = finite_el_pt->nnode();
                for (unsigned j = 0; j < nnod_this_el; j++)
                {
                  Node* nod_pt = finite_el_pt->node_pt(j);
                  // unsigned nod_dim = nod_pt->ndim();

                  // Throw error if node doesn't exist
                  if (nod_pt == 0)
                  {
                    // Print our nodes in element
                    oomph_info << "element: " << finite_el_pt << std::endl;
                    for (unsigned i = 0; i < finite_el_pt->nnode(); i++)
                    {
                      oomph_info << finite_el_pt->node_pt(i) << std::endl;
                    }

                    // Throw error
                    throw OomphLibError("Node doesn't exist!",
                                        OOMPH_CURRENT_FUNCTION,
                                        OOMPH_EXCEPTION_LOCATION);
                  }

                  // Testing POSITIONS, not x location (cf hanging nodes,
                  // nodes.h)
                  for (unsigned i = 0; i < nod_dim; i++)
                  {
                    // nodal_positions[count]=nod_pt->position(i);
                    // count++;
                    nodal_positions.push_back(nod_pt->position(i));
                  }

                  unsigned nval = nod_pt->nvalue();
                  nodal_hangings.push_back(nval);
                  for (int i = -1; i < int(nval); i++)
                  {
                    if (nod_pt->is_hanging(i))
                    {
                      unsigned nmaster = nod_pt->hanging_pt(i)->nmaster();
                      nodal_hangings.push_back(nmaster);
                    }
                    else
                    {
                      nodal_hangings.push_back(0);
                    }
                  }
                }
              }
            }

            // Total number of nodal positions to be checked
            unsigned n_nodal_positions = nodal_positions.size();

            // Total number of nodal hang information to be checked
            unsigned n_nodal_hangings = nodal_hangings.size();

            // Send it across to the processor whose haloed elements are being
            // checked
            // MPI_Send(&nodal_positions[0],nod_dim*nnod_per_el*nelem_halo,
            //         MPI_DOUBLE,d,0,comm_pt->mpi_comm());
            MPI_Send(
              &n_nodal_positions, 1, MPI_UNSIGNED, d, 0, Comm_pt->mpi_comm());
            if (n_nodal_positions > 0)
            {
              MPI_Send(&nodal_positions[0],
                       n_nodal_positions,
                       MPI_DOUBLE,
                       d,
                       0,
                       Comm_pt->mpi_comm());
            }
            MPI_Send(
              &n_nodal_hangings, 1, MPI_UNSIGNED, d, 1, Comm_pt->mpi_comm());
            if (n_nodal_hangings > 0)
            {
              MPI_Send(&nodal_hangings[0],
                       n_nodal_hangings,
                       MPI_INT,
                       d,
                       1,
                       Comm_pt->mpi_comm());
            }
          }
        }
      }
    }

    oomph_info << "Max. error for halo/haloed elements " << max_error
               << std::endl;
    if (max_error > max_permitted_error_for_halo_check)
    {
      shout_and_terminate = true;
      oomph_info << "This is bigger than the permitted threshold "
                 << max_permitted_error_for_halo_check << std::endl;
      oomph_info
        << "If you believe this to be acceptable for your problem\n"
        << "increase Problem::Max_permitted_error_for_halo_check and re-run \n";
    }

    if (shout_and_terminate)
    {
      throw OomphLibError("Error in halo checking",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Now check the halo/haloed nodes lookup schemes

    // Doc halo/haloed nodes lookup schemes
    //-------------------------------------
    if (doc_info.is_doc_enabled())
    {
      // Loop over domains for halo nodes
      for (int dd = 0; dd < n_proc; dd++)
      {
        filename.str("");
        filename << doc_info.directory() << "/halo_node_check"
                 << doc_info.label() << "_on_proc" << my_rank << "_with_proc"
                 << dd << "_" << doc_info.number() << ".dat";
        halo_file.open(filename.str().c_str());
        halo_file << "ZONE " << std::endl;

        unsigned nnod = nhalo_node(dd);
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = halo_node_pt(dd, j);
          unsigned ndim = nod_pt->ndim();
          for (unsigned i = 0; i < ndim; i++)
          {
            halo_file << nod_pt->position(i) << " ";
          }
          halo_file << nod_pt->is_hanging() << std::endl;
        }
        // Dummy output for processor that doesn't share halo nodes
        // (needed for tecplot)
        // (This will only work if there are elements on this processor...)
        if ((nnod == 0) && (nelement() != 0))
        {
          FiniteElement* f_el_pt = dynamic_cast<FiniteElement*>(element_pt(0));
          // If it's a generalised element mesh dummy output
          if (f_el_pt == 0)
          {
            halo_file << "0.0" << std::endl;
          }
          else
          {
            unsigned ndim = f_el_pt->node_pt(0)->ndim();
            if (ndim == 2)
            {
              halo_file << " 1.0 1.1 " << std::endl;
            }
            else
            {
              halo_file << " 1.0 1.1 1.1" << std::endl;
            }
          }
        }
        halo_file.close();
      }


      // Loop over domains for haloed nodes
      for (int d = 0; d < n_proc; d++)
      {
        filename.str("");
        filename << doc_info.directory() << "/haloed_node_check"
                 << doc_info.label() << "_on_proc" << my_rank << "_with_proc"
                 << d << "_" << doc_info.number() << ".dat";
        haloed_file.open(filename.str().c_str());
        haloed_file << "ZONE " << std::endl;

        unsigned nnod = nhaloed_node(d);
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = haloed_node_pt(d, j);
          unsigned ndim = nod_pt->ndim();
          for (unsigned i = 0; i < ndim; i++)
          {
            haloed_file << nod_pt->position(i) << " ";
          }
          haloed_file << nod_pt->is_hanging() << std::endl;
        }
        // Dummy output for processor that doesn't share halo nodes
        // (needed for tecplot)
        if ((nnod == 0) && (nelement() != 0))
        {
          FiniteElement* f_el_pt = dynamic_cast<FiniteElement*>(element_pt(0));
          // If it's a generalised element mesh dummy output
          if (f_el_pt == 0)
          {
            haloed_file << "0.0" << std::endl;
          }
          else
          {
            unsigned ndim = f_el_pt->node_pt(0)->ndim();
            if (ndim == 2)
            {
              haloed_file << " 1.0 1.1 " << std::endl;
            }
            else
            {
              haloed_file << " 1.0 1.1 1.1" << std::endl;
            }
          }
        }
        haloed_file.close();
      }
    }

    // Check halo/haloed nodes lookup schemes
    //---------------------------------------
    max_error = 0.0;

    // Loop over domains for haloed nodes
    for (int d = 0; d < n_proc; d++)
    {
      // Are my haloed nodes being checked?
      if (d == my_rank)
      {
        // Loop over domains for halo nodes
        for (int dd = 0; dd < n_proc; dd++)
        {
          // Don't talk to yourself
          if (dd != d)
          {
            // How many of my nodes are haloed nodes whose halo
            // counterpart is located on processor dd?
            int nnod_haloed = nhaloed_node(dd);

            if (nnod_haloed != 0)
            {
              // Receive from processor dd how many of his nodes are halo
              // nodes whose non-halo counterparts are located here
              int nnod_halo = 0;
              MPI_Recv(
                &nnod_halo, 1, MPI_INT, dd, 0, Comm_pt->mpi_comm(), &status);

              if (nnod_haloed != nnod_halo)
              {
                std::ostringstream error_message;

                error_message << "Clash in numbers of halo and haloed nodes! "
                              << std::endl;
                error_message
                  << "# of haloed nodes whose halo counterpart lives on proc "
                  << dd << ": " << nnod_haloed << std::endl;
                error_message
                  << "# of halo nodes whose non-halo counterpart lives on proc "
                  << d << ": " << nnod_halo << std::endl;
                error_message
                  << "(Re-)run Mesh::check_halo_schemes() with DocInfo object"
                  << std::endl;
                error_message << "to identify the problem" << std::endl;
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }


              unsigned nod_dim = finite_element_pt(0)->node_pt(0)->ndim();

              // Get strung-together nodal positions from other processor
              Vector<double> other_nodal_positions(nod_dim * nnod_halo);
              MPI_Recv(&other_nodal_positions[0],
                       nod_dim * nnod_halo,
                       MPI_DOUBLE,
                       dd,
                       0,
                       Comm_pt->mpi_comm(),
                       &status);

              // Check
              unsigned count = 0;
              for (int j = 0; j < nnod_halo; j++)
              {
                Vector<double> x_haloed(nod_dim);
                for (unsigned i = 0; i < nod_dim; i++)
                {
                  x_haloed[i] = haloed_node_pt(dd, j)->position(i);
                }
                Vector<double> x_halo(nod_dim);
                for (unsigned i = 0; i < nod_dim; i++)
                {
                  x_halo[i] = other_nodal_positions[count];
                  count++;
                }
                double error = 0.0;
                for (unsigned i = 0; i < nod_dim; i++)
                {
                  error +=
                    (x_haloed[i] - x_halo[i]) * (x_haloed[i] - x_halo[i]);
                }
                error = sqrt(error);
                if (error > max_error)
                {
                  max_error = error;
                }
              }
            }
          }
        }
      }
      // My haloed nodes are not being checked: Send my halo nodes
      // whose non-halo counterparts are located on processor d
      else
      {
        int nnod_halo = nhalo_node(d);

        if (nnod_halo != 0)
        {
          // Send it across to the processor whose haloed nodes are being
          // checked
          MPI_Send(&nnod_halo, 1, MPI_INT, d, 0, Comm_pt->mpi_comm());

          unsigned nod_dim = finite_element_pt(0)->node_pt(0)->ndim();

          // Now string together the nodal positions of all halo nodes
          Vector<double> nodal_positions(nod_dim * nnod_halo);
          unsigned count = 0;
          for (int j = 0; j < nnod_halo; j++)
          {
            for (unsigned i = 0; i < nod_dim; i++)
            {
              nodal_positions[count] = halo_node_pt(d, j)->position(i);
              count++;
            }
          }
          // Send it across to the processor whose haloed nodes are being
          // checked
          MPI_Send(&nodal_positions[0],
                   nod_dim * nnod_halo,
                   MPI_DOUBLE,
                   d,
                   0,
                   Comm_pt->mpi_comm());
        }
      }
    }

    oomph_info << "Max. error for halo/haloed nodes " << max_error << std::endl;

    if (max_error > max_permitted_error_for_halo_check)
    {
      std::ostringstream error_message;
      error_message << "This is bigger than the permitted threshold "
                    << max_permitted_error_for_halo_check << std::endl;
      error_message
        << "If you believe this to be acceptable for your problem\n"
        << "increase Problem::Max_permitted_error_for_halo_check and re-run \n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


    // Now check the external halo/haloed element lookup scheme

    // Doc external halo/haoloed element lookup schemes
    //--------------------------------------------------
    if (doc_info.is_doc_enabled())
    {
      // Loop over domains for external halo elements
      for (int dd = 0; dd < n_proc; dd++)
      {
        filename.str("");
        filename << doc_info.directory() << "/ext_halo_element_check"
                 << doc_info.label() << "_on_proc" << my_rank << "_with_proc"
                 << dd << "_" << doc_info.number() << ".dat";
        ext_halo_file.open(filename.str().c_str());
        output_external_halo_elements(dd, ext_halo_file);
        ext_halo_file.close();


        filename.str("");
        filename << doc_info.directory() << "/ext_halo_node_check"
                 << doc_info.label() << "_on_proc" << my_rank << "_with_proc"
                 << dd << "_" << doc_info.number() << ".dat";
        ext_halo_file.open(filename.str().c_str());

        // Get vectors of external halo/haloed elements by copy operation
        Vector<GeneralisedElement*> ext_halo_elem_pt(
          External_halo_element_pt[dd]);

        unsigned nelem = ext_halo_elem_pt.size();

        for (unsigned e = 0; e < nelem; e++)
        {
          ext_halo_file << "ZONE " << std::endl;
          // Can I cast to a finite element
          FiniteElement* finite_el_pt =
            dynamic_cast<FiniteElement*>(ext_halo_elem_pt[e]);
          if (finite_el_pt != 0)
          {
            unsigned nnod = finite_el_pt->nnode();
            for (unsigned j = 0; j < nnod; j++)
            {
              Node* nod_pt = finite_el_pt->node_pt(j);
              unsigned ndim = nod_pt->ndim();
              for (unsigned i = 0; i < ndim; i++)
              {
                ext_halo_file << nod_pt->position(i) << " ";
              }
              ext_halo_file << std::endl;
            }
          }
        }
        ext_halo_file.close();
      }

      // Loop over domains for external halo elements
      for (int d = 0; d < n_proc; d++)
      {
        filename.str("");
        filename << doc_info.directory() << "/ext_haloed_element_check"
                 << doc_info.label() << "_on_proc" << my_rank << "_with_proc"
                 << d << "_" << doc_info.number() << ".dat";
        ext_haloed_file.open(filename.str().c_str());
        output_external_haloed_elements(d, ext_haloed_file);
        ext_haloed_file.close();


        filename.str("");
        filename << doc_info.directory() << "/ext_haloed_node_check"
                 << doc_info.label() << "_on_proc" << my_rank << "_with_proc"
                 << d << "_" << doc_info.number() << ".dat";
        ext_haloed_file.open(filename.str().c_str());

        // Get vectors of external halo/haloed elements by copy operation
        Vector<GeneralisedElement*> ext_haloed_elem_pt(
          External_haloed_element_pt[d]);

        unsigned nelem2 = ext_haloed_elem_pt.size();
        for (unsigned e = 0; e < nelem2; e++)
        {
          ext_haloed_file << "ZONE " << std::endl;
          // Is it a finite element
          FiniteElement* finite_el_pt =
            dynamic_cast<FiniteElement*>(ext_haloed_elem_pt[e]);
          if (finite_el_pt != 0)
          {
            unsigned nnod2 = finite_el_pt->nnode();
            for (unsigned j = 0; j < nnod2; j++)
            {
              Node* nod_pt = finite_el_pt->node_pt(j);
              unsigned ndim = nod_pt->ndim();
              for (unsigned i = 0; i < ndim; i++)
              {
                ext_haloed_file << nod_pt->position(i) << " ";
              }
              ext_haloed_file << std::endl;
            }
          }
        }
        ext_haloed_file.close();
      }
    } // end of if doc flag

    // Check external halo/haloed element lookup schemes
    //--------------------------------------------------
    max_error = 0.0;
    shout = false;
    shout_and_terminate = false;

    // Loop over domains for external haloed elements
    for (int d = 0; d < n_proc; d++)
    {
      // Are my external haloed elements being checked?
      if (d == my_rank)
      {
        // Loop over domains for external halo elements
        for (int dd = 0; dd < n_proc; dd++)
        {
          // Don't talk to yourself
          if (dd != d)
          {
            // Get vectors of external haloed elements by copy operation
            Vector<GeneralisedElement*> ext_haloed_elem_pt(
              External_haloed_element_pt[dd]);

            // How many of my elements are external haloed elements whose halo
            // counterpart is located on processor dd?
            int nelem_haloed = ext_haloed_elem_pt.size();

            if (nelem_haloed != 0)
            {
              // Receive from processor dd how many of his elements are halo
              // nodes whose non-halo counterparts are located here
              int nelem_halo = 0;
              MPI_Recv(
                &nelem_halo, 1, MPI_INT, dd, 0, Comm_pt->mpi_comm(), &status);
              if (nelem_halo != nelem_haloed)
              {
                std::ostringstream error_message;
                error_message
                  << "Clash in numbers of external halo and haloed elements! "
                  << std::endl;
                error_message << "# of external haloed elements whose halo "
                                 "counterpart lives on proc "
                              << dd << ": " << nelem_haloed << std::endl;
                error_message << "# of external halo elements whose non-halo "
                                 "counterpart lives on proc "
                              << d << ": " << nelem_halo << std::endl;
                error_message << "(Re-)run Problem::check_halo_schemes() with "
                                 "DocInfo object"
                              << std::endl;
                error_message << "to identify the problem" << std::endl;
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }


              // We can only check nodal stuff for meshes of finite elements
              FiniteElement* fe_pt =
                dynamic_cast<FiniteElement*>(ext_haloed_elem_pt[0]);
              if (fe_pt != 0)
              {
                // Get strung-together elemental nodal positions
                // from other processor
                // unsigned nnod_per_el=fe_pt->nnode();
                unsigned nod_dim = fe_pt->node_pt(0)->ndim();
                // Vector<double> other_nodal_positions
                // (nod_dim*nnod_per_el*nelem_halo);
                // MPI_Recv(&other_nodal_positions[0],nod_dim*nnod_per_el*nelem_halo,
                //         MPI_DOUBLE,dd,0,comm_pt->mpi_comm(),&status);
                unsigned n_nodal_positions = 0;
                MPI_Recv(&n_nodal_positions,
                         1,
                         MPI_UNSIGNED,
                         dd,
                         0,
                         Comm_pt->mpi_comm(),
                         &status);
                Vector<double> other_nodal_positions(n_nodal_positions);
                if (n_nodal_positions > 0)
                {
                  MPI_Recv(&other_nodal_positions[0],
                           n_nodal_positions,
                           MPI_DOUBLE,
                           dd,
                           0,
                           Comm_pt->mpi_comm(),
                           &status);
                }


                // Receive hanging info to be checked
                Vector<int> other_nodal_hangings;
                unsigned n_other_nodal_hangings = 0;
                MPI_Recv(&n_other_nodal_hangings,
                         1,
                         MPI_UNSIGNED,
                         dd,
                         1,
                         Comm_pt->mpi_comm(),
                         &status);
                if (n_other_nodal_hangings > 0)
                {
                  other_nodal_hangings.resize(n_other_nodal_hangings);
                  MPI_Recv(&other_nodal_hangings[0],
                           n_other_nodal_hangings,
                           MPI_INT,
                           dd,
                           1,
                           Comm_pt->mpi_comm(),
                           &status);
                }

                // If documenting, open the output files
                if (doc_info.is_doc_enabled())
                {
                  filename.str("");
                  filename << doc_info.directory() << "/error_ext_haloed_check"
                           << doc_info.label() << "_on_proc" << my_rank
                           << "_with_proc" << dd << "_" << doc_info.number()
                           << ".dat";
                  ext_haloed_file.open(filename.str().c_str());
                  filename.str("");
                  filename << doc_info.directory() << "/error_ext_halo_check"
                           << doc_info.label() << "_on_proc" << my_rank
                           << "_with_proc" << dd << "_" << doc_info.number()
                           << ".dat";
                  ext_halo_file.open(filename.str().c_str());
                }

                unsigned count = 0;
                unsigned count_hanging = 0;
                for (int e = 0; e < nelem_haloed; e++)
                {
                  FiniteElement* finite_el_pt =
                    dynamic_cast<FiniteElement*>(ext_haloed_elem_pt[e]);

                  if (finite_el_pt != 0)
                  {
                    unsigned nnod_this_el = finite_el_pt->nnode();
                    for (unsigned j = 0; j < nnod_this_el; j++)
                    {
                      Node* nod_pt = finite_el_pt->node_pt(j);
                      // unsigned nod_dim = mod_pt->ndim();

                      // Testing POSITIONS, not x location
                      // (cf hanging nodes, nodes.h)
                      Vector<double> x_haloed(nod_dim);
                      for (unsigned i = 0; i < nod_dim; i++)
                      {
                        x_haloed[i] = nod_pt->position(i);
                      }

                      Vector<double> x_halo(nod_dim);
                      for (unsigned i = 0; i < nod_dim; i++)
                      {
                        x_halo[i] = other_nodal_positions[count];
                        ++count;
                      }

                      double error = 0.0;
                      shout = false;
                      for (unsigned i = 0; i < nod_dim; i++)
                      {
                        error +=
                          (x_haloed[i] - x_halo[i]) * (x_haloed[i] - x_halo[i]);
                      }
                      error = sqrt(error);

                      if (error > max_error)
                      {
                        max_error = error;
                      }
                      double tol = 1.0e-12;
                      if (error > tol)
                      {
                        oomph_info << "Discrepancy between nodal coordinates "
                                      "of external halo(ed)"
                                   << "element larger than tolerance (" << tol
                                   << ")\n  Error: " << error << "\n";
                        shout = true;
                      }

                      unsigned nval = nod_pt->nvalue();
                      int nval_other = other_nodal_hangings[count_hanging];
                      count_hanging++;
                      if (int(nval) != nval_other)
                      {
                        oomph_info
                          << "Number of values of node, " << nval
                          << ", does not match number of values on other proc, "
                          << nval_other << std::endl;
                        shout = true;
                        shout_and_terminate = true;
                      }

                      // Is other node geometrically hanging?
                      int other_geom_hanging = 0;

                      // Check hangingness/number of master nodes
                      for (int i = -1; i < int(nval); i++)
                      {
                        int nmaster_other = other_nodal_hangings[count_hanging];
                        count_hanging++;

                        // Record geom hang status of other node
                        if (i == -1) other_geom_hanging = nmaster_other;

                        // Value is hanging on local proc: Does it have the same
                        // number of masters as its counterpart on other proc?
                        if (nod_pt->is_hanging(i))
                        {
                          unsigned nmaster = nod_pt->hanging_pt(i)->nmaster();
                          if (int(nmaster) != nmaster_other)
                          {
                            oomph_info
                              << "Number of master nodes for hanging value "
                              << i << " of node, " << nmaster
                              << ", does not match number of master "
                              << "nodes on other proc, " << nmaster_other
                              << std::endl;
                            shout = true;
                            shout_and_terminate = true;
                          }
                        }
                        // Value is not hanging on local proc: It had better
                        // not have any masters (i.e. be hanging) on the other
                        // proc
                        else
                        {
                          if (nmaster_other != 0)
                          {
                            oomph_info
                              << "Value " << i
                              << " of node is not hanging whereas "
                              << " node on other proc has " << nmaster_other
                              << " masters and therefore is hanging. \n";
                            shout = true;
                            shout_and_terminate = true;
                          }
                        }
                      }

                      if (shout)
                      {
                        // Report error. NOTE: ERROR IS THROWN BELOW ONCE
                        // ALL THIS HAS BEEN PROCESSED.

                        oomph_info << "Error(s) displayed above are for "
                                   << "domain with external non-halo (i.e. "
                                      "haloed) elem: "
                                   << dd << "\n";
                        oomph_info
                          << "Domain with    halo                elem: " << d
                          << "\n";
                        switch (nod_dim)
                        {
                          case 1:
                            oomph_info
                              << "Current processor is " << my_rank << "\n"
                              << "Nodal positions: " << x_halo[0] << "\n"
                              << "and haloed:      " << x_haloed[0]
                              << "\n"
                              //<< "Node pointer: " << finite_el_pt->node_pt(j)
                              << "\n";
                            break;
                          case 2:
                            oomph_info
                              << "Current processor is " << my_rank << "\n"
                              << "Nodal positions: " << x_halo[0] << " "
                              << x_halo[1] << "\n"
                              << "and haloed:      " << x_haloed[0] << " "
                              << x_haloed[1]
                              << std::endl
                              //<< "Node pointer: " << finite_el_pt->node_pt(j)
                              << "\n";
                            break;
                          case 3:
                            oomph_info
                              << "Current processor is " << my_rank << "\n"
                              << "Nodal positions: " << x_halo[0] << " "
                              << x_halo[1] << " " << x_halo[2] << "\n"
                              << "and haloed:      " << x_haloed[0] << " "
                              << x_haloed[1] << " " << x_haloed[2]
                              << std::endl
                              //<< "Node pointer: " << finite_el_pt->node_pt(j)
                              << "\n";
                            break;
                          default:
                            throw OomphLibError(
                              "Nodal dimension not equal to 1, 2 or 3\n",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
                        }


                        // If documenting, write to output files
                        if (doc_info.is_doc_enabled())
                        {
                          for (unsigned i = 0; i < nod_dim; i++)
                          {
                            ext_haloed_file << x_haloed[i] << " ";
                            ext_halo_file << x_halo[i] << "  ";
                          }
                          ext_haloed_file
                            << error << " " << my_rank << " " << dd << " "
                            << finite_el_pt->node_pt(j)->is_hanging()
                            << std::endl;
                          ext_halo_file << error << " " << my_rank << " " << dd
                                        << " " << other_geom_hanging
                                        << std::endl;
                        }
                      }
                    } // j<nnod_per_el
                  }
                } // e<nelem_haloed

                // If documenting, close output files
                if (doc_info.is_doc_enabled())
                {
                  ext_haloed_file.close();
                  ext_halo_file.close();
                }
              }
            }
          }
        }
      }
      // My external haloed elements are not being checked: Send my halo
      // elements whose non-halo counterparts are located on processor d
      else
      {
        // Get vectors of external halo elements by copy operation
        Vector<GeneralisedElement*> ext_halo_elem_pt(
          External_halo_element_pt[d]);

        // How many of my elements are external halo elements whose non-halo
        // counterpart is located on processor d?
        unsigned nelem_halo = ext_halo_elem_pt.size();

        if (nelem_halo != 0)
        {
          // Send it across to the processor whose external haloed nodes
          // are being checked
          MPI_Send(&nelem_halo, 1, MPI_UNSIGNED, d, 0, Comm_pt->mpi_comm());

          // Only bother if the mesh consists of finite elements
          FiniteElement* fe_pt =
            dynamic_cast<FiniteElement*>(ext_halo_elem_pt[0]);
          if (fe_pt != 0)
          {
            // Now string together the nodal positions of all halo nodes
            unsigned nnod_first_el = fe_pt->nnode();
            unsigned nod_dim = fe_pt->node_pt(0)->ndim();
            Vector<double> nodal_positions;
            nodal_positions.reserve(nod_dim * nnod_first_el * nelem_halo);

            // Storage for hang information
            Vector<int> nodal_hangings;

            unsigned count = 0;
            for (unsigned e = 0; e < nelem_halo; e++)
            {
              FiniteElement* finite_el_pt =
                dynamic_cast<FiniteElement*>(ext_halo_elem_pt[e]);
              if (finite_el_pt != 0)
              {
                unsigned nnod_this_el = finite_el_pt->nnode();
                for (unsigned j = 0; j < nnod_this_el; j++)
                {
                  Node* nod_pt = finite_el_pt->node_pt(j);

                  // Testing POSITIONS, not x location (cf hanging nodes,
                  // nodes.h)
                  for (unsigned i = 0; i < nod_dim; i++)
                  {
                    nodal_positions.push_back(nod_pt->position(i));
                    count++;
                  }

                  unsigned nval = nod_pt->nvalue();
                  nodal_hangings.push_back(nval);
                  for (int i = -1; i < int(nval); i++)
                  {
                    if (nod_pt->is_hanging(i))
                    {
                      unsigned nmaster = nod_pt->hanging_pt(i)->nmaster();
                      nodal_hangings.push_back(nmaster);
                    }
                    else
                    {
                      nodal_hangings.push_back(0);
                    }
                  }
                }
              }
            }

            // Total number of nodal positions to be checked
            unsigned n_nodal_positions = nodal_positions.size();

            // Total number of nodal hang information to be checked
            unsigned n_nodal_hangings = nodal_hangings.size();

            // Send it across to the processor whose external haloed
            // elements are being checked
            // MPI_Send(&nodal_positions[0],nod_dim*nnod_per_el*nelem_halo,
            //         MPI_DOUBLE,d,0,comm_pt->mpi_comm());
            MPI_Send(
              &n_nodal_positions, 1, MPI_UNSIGNED, d, 0, Comm_pt->mpi_comm());
            if (n_nodal_positions > 0)
            {
              MPI_Send(&nodal_positions[0],
                       n_nodal_positions,
                       MPI_DOUBLE,
                       d,
                       0,
                       Comm_pt->mpi_comm());
            }
            MPI_Send(
              &n_nodal_hangings, 1, MPI_UNSIGNED, d, 1, Comm_pt->mpi_comm());
            if (n_nodal_hangings > 0)
            {
              MPI_Send(&nodal_hangings[0],
                       n_nodal_hangings,
                       MPI_INT,
                       d,
                       1,
                       Comm_pt->mpi_comm());
            }
          }
        }
      }
    }

    oomph_info << "Max. error for external halo/haloed elements " << max_error
               << std::endl;
    if (max_error > max_permitted_error_for_halo_check)
    {
      shout_and_terminate = true;
      oomph_info << "This is bigger than the permitted threshold "
                 << max_permitted_error_for_halo_check << std::endl;
      oomph_info
        << "If you believe this to be acceptable for your problem\n"
        << "increase Problem::Max_permitted_error_for_halo_check and re-run \n";
    }

    if (shout_and_terminate)
    {
      throw OomphLibError("Error in halo checking",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  }


  //========================================================================
  /// Null out specified external halo node (used when deleting duplicates)
  //========================================================================
  void Mesh::null_external_halo_node(const unsigned& p, Node* nod_pt)
  {
    // Loop over all external halo nodes with specified processor
    Vector<Node*> ext_halo_node_pt = External_halo_node_pt[p];
    unsigned n = ext_halo_node_pt.size();
    for (unsigned j = 0; j < n; j++)
    {
      if (ext_halo_node_pt[j] == nod_pt)
      {
        External_halo_node_pt[p][j] = 0;
        break;
      }
    }
  }


  //========================================================================
  /// Consolidate external halo node storage by removing nulled out
  /// pointes in external halo and haloed schemes
  //========================================================================
  void Mesh::remove_null_pointers_from_external_halo_node_storage()
  {
    // Storage for number of processors and current processor
    int n_proc = Comm_pt->nproc();
    int my_rank = Comm_pt->my_rank();

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
        // Make backup of external halo node pointers with this domain
        Vector<Node*> backup_pt = External_halo_node_pt[domain];

        // Wipe
        External_halo_node_pt[domain].clear();

        // How many do we have currently?
        unsigned nnod = backup_pt.size();
        External_halo_node_pt[domain].reserve(nnod);

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
            External_halo_node_pt[domain].push_back(nod_pt);
          }
        }
      }

      // End of data
      send_data.push_back(-1);

      // Find the number of data added to the vector
      send_n[domain] = send_data.size() - send_displacement[domain];
    }


    // Storage for the number of data to be received from each processor
    Vector<int> receive_n(n_proc, 0);

    // Now send numbers of data to be sent between all processors
    MPI_Alltoall(
      &send_n[0], 1, MPI_INT, &receive_n[0], 1, MPI_INT, Comm_pt->mpi_comm());


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
                  Comm_pt->mpi_comm());

    // Now use the received data
    for (int send_rank = 0; send_rank < n_proc; send_rank++)
    {
      // Don't bother to do anything for the processor corresponding to the
      // current processor or if no data were received from this processor
      if ((send_rank != my_rank) && (receive_n[send_rank] != 0))
      {
        // Counter for the data within the large array
        unsigned count = receive_displacement[send_rank];

        // Unpack until we reach "end of data" indicator (-1)
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
            External_haloed_node_pt[send_rank][next_one] = 0;
          }
        }

        // Make backup of external haloed node pointers with this domain
        Vector<Node*> backup_pt = External_haloed_node_pt[send_rank];

        // Wipe
        External_haloed_node_pt[send_rank].clear();

        // How many do we have currently?
        unsigned nnod = backup_pt.size();
        External_haloed_node_pt[send_rank].reserve(nnod);

        // Loop over external haloed nodes with this domain
        for (unsigned j = 0; j < nnod; j++)
        {
          // Get pointer to node
          Node* nod_pt = backup_pt[j];

          // Has it been nulled out?
          if (nod_pt != 0)
          {
            // Still alive: Copy across
            External_haloed_node_pt[send_rank].push_back(nod_pt);
          }
        }
      }

    } // End of data is received
  }

#endif


  // =================================================================
  /// Get the number of dof types in the mesh from the first element of the
  /// mesh. If MPI is on then also do some consistency checks between
  /// processors. \b Careful: Involves MPI Broadcasts and must therefore be
  /// called on all processors!
  // =================================================================
  unsigned Mesh::ndof_types() const
  {
    // Remains -1 if we don't have any elements on this processor.
    int int_ndof_types = -1;
    unsigned nel = nelement();
    if (nel > 0)
    {
      int_ndof_types = element_pt(0)->ndof_types();
#ifdef PARANOID
      // Check that every element in this mesh has the same number of
      // types of DOF.
      for (unsigned i = 1; i < nel; i++)
      {
        if (int_ndof_types != int(element_pt(i)->ndof_types()))
        {
          std::ostringstream error_message;
          error_message
            << "Every element in the mesh must have the same number of "
            << "types of DOF for ndof_types() to work\n"
            << "Element 0 has " << int_ndof_types << " DOF types\n"
            << "Element " << i << " [out of a total of " << nel << " ] has "
            << element_pt(i)->ndof_types() << " DOF types"
            << "Element types are: Element 0:" << typeid(*element_pt(0)).name()
            << "\n"
            << "           Current Element  :" << typeid(*element_pt(i)).name()
            << "\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
#endif
    }

#ifdef OOMPH_HAS_MPI

    // If mesh is distributed
    if (is_mesh_distributed())
    {
      // if more than one processor then
      // + ensure number of DOFs is consistent on each processor (PARANOID)
      // + ensure processors with no elements in this mesh have the
      //   correct number of DOF types
      if (Comm_pt->nproc() > 1)
      {
        unsigned nproc = Comm_pt->nproc();
        unsigned my_rank = Comm_pt->my_rank();

        // Collect on root the number of dofs types determined independently
        // on all processors (-1 indicates that the processor didn't have
        // any elements and therefore doesn't know!)
        int* ndof_types_recv = 0;
        if (my_rank == 0)
        {
          ndof_types_recv = new int[nproc];
        }

        MPI_Gather(&int_ndof_types,
                   1,
                   MPI_INT,
                   ndof_types_recv,
                   1,
                   MPI_INT,
                   0,
                   Comm_pt->mpi_comm());

        // Root: Update own number of dof types, check consistency amongst
        // all processors (in paranoid mode) and send out the actual
        // number of dof types to those processors who couldn't figure this
        // out themselves
        if (my_rank == 0)
        {
          // Check number of types of all non-root processors
          for (unsigned p = 1; p < nproc; p++)
          {
            if (ndof_types_recv[p] != -1)
            {
              // Processor p was able to figure out how many
              // dof types there are, so I root can update
              // its own (if required)
              if (int_ndof_types == -1)
              {
                int_ndof_types = ndof_types_recv[p];
              }
#ifdef PARANOID
              // Check consistency
              else if (int_ndof_types != ndof_types_recv[p])
              {
                std::ostringstream error_message;
                error_message
                  << "The elements in this mesh must have the same number "
                  << "of types of DOF on each processor";
                for (unsigned p = 0; p < nproc; p++)
                {
                  if (ndof_types_recv[p] != -1)
                  {
                    error_message << "Processor " << p << " : "
                                  << ndof_types_recv[p] << "\n";
                  }
                  else
                  {
                    error_message << "Processor " << p << " : (no elements)\n";
                  }
                }
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif
            }
          }

          // Now send ndof types to non-root processors that don't have it
          for (unsigned p = 1; p < nproc; p++)
          {
            if (ndof_types_recv[p] == -1)
            {
              MPI_Send(&int_ndof_types, 1, MPI_INT, p, 0, Comm_pt->mpi_comm());
            }
          }
          // clean up
          delete[] ndof_types_recv;
        }
        // "else if": "else" for non-root; "if" for checking if current
        // (non-root) processor does not know ndof type and is therefore
        // about to receive it from root.
        else if (int_ndof_types == -1)
        {
          MPI_Recv(&int_ndof_types,
                   1,
                   MPI_INT,
                   0,
                   0,
                   Comm_pt->mpi_comm(),
                   MPI_STATUS_IGNORE);
        }
      }
    }
#endif

    // If int_ndof_types if still -1 then no elements were found for this mesh,
    // so it has no dofs.
    if (int_ndof_types == -1) int_ndof_types = 0;

    return unsigned(int_ndof_types);
  }

  // =================================================================
  /// Get the number of elemental dimension in the mesh from the first
  /// element of the mesh. If MPI is on then also do some consistency
  /// checks between processors. \b Careful: Involves MPI Broadcasts
  /// and must therefore be called on all processors!
  // =================================================================
  unsigned Mesh::elemental_dimension() const
  {
    // Remains -1 if we don't have any elements on this processor.
    int int_dim = -1;
    if (nelement() > 0)
    {
      int_dim = finite_element_pt(0)->dim();
#ifdef PARANOID
      // Check that every element in this mesh has the same number of
      // types of elemental dimension.
      for (unsigned i = 1; i < nelement(); i++)
      {
        if (int_dim != int(finite_element_pt(i)->dim()))
        {
          std::ostringstream error_message;
          error_message
            << "Every element in the mesh must have the same number of "
            << "elemental dimension for elemental_dimension() to work.\n"
            << "Element 0 has elemental dimension " << int_dim << "\n"
            << "Element " << i << " has elemental dimension "
            << finite_element_pt(i)->dim() << ".";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
#endif
    }

#ifdef OOMPH_HAS_MPI

    // If mesh is distributed
    if (Comm_pt != 0)
    {
      // if more than one processor then
      // + ensure dimension number is consistent on each processor (PARANOID)
      // + ensure processors with no elements in this mesh have the
      //   correct dimension number.
      if (Comm_pt->nproc() > 1)
      {
        unsigned nproc = Comm_pt->nproc();
        unsigned my_rank = Comm_pt->my_rank();

        // Collect on root the dimension number determined independently
        // on all processors (-1 indicates that the processor didn't have
        // any elements and therefore doesn't know!)
        int* dim_recv = 0;
        if (my_rank == 0)
        {
          dim_recv = new int[nproc];
        }

        MPI_Gather(
          &int_dim, 1, MPI_INT, dim_recv, 1, MPI_INT, 0, Comm_pt->mpi_comm());

        // Root: Update own dimension, check consistency amongst
        // all processors (in paranoid mode) and send out the actual
        // dimension number to those processors who couldn't figure this
        // out themselves
        if (my_rank == 0)
        {
          // Check number of types of all non-root processors
          for (unsigned p = 1; p < nproc; p++)
          {
            if (dim_recv[p] != -1)
            {
              // Processor p was able to figure out the elemental
              // dimension, so I root can update
              // its own (if required)
              if (int_dim == -1)
              {
                int_dim = dim_recv[p];
              }
#ifdef PARANOID
              // Check consistency
              else if (int_dim != dim_recv[p])
              {
                std::ostringstream error_message;
                error_message
                  << "The elements in this mesh must have the same elemental "
                  << "dimension number on each processor";
                for (unsigned p = 0; p < nproc; p++)
                {
                  if (dim_recv[p] != -1)
                  {
                    error_message << "Processor " << p << " : " << dim_recv[p]
                                  << "\n";
                  }
                  else
                  {
                    error_message << "Processor " << p << " : (no elements)\n";
                  }
                }
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif
            }
          }

          // Now send the elemental dimension to non-root processors that
          // don't have it
          for (unsigned p = 1; p < nproc; p++)
          {
            if (dim_recv[p] == -1)
            {
              MPI_Send(&int_dim, 1, MPI_INT, p, 0, Comm_pt->mpi_comm());
            }
          }
          // clean up
          delete[] dim_recv;
        }
        // "else if": "else" for non-root; "if" for checking if current
        // (non-root) processor does not know elemental dimension and is
        // therefore about to receive it from root.
        else if (int_dim == -1)
        {
          MPI_Recv(
            &int_dim, 1, MPI_INT, 0, 0, Comm_pt->mpi_comm(), MPI_STATUS_IGNORE);
        }
      }
    }
#endif

    // If int_dim if still -1 then no elements were found for this mesh, so it
    // has no elemental dimension.
    if (int_dim == -1) int_dim = 0;

    return unsigned(int_dim);
  }

  // =================================================================
  /// Get the number of nodal dimension in the mesh from the first
  /// element of the mesh. If MPI is on then also do some consistency
  /// checks between processors. \b Careful: Involves MPI Broadcasts
  /// and must therefore be called on all processors!
  // =================================================================
  unsigned Mesh::nodal_dimension() const
  {
    // Remains -1 if we don't have any elements on this processor.
    int int_dim = -1;
    if (nelement() > 0)
    {
      int_dim = finite_element_pt(0)->nodal_dimension();
#ifdef PARANOID
      // Check that every element in this mesh has the same number of
      // types of nodal dimension.
      for (unsigned i = 1; i < nelement(); i++)
      {
        if (int_dim != int(finite_element_pt(i)->nodal_dimension()))
        {
          std::ostringstream error_message;
          error_message
            << "Every element in the mesh must have the same number of "
            << "nodal dimension for nodal_dimension() to work.\n"
            << "Element 0 has nodal dimension " << int_dim << "\n"
            << "Element " << i << " has nodal dimension "
            << finite_element_pt(i)->nodal_dimension() << ".";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
#endif
    }

#ifdef OOMPH_HAS_MPI

    // If mesh is distributed
    if (Comm_pt != 0)
    {
      // if more than one processor then
      // + ensure dimension number is consistent on each processor (PARANOID)
      // + ensure processors with no elements in this mesh have the
      //   correct dimension number.
      if (Comm_pt->nproc() > 1)
      {
        unsigned nproc = Comm_pt->nproc();
        unsigned my_rank = Comm_pt->my_rank();

        // Collect on root the dimension number determined independently
        // on all processors (-1 indicates that the processor didn't have
        // any elements and therefore doesn't know!)
        int* dim_recv = 0;
        if (my_rank == 0)
        {
          dim_recv = new int[nproc];
        }

        MPI_Gather(
          &int_dim, 1, MPI_INT, dim_recv, 1, MPI_INT, 0, Comm_pt->mpi_comm());

        // Root: Update own dimension, check consistency amongst
        // all processors (in paranoid mode) and send out the actual
        // dimension number to those processors who couldn't figure this
        // out themselves
        if (my_rank == 0)
        {
          // Check number of types of all non-root processors
          for (unsigned p = 1; p < nproc; p++)
          {
            if (dim_recv[p] != -1)
            {
              // Processor p was able to figure out the nodal
              // dimension, so I root can update
              // its own (if required)
              if (int_dim == -1)
              {
                int_dim = dim_recv[p];
              }
#ifdef PARANOID
              // Check consistency
              else if (int_dim != dim_recv[p])
              {
                std::ostringstream error_message;
                error_message
                  << "The elements in this mesh must have the same nodal "
                  << "dimension number on each processor";
                for (unsigned p = 0; p < nproc; p++)
                {
                  if (dim_recv[p] != -1)
                  {
                    error_message << "Processor " << p << " : " << dim_recv[p]
                                  << "\n";
                  }
                  else
                  {
                    error_message << "Processor " << p << " : (no elements)\n";
                  }
                }
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif
            }
          }

          // Now send the nodal dimension to non-root processors that
          // don't have it
          for (unsigned p = 1; p < nproc; p++)
          {
            if (dim_recv[p] == -1)
            {
              MPI_Send(&int_dim, 1, MPI_INT, p, 0, Comm_pt->mpi_comm());
            }
          }
          // clean up
          delete[] dim_recv;
        }
        // "else if": "else" for non-root; "if" for checking if current
        // (non-root) processor does not know nodal dimension and is therefore
        // about to receive it from root.
        else if (int_dim == -1)
        {
          MPI_Recv(
            &int_dim, 1, MPI_INT, 0, 0, Comm_pt->mpi_comm(), MPI_STATUS_IGNORE);
        }
      }
    }
#endif

    // If int_dim if still -1 then no elements were found for this mesh, so it
    // has no nodal dimension.
    if (int_dim == -1) int_dim = 0;

    return unsigned(int_dim);
  }

  //========================================================================
  /// Wipe the storage for all externally-based elements and delete halos
  //========================================================================
  void Mesh::delete_all_external_storage()
  {
#ifdef OOMPH_HAS_MPI

    // Only do for distributed meshes
    if (is_mesh_distributed())
    {
      // Some of the external halo/haloed nodes are masters of nodes
      // in this mesh. We must set to be non-hanging any nodes whose
      // masters we are about to delete, to remove any dependencies.

      // Loop over all the mesh nodes and check their masters
      for (unsigned i = 0; i < nnode(); i++)
      {
        // Get pointer to the node
        Node* nod_pt = node_pt(i);

        // Check if the node exists
        if (nod_pt != 0)
        {
          // Check if the node is hanging
          if (nod_pt->is_hanging())
          {
            // Get pointer to the hang info
            HangInfo* hang_pt = nod_pt->hanging_pt();

            // Check if any master is in the external halo storage
            // External haloed nodes don't get deleted, so we don't need to
            //(and shouldn't) un-hang their dependents
            bool found_a_master_in_external_halo_storage = false;
            for (unsigned m = 0; m < hang_pt->nmaster(); m++)
            {
              // Iterator for vector of nodes
              Vector<Node*>::iterator it;

              // Loop over external halo storage with all processors
              bool found_this_master_in_external_halo_storage = false;
              for (int d = 0; d < Comm_pt->nproc(); d++)
              {
                // Find master in map of external halo nodes
                it = std::find(External_halo_node_pt[d].begin(),
                               External_halo_node_pt[d].end(),
                               hang_pt->master_node_pt(m));

                // Check if it was found
                if (it != External_halo_node_pt[d].end())
                {
                  // Mark as found
                  found_this_master_in_external_halo_storage = true;
                  // Don't need to search remaining processors
                  break;
                }
              }

              // Check if any have been found
              if (found_this_master_in_external_halo_storage)
              {
                // Mark as found
                found_a_master_in_external_halo_storage = true;
                // Don't need to search remaining masters
                break;
              }
            }

            // If it was found...
            if (found_a_master_in_external_halo_storage)
            {
              // Master is in external halo storage and is about to be deleted,
              // so we'd better make this node non-hanging. In case the node
              // does not become hanging again, we must get all the required
              // information from its masters to make it a 'proper' node again.

              // Reconstruct the nodal values/position from the node's
              // hanging node representation
              unsigned nt = nod_pt->ntstorage();
              unsigned n_value = nod_pt->nvalue();
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

              // If it's an algebraic node: Update its previous nodal positions
              // too
              AlgebraicNode* alg_node_pt = dynamic_cast<AlgebraicNode*>(nod_pt);
              if (alg_node_pt != 0)
              {
                bool update_all_time_levels = true;
                alg_node_pt->node_update(update_all_time_levels);
              }


              // If it's a Solid node, update Lagrangian coordinates
              // from its hanging node representation
              SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);
              if (solid_node_pt != 0)
              {
                unsigned n_lagrangian = solid_node_pt->nlagrangian();
                for (unsigned i = 0; i < n_lagrangian; i++)
                {
                  solid_node_pt->xi(i) = solid_node_pt->lagrangian_position(i);
                }
              }

              // No need to worry about geometrically hanging nodes
              // on boundaries (as in (p_)adapt_mesh())
              /// /Now store geometrically hanging nodes on boundaries that
              /// /may need updating after refinement.
              /// /There will only be a problem if we have 3 spatial dimensions
              // if((mesh_dim > 2) && (nod_pt->is_hanging()))
              // {
              //  //If the node is on a boundary then add a pointer to the node
              //  //to our lookup scheme
              //  if(nod_pt->is_on_boundary())
              //   {
              //    //Storage for the boundaries on which the Node is located
              //    std::set<unsigned>* boundaries_pt;
              //    nod_pt->get_boundaries_pt(boundaries_pt);
              //    if(boundaries_pt!=0)
              //     {
              //      //Loop over the boundaries and add a pointer to the node
              //      //to the appropriate storage scheme
              //      for(std::set<unsigned>::iterator
              //      it=boundaries_pt->begin();
              //          it!=boundaries_pt->end();++it)
              //       {
              //        hanging_nodes_on_boundary_pt[*it].insert(nod_pt);
              //       }
              //     }
              //   }
              // }

              // Finally set nonhanging
              nod_pt->set_nonhanging();
            }
          }
        }
        else
        {
          // Node doesn't exist!
        }
      }
    }

    // Careful: some of the external halo nodes are also in boundary
    //          node storage and should be removed from this first
    for (std::map<unsigned, Vector<Node*>>::iterator it =
           External_halo_node_pt.begin();
         it != External_halo_node_pt.end();
         it++)
    {
      // Processor ID
      int d = (*it).first;

      // How many external haloes with this process?
      unsigned n_ext_halo_nod = nexternal_halo_node(d);
      for (unsigned j = 0; j < n_ext_halo_nod; j++)
      {
        Node* ext_halo_nod_pt = external_halo_node_pt(d, j);
        unsigned n_bnd = nboundary();
        for (unsigned i_bnd = 0; i_bnd < n_bnd; i_bnd++)
        {
          // Call this for all boundaries; it will do nothing
          // if the node is not on the current boundary
          remove_boundary_node(i_bnd, ext_halo_nod_pt);
        }
      }
    }

    // A loop to delete external halo nodes
    for (std::map<unsigned, Vector<Node*>>::iterator it =
           External_halo_node_pt.begin();
         it != External_halo_node_pt.end();
         it++)
    {
      // Processor ID
      int d = (*it).first;

      unsigned n_ext_halo_nod = nexternal_halo_node(d);
      for (unsigned j = 0; j < n_ext_halo_nod; j++)
      {
        // Only delete if it's not a node stored in the current mesh
        bool is_a_mesh_node = false;
        unsigned n_node = nnode();
        for (unsigned jj = 0; jj < n_node; jj++)
        {
          if (Node_pt[jj] == External_halo_node_pt[d][j])
          {
            is_a_mesh_node = true;
          }
        }

        // There will also be duplications between multiple processors,
        // so make sure that we don't try to delete these twice
        if (!is_a_mesh_node)
        {
          // Loop over all other higher-numbered processors and check
          // for duplicated external halo nodes
          // (The highest numbered processor should delete all its ext halos)
          for (std::map<unsigned, Vector<Node*>>::iterator itt =
                 External_halo_node_pt.begin();
               itt != External_halo_node_pt.end();
               itt++)
          {
            // Processor ID
            int dd = (*itt).first;

            if (dd > d)
            {
              unsigned n_ext_halo = nexternal_halo_node(dd);
              for (unsigned jjj = 0; jjj < n_ext_halo; jjj++)
              {
                if (External_halo_node_pt[dd][jjj] ==
                    External_halo_node_pt[d][j])
                {
                  is_a_mesh_node = true;
                }
              }
            }
          }
        }

        // Only now if no duplicates exist can the node be safely deleted
        if (!is_a_mesh_node)
        {
          delete External_halo_node_pt[d][j];
        }
      }
    }

    // Another loop to delete external halo elements (which are distinct)
    for (std::map<unsigned, Vector<GeneralisedElement*>>::iterator it =
           External_halo_element_pt.begin();
         it != External_halo_element_pt.end();
         it++)
    {
      // Processor ID
      int d = (*it).first;

      unsigned n_ext_halo_el = nexternal_halo_element(d);
      for (unsigned e = 0; e < n_ext_halo_el; e++)
      {
        delete External_halo_element_pt[d][e];
      }
    }

    // Now we are okay to clear the external halo node storage
    External_halo_node_pt.clear();
    External_halo_element_pt.clear();

    // External haloed nodes and elements are actual members
    // of the external mesh and should not be deleted
    External_haloed_node_pt.clear();
    External_haloed_element_pt.clear();
#endif
  }


#ifdef OOMPH_HAS_MPI

  // NOTE: the add_external_haloed_node_pt and add_external_haloed_element_pt
  //       functions need to check whether the Node/FiniteElement argument
  //       has been added to the storage already; this is not the case
  //       for the add_external_halo_node_pt and add_external_halo_element_pt
  //       functions as these are newly-created elements that are created and
  //       added to the storage based on the knowledge of when their haloed
  //       counterparts were created and whether they were newly added

  //========================================================================
  /// Add external haloed element whose non-halo counterpart is held
  /// on processor p to the storage scheme for external haloed elements.
  /// If the element is already in the storage scheme then return its index
  //========================================================================
  unsigned Mesh::add_external_haloed_element_pt(const unsigned& p,
                                                GeneralisedElement*& el_pt)
  {
    // Loop over current storage
    unsigned n_extern_haloed = nexternal_haloed_element(p);

    // Is this already an external haloed element?
    bool already_external_haloed_element = false;
    unsigned external_haloed_el_index = 0;
    for (unsigned eh = 0; eh < n_extern_haloed; eh++)
    {
      if (el_pt == External_haloed_element_pt[p][eh])
      {
        // It's already there, so...
        already_external_haloed_element = true;
        // ...set the index of this element
        external_haloed_el_index = eh;
        break;
      }
    }

    // Has it been found?
    if (!already_external_haloed_element)
    {
      // Not found, so add it:
      External_haloed_element_pt[p].push_back(el_pt);
      // Return the index where it's just been added
      return n_extern_haloed;
    }
    else
    {
      // Return the index where it was found
      return external_haloed_el_index;
    }
  }

  //========================================================================
  /// Add external haloed node whose halo (external) counterpart
  /// is held on processor p to the storage scheme for external haloed nodes.
  /// If the node is already in the storage scheme then return its index
  //========================================================================
  unsigned Mesh::add_external_haloed_node_pt(const unsigned& p, Node*& nod_pt)
  {
    // Loop over current storage
    unsigned n_ext_haloed_nod = nexternal_haloed_node(p);

    // Is this already an external haloed node?
    bool is_an_external_haloed_node = false;
    unsigned external_haloed_node_index = 0;
    for (unsigned k = 0; k < n_ext_haloed_nod; k++)
    {
      if (nod_pt == External_haloed_node_pt[p][k])
      {
        is_an_external_haloed_node = true;
        external_haloed_node_index = k;
        break;
      }
    }

    // Has it been found?
    if (!is_an_external_haloed_node)
    {
      // Not found, so add it
      External_haloed_node_pt[p].push_back(nod_pt);
      // Return the index where it's just been added
      return n_ext_haloed_nod;
    }
    else
    {
      // Return the index where it was found
      return external_haloed_node_index;
    }
  }

#endif


  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  // Functions for solid meshes
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////


  //========================================================================
  /// Make the current configuration the undeformed one by
  /// setting the nodal Lagrangian coordinates to their current
  /// Eulerian ones
  //========================================================================
  void SolidMesh::set_lagrangian_nodal_coordinates()
  {
    // Find out how many nodes there are
    unsigned long n_node = nnode();

    // Loop over all the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // Cast node to solid node (can safely be done because
      // SolidMeshes consist of SolidNodes
      SolidNode* node_pt = static_cast<SolidNode*>(Node_pt[n]);

      // Number of Lagrangian coordinates
      unsigned n_lagrangian = node_pt->nlagrangian();

      // Number of generalised Lagrangian coordinates
      unsigned n_lagrangian_type = node_pt->nlagrangian_type();

      // The assumption here is that there must be fewer lagrangian coordinates
      // than eulerian (which must be true?)
      // Set (generalised) Lagrangian coords = (generalised) Eulerian coords
      for (unsigned k = 0; k < n_lagrangian_type; k++)
      {
        // Loop over lagrangian coordinates and set their values
        for (unsigned j = 0; j < n_lagrangian; j++)
        {
          node_pt->xi_gen(k, j) = node_pt->x_gen(k, j);
        }
      }
    }
  }


  //=======================================================================
  /// Static problem that can be used to assign initial conditions
  /// on a given mesh.
  //=======================================================================
  SolidICProblem SolidMesh::Solid_IC_problem;


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Namespace for paraview-style output helper functions
  //=================================================================
  namespace ParaviewHelper
  {
    /// Write the pvd file header
    void write_pvd_header(std::ofstream& pvd_file)
    {
      pvd_file << "<?xml version=\"1.0\"?>" << std::endl
               << "<VTKFile type=\"Collection\" version=\"0.1\">" << std::endl
               << "<Collection>" << std::endl;
    }

    /// Add name of output file and associated continuous time
    /// to pvd file.
    void write_pvd_information(std::ofstream& pvd_file,
                               const std::string& output_filename,
                               const double& time)
    {
      // Output the actual time values
      pvd_file << "<DataSet timestep=\"" << time << "\" ";

      // Apparently this has to go in
      pvd_file << "part=\"0\" ";

      // Add the name of the file, so that the pvd file knows what it is called
      pvd_file << "file=\"" << output_filename << "\"/>" << std::endl;
    }

    /// Write the pvd file footer
    void write_pvd_footer(std::ofstream& pvd_file)
    {
      pvd_file << "</Collection>" << std::endl << "</VTKFile>";
    }

  } // namespace ParaviewHelper

  /// /////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////


} // namespace oomph
