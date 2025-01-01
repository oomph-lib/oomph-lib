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
// Non-inline member functions for hp-refineable elements

// oomph-lib includes
#include "algebraic_elements.h"
#include "macro_element_node_update_element.h"
#include "hp_refineable_elements.h"
//#include "shape.h"

namespace oomph
{
  /// /////////////////////////////////////////////////////////////
  //       1D PRefineableQElements
  /// /////////////////////////////////////////////////////////////

  /// Get local coordinates of node j in the element; vector sets its own size
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<1, INITIAL_NNODE_1D>::local_coordinate_of_node(
    const unsigned& n, Vector<double>& s) const
  {
    s.resize(1);

    switch (this->nnode_1d())
    {
      case 2:
        OneDimensionalLegendreShape<2>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<2>::nodal_position(n);
        break;
      case 3:
        OneDimensionalLegendreShape<3>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<3>::nodal_position(n);
        break;
      case 4:
        OneDimensionalLegendreShape<4>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<4>::nodal_position(n);
        break;
      case 5:
        OneDimensionalLegendreShape<5>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<5>::nodal_position(n);
        break;
      case 6:
        OneDimensionalLegendreShape<6>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<6>::nodal_position(n);
        break;
      case 7:
        OneDimensionalLegendreShape<7>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<7>::nodal_position(n);
        break;
      default:
        oomph_info << "\n ERROR: Exceeded maximum polynomial order for";
        oomph_info << "\n        shape functions." << std::endl;
        break;
    }
  }

  /// Get the local fractino of node j in the element
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<1, INITIAL_NNODE_1D>::local_fraction_of_node(
    const unsigned& n, Vector<double>& s_fraction)
  {
    this->local_coordinate_of_node(n, s_fraction);
    s_fraction[0] = 0.5 * (s_fraction[0] + 1.0);
  }

  template<unsigned INITIAL_NNODE_1D>
  double PRefineableQElement<1, INITIAL_NNODE_1D>::local_one_d_fraction_of_node(
    const unsigned& n1d, const unsigned& i)
  {
    switch (this->nnode_1d())
    {
      case 2:
        OneDimensionalLegendreShape<2>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<2>::nodal_position(n1d) + 1.0);
      case 3:
        OneDimensionalLegendreShape<3>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<3>::nodal_position(n1d) + 1.0);
      case 4:
        OneDimensionalLegendreShape<4>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<4>::nodal_position(n1d) + 1.0);
      case 5:
        OneDimensionalLegendreShape<5>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<5>::nodal_position(n1d) + 1.0);
      case 6:
        OneDimensionalLegendreShape<6>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<6>::nodal_position(n1d) + 1.0);
      case 7:
        OneDimensionalLegendreShape<7>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<7>::nodal_position(n1d) + 1.0);
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       shape functions.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        return 0.0;
    }
  }


  //==================================================================
  /// Return the node at the specified local coordinate
  //==================================================================
  template<unsigned INITIAL_NNODE_1D>
  Node* PRefineableQElement<1, INITIAL_NNODE_1D>::get_node_at_local_coordinate(
    const Vector<double>& s) const
  {
    // Load the tolerance into a local variable
    double tol = FiniteElement::Node_location_tolerance;
    // There is one possible index.
    Vector<int> index(1);

    // Determine the index
    // -------------------

    // If we are at the lower limit, the index is zero
    if (std::fabs(s[0] + 1.0) < tol)
    {
      index[0] = 0;
    }
    // If we are at the upper limit, the index is the number of nodes minus 1
    else if (std::fabs(s[0] - 1.0) < tol)
    {
      index[0] = this->nnode_1d() - 1;
    }
    // Otherwise, we have to calculate the index in general
    else
    {
      // Compute Gauss-Lobatto-Legendre node positions
      Vector<double> z;
      Orthpoly::gll_nodes(this->nnode_1d(), z);
      // Loop over possible internal nodal positions
      for (unsigned n = 1; n < this->nnode_1d() - 1; n++)
      {
        if (std::fabs(z[n] - s[0]) < tol)
        {
          index[0] = n;
          break;
        }
      }
      // No matching nodes
      return 0;
    }
    // If we've got here we have a node, so let's return a pointer to it
    return this->node_pt(index[0]);
  }

  //===================================================================
  /// If a neighbouring element's son has already created a node at
  /// a position corresponding to the local fractional position within the
  /// present element, s_fraction, return
  /// a pointer to that node. If not, return NULL (0). If the node is
  /// periodic the flag is_periodic will be true
  //===================================================================
  template<unsigned INITIAL_NNODE_1D>
  Node* PRefineableQElement<1, INITIAL_NNODE_1D>::
    node_created_by_son_of_neighbour(const Vector<double>& s_fraction,
                                     bool& is_periodic)
  {
    // Not possible in 1D case, so return null pointer
    return 0;
  }

  //==================================================================
  /// Set the correct p-order of the element based on that of its
  /// father. Then construct an integration scheme of the correct order.
  /// If an adopted father is specified, information from this is
  /// used instead of using the father found from the tree.
  //==================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<1, INITIAL_NNODE_1D>::initial_setup(
    Tree* const& adopted_father_pt, const unsigned& initial_p_order)
  {
    // Storage for pointer to my father (in binarytree impersonation)
    BinaryTree* father_pt = 0;

    // Check if an adopted father has been specified
    if (adopted_father_pt != 0)
    {
      // Get pointer to my father (in binarytree impersonation)
      father_pt = dynamic_cast<BinaryTree*>(adopted_father_pt);
    }
    // Check if element is in a tree
    else if (Tree_pt != 0)
    {
      // Get pointer to my father (in binarytree impersonation)
      father_pt = dynamic_cast<BinaryTree*>(binary_tree_pt()->father_pt());
    }
    // else
    // {
    //  throw OomphLibError(
    //         "Element not in a tree, and no adopted father has been
    //         specified!", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    // }

    // Check if element has father
    if (father_pt != 0 || initial_p_order != 0)
    {
      if (father_pt != 0)
      {
        PRefineableQElement<1, INITIAL_NNODE_1D>* father_el_pt =
          dynamic_cast<PRefineableQElement<1, INITIAL_NNODE_1D>*>(
            father_pt->object_pt());
        if (father_el_pt != 0)
        {
          unsigned father_p_order = father_el_pt->p_order();
          // Set p-order to that of father
          P_order = father_p_order;
        }
      }
      else
      {
        P_order = initial_p_order;
      }

      // Now sort out the element...
      // (has p nodes)
      unsigned new_n_node = P_order;

      // Allocate new space for Nodes (at the element level)
      this->set_n_node(new_n_node);

      // Set integration scheme
      delete this->integral_pt();
      switch (P_order)
      {
        case 2:
          this->set_integration_scheme(new GaussLobattoLegendre<1, 2>);
          break;
        case 3:
          this->set_integration_scheme(new GaussLobattoLegendre<1, 3>);
          break;
        case 4:
          this->set_integration_scheme(new GaussLobattoLegendre<1, 4>);
          break;
        case 5:
          this->set_integration_scheme(new GaussLobattoLegendre<1, 5>);
          break;
        case 6:
          this->set_integration_scheme(new GaussLobattoLegendre<1, 6>);
          break;
        case 7:
          this->set_integration_scheme(new GaussLobattoLegendre<1, 7>);
          break;
        default:
          std::ostringstream error_message;
          error_message << "\nERROR: Exceeded maximum polynomial order for";
          error_message << "\n       integration scheme.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }
  }

  //==================================================================
  /// Check the father element for any required nodes which
  /// already exist
  //==================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<1, INITIAL_NNODE_1D>::pre_build(
    Mesh*& mesh_pt, Vector<Node*>& new_node_pt)
  {
    /*
      //Pointer to my father (in binarytree impersonation)
      BinaryTree* father_pt =
      dynamic_cast<BinaryTree*>(binary_tree_pt()->father_pt());

      // Check if element has father
      if (father_pt!=0)
      {
      PRefineableQElement<1>* father_el_pt =
      dynamic_cast<PRefineableQElement<1>*>
      (this->tree_pt()->father_pt()->object_pt());
      if (father_el_pt!=0)
      {
      // Pre-build actions
      //??
      }
      else
      {
      std::ostringstream error_message;
      error_message <<"\nERROR: Dynamic cast failed!\n";
      throw OomphLibError(error_message.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
      }
      }
    */
  }

  //==================================================================
  /// p-refine the element inc times. (If inc<0 then p-unrefinement
  /// is performed.)
  //==================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<1, INITIAL_NNODE_1D>::p_refine(
    const int& inc, Mesh* const& mesh_pt, GeneralisedElement* const& clone_pt)
  {
    // Cast clone to correct type
    PRefineableQElement<1, INITIAL_NNODE_1D>* clone_el_pt =
      dynamic_cast<PRefineableQElement<1, INITIAL_NNODE_1D>*>(clone_pt);

    // Check if we can use it
    if (clone_el_pt == 0)
    {
      throw OomphLibError(
        "Cloned copy must be a PRefineableQElement<1,INITIAL_NNODE_1D>!",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#ifdef PARANOID
    // Clone exists, so check that it is infact a copy of me
    else
    {
      // Flag to keep track of check
      bool clone_is_ok = true;

      // Does the clone have the correct p-order?
      clone_is_ok = clone_is_ok && (clone_el_pt->p_order() == this->p_order());

      if (!clone_is_ok)
      {
        std::ostringstream err_stream;
        err_stream << "Clone element has a different p-order from me,"
                   << std::endl
                   << "but it is supposed to be a copy!" << std::endl;

        throw OomphLibError(
          err_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Does the clone have the same number of nodes as me?
      clone_is_ok = clone_is_ok && (clone_el_pt->nnode() == this->nnode());

      if (!clone_is_ok)
      {
        std::ostringstream err_stream;
        err_stream << "Clone element has a different number of nodes from me,"
                   << std::endl
                   << "but it is supposed to be a copy!" << std::endl;

        throw OomphLibError(
          err_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Does the clone have the same nodes as me?
      for (unsigned n = 0; n < this->nnode(); n++)
      {
        clone_is_ok =
          clone_is_ok && (clone_el_pt->node_pt(n) == this->node_pt(n));
      }

      if (!clone_is_ok)
      {
        std::ostringstream err_stream;
        err_stream << "The nodes belonging to the clone element are different"
                   << std::endl
                   << "from mine, but it is supposed to be a copy!"
                   << std::endl;

        throw OomphLibError(
          err_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // If we get to here then the clone has all the information we require
    }
#endif

    // Currently we can't handle the case of generalised coordinates
    // since we haven't established how they should be interpolated.
    // Buffer this case:
    if (clone_el_pt->node_pt(0)->nposition_type() != 1)
    {
      throw OomphLibError("Can't handle generalised nodal positions (yet).",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Timestepper should be the same for all nodes -- use it
    // to create timesteppers for new nodes
    TimeStepper* time_stepper_pt = this->node_pt(0)->time_stepper_pt();

    // Get number of history values (incl. present)
    unsigned ntstorage = time_stepper_pt->ntstorage();

    // Increment p-order of the element
    P_order += inc;

    // Change integration scheme
    delete this->integral_pt();
    switch (P_order)
    {
      case 2:
        this->set_integration_scheme(new GaussLobattoLegendre<1, 2>);
        break;
      case 3:
        this->set_integration_scheme(new GaussLobattoLegendre<1, 3>);
        break;
      case 4:
        this->set_integration_scheme(new GaussLobattoLegendre<1, 4>);
        break;
      case 5:
        this->set_integration_scheme(new GaussLobattoLegendre<1, 5>);
        break;
      case 6:
        this->set_integration_scheme(new GaussLobattoLegendre<1, 6>);
        break;
      case 7:
        this->set_integration_scheme(new GaussLobattoLegendre<1, 7>);
        break;
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       integration scheme.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }

    // Allocate new space for Nodes (at the element level)
    this->set_n_node(P_order);

    // Copy vertex nodes and create new internal nodes
    //------------------------------------------------

    // Setup vertex coordinates in element:
    //-------------------------------------
    Vector<double> s_left(1);
    Vector<double> s_right(1);
    s_left[0] = -1.0;
    s_right[0] = 1.0;

    // Local coordinate in element
    Vector<double> s(1);

    Vector<double> s_fraction(1);
    // Loop over all nodes in the element
    for (unsigned n = 0; n < P_order; n++)
    {
      // Get the fractional position (in the current element) of the node
      // in the direction of s[0]
      s_fraction[0] = this->local_one_d_fraction_of_node(n, 0);

      // Evaluate the local coordinate of the node in the father element
      s[0] = s_left[0] + (s_right[0] - s_left[0]) * s_fraction[0];

      // Initialise flag: So far, this node hasn't been built or copied yet
      bool node_done = false;

      // Get the pointer to the node in this element (or rather, its clone),
      // returns NULL if there is not node
      Node* created_node_pt = clone_el_pt->get_node_at_local_coordinate(s);

      // Does this node already exist in this element?
      //----------------------------------------------
      if (created_node_pt != 0)
      {
        // Copy node across
        this->node_pt(n) = created_node_pt;

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
          clone_el_pt->get_interpolated_values(t, s, prev_values);
          // Find the minimum number of values
          //(either those stored at the node, or those returned by
          // the function)
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

      // Node does not exist in this element
      //------------------------------------

      // If the node has not been built anywhere ---> build it here
      if (!node_done)
      {
        // In a 1D mesh any node which lies on the boundary must exist in
        // the initial (coarse) mesh, so any newly-built nodes cannot be
        // boundary nodes. Therefore we always create a normal "bulk" node.

        // Create node and set the pointer to it from the element
        created_node_pt = this->construct_node(n, time_stepper_pt);

        // Now we set the position and values at the newly created node

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
          clone_el_pt->get_x(t, s, x_prev);

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
          clone_el_pt->get_interpolated_values(t, s, prev_values);

          // Determine the number of values at the new node
          const unsigned n_value = created_node_pt->nvalue();

          // Loop over all values and set the previous values
          for (unsigned v = 0; v < n_value; v++)
          {
            created_node_pt->set_value(t, v, prev_values[v]);
          }
        } // End of loop over history values

        // Add new node to mesh (if requested)
        if (mesh_pt != 0)
        {
          mesh_pt->add_node_pt(created_node_pt);
        }

        AlgebraicElementBase* alg_el_pt =
          dynamic_cast<AlgebraicElementBase*>(this);

        // If we do have an algebraic element
        if (alg_el_pt != 0)
        {
          std::string error_message = "Have not implemented p-refinement for";
          error_message += "Algebraic p-refineable elements yet\n";

          throw OomphLibError(
            error_message,
            "PRefineableQElement<1,INITIAL_NNODE_1D>::p_refine()",
            OOMPH_EXCEPTION_LOCATION);
        }

      } // End of case where we build the node ourselves

      // Check if the element is an algebraic element
      AlgebraicElementBase* alg_el_pt =
        dynamic_cast<AlgebraicElementBase*>(this);

      // If the element is an algebraic element, setup
      // node position (past and present) from algebraic node update
      // function. This over-writes previous assingments that
      // were made based on the macro-element/FE representation.
      // NOTE: YES, THIS NEEDS TO BE CALLED REPEATEDLY IF THE
      // NODE IS MEMBER OF MULTIPLE ELEMENTS: THEY ALL ASSIGN
      // THE SAME NODAL POSITIONS BUT WE NEED TO ADD THE REMESH
      // INFO FOR *ALL* ROOT ELEMENTS!
      if (alg_el_pt != 0)
      {
        // Build algebraic node update info for new node
        // This sets up the node update data for all node update
        // functions that are shared by all nodes in the father
        // element
        alg_el_pt->setup_algebraic_node_update(
          this->node_pt(n), s, clone_el_pt);
      }

    } // End of loop over all nodes in element


    // If the element is a MacroElementNodeUpdateElement, set
    // the update parameters for the current element's nodes --
    // all this needs is the vector of (pointers to the)
    // geometric objects that affect the MacroElement-based
    // node update -- this needs to be done to set the node
    // update info for newly created nodes
    MacroElementNodeUpdateElementBase* clone_m_el_pt =
      dynamic_cast<MacroElementNodeUpdateElementBase*>(clone_el_pt);
    if (clone_m_el_pt != 0)
    {
      // Get vector of geometric objects from father (construct vector
      // via copy operation)
      Vector<GeomObject*> geom_object_pt(clone_m_el_pt->geom_object_pt());

      // Cast current element to MacroElementNodeUpdateElement:
      MacroElementNodeUpdateElementBase* m_el_pt =
        dynamic_cast<MacroElementNodeUpdateElementBase*>(this);

#ifdef PARANOID
      if (m_el_pt == 0)
      {
        std::string error_message =
          "Failed to cast to MacroElementNodeUpdateElementBase*\n";
        error_message +=
          "Strange -- if my clone is a MacroElementNodeUpdateElement\n";
        error_message += "then I should be too....\n";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // Build update info by passing vector of geometric objects:
      // This sets the current element to be the update element
      // for all of the element's nodes -- this is reversed
      // if the element is ever un-refined in the father element's
      // rebuild_from_sons() function which overwrites this
      // assignment to avoid nasty segmentation faults that occur
      // when a node tries to update itself via an element that no
      // longer exists...
      m_el_pt->set_node_update_info(geom_object_pt);
    }


    // Loop over all nodes in element again, to re-set the positions
    // This must be done using the new element's macro-element
    // representation, rather than the old version which may be
    // of a different p-order!
    for (unsigned n = 0; n < P_order; n++)
    {
      // Get the fractional position of the node in the direction of s[0]
      s_fraction[0] = this->local_one_d_fraction_of_node(n, 0);
      // Local coordinate
      s[0] = s_left[0] + (s_right[0] - s_left[0]) * s_fraction[0];

      // Loop over # of history values
      for (unsigned t = 0; t < ntstorage; t++)
      {
        // Get position from father element -- this uses the macro
        // element representation if appropriate. If the node
        // turns out to be a hanging node later on, then
        // its position gets adjusted in line with its
        // hanging node interpolation.
        Vector<double> x_prev(1);
        this->get_x(t, s, x_prev);

        // Set previous positions of the new node
        this->node_pt(n)->x(t, 0) = x_prev[0];
      }
    }

    // Not necessary to delete the old nodes since all original nodes are in the
    // current mesh and so will be pruned as part of the mesh adaption process.

    // Do any further-build required
    this->further_build();
  }

  //=======================================================================
  /// Shape functions for PRefineableQElement<DIM>
  //=======================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<1, INITIAL_NNODE_1D>::shape(const Vector<double>& s,
                                                       Shape& psi) const
  {
    switch (p_order())
    {
      case 2:
      {
        // Calculate nodal positions
        OneDimensionalLegendreShape<2>::calculate_nodal_positions();
        // Create one-dim shape functions
        OneDimensionalLegendreShape<2> psi1(s[0]);
        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < p_order(); i++)
        {
          psi(i) = psi1[i];
        }
        break;
      }
      case 3:
      {
        // Calculate nodal positions
        OneDimensionalLegendreShape<3>::calculate_nodal_positions();
        // Create one-dim shape functions
        OneDimensionalLegendreShape<3> psi1(s[0]);
        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < p_order(); i++)
        {
          psi(i) = psi1[i];
        }
        break;
      }
      case 4:
      {
        // Calculate nodal positions
        OneDimensionalLegendreShape<4>::calculate_nodal_positions();
        // Create one-dim shape functions
        OneDimensionalLegendreShape<4> psi1(s[0]);
        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < p_order(); i++)
        {
          psi(i) = psi1[i];
        }
        break;
      }
      case 5:
      {
        // Calculate nodal positions
        OneDimensionalLegendreShape<5>::calculate_nodal_positions();
        // Create one-dim shape functions
        OneDimensionalLegendreShape<5> psi1(s[0]);
        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < p_order(); i++)
        {
          psi(i) = psi1[i];
        }
        break;
      }
      case 6:
      {
        // Calculate nodal positions
        OneDimensionalLegendreShape<6>::calculate_nodal_positions();
        // Create one-dim shape functions
        OneDimensionalLegendreShape<6> psi1(s[0]);
        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < p_order(); i++)
        {
          psi(i) = psi1[i];
        }
        break;
      }
      case 7:
      {
        // Calculate nodal positions
        OneDimensionalLegendreShape<7>::calculate_nodal_positions();
        // Create one-dim shape functions
        OneDimensionalLegendreShape<7> psi1(s[0]);
        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < p_order(); i++)
        {
          psi(i) = psi1[i];
        }
        break;
      }
      default:
        oomph_info << "\n ERROR: PRefineableQElement::shape() exceeded maximum";
        oomph_info << "\n        polynomial order for shape functions."
                   << std::endl;
    }
  }

  //=======================================================================
  /// Derivatives of shape functions for PRefineableQElement<DIM>
  //=======================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<1, INITIAL_NNODE_1D>::dshape_local(
    const Vector<double>& s, Shape& psi, DShape& dpsi) const
  {
    switch (p_order())
    {
      case 2:
      {
        // Calculate nodal positions
        OneDimensionalLegendreShape<2>::calculate_nodal_positions();
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<2> psi1(s[0]);
        OneDimensionalLegendreDShape<2> dpsi1ds(s[0]);
        // Loop over shapes and copy across
        for (unsigned i = 0; i < p_order(); i++)
        {
          psi(i) = psi1[i];
          dpsi(i, 0) = dpsi1ds[i];
        }

        break;
      }
      case 3:
      {
        // Calculate nodal positions
        OneDimensionalLegendreShape<3>::calculate_nodal_positions();
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<3> psi1(s[0]);
        OneDimensionalLegendreDShape<3> dpsi1ds(s[0]);
        // Loop over shapes and copy across
        for (unsigned i = 0; i < p_order(); i++)
        {
          psi(i) = psi1[i];
          dpsi(i, 0) = dpsi1ds[i];
        }
        break;
      }
      case 4:
      {
        // Calculate nodal positions
        OneDimensionalLegendreShape<4>::calculate_nodal_positions();
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<4> psi1(s[0]);
        OneDimensionalLegendreDShape<4> dpsi1ds(s[0]);
        // Loop over shapes and copy across
        for (unsigned i = 0; i < p_order(); i++)
        {
          psi(i) = psi1[i];
          dpsi(i, 0) = dpsi1ds[i];
        }
        break;
      }
      case 5:
      {
        // Calculate nodal positions
        OneDimensionalLegendreShape<5>::calculate_nodal_positions();
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<5> psi1(s[0]);
        OneDimensionalLegendreDShape<5> dpsi1ds(s[0]);
        // Loop over shapes and copy across
        for (unsigned i = 0; i < p_order(); i++)
        {
          psi(i) = psi1[i];
          dpsi(i, 0) = dpsi1ds[i];
        }
        break;
      }
      case 6:
      {
        // Calculate nodal positions
        OneDimensionalLegendreShape<6>::calculate_nodal_positions();
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<6> psi1(s[0]);
        OneDimensionalLegendreDShape<6> dpsi1ds(s[0]);
        // Loop over shapes and copy across
        for (unsigned i = 0; i < p_order(); i++)
        {
          psi(i) = psi1[i];
          dpsi(i, 0) = dpsi1ds[i];
        }
        break;
      }
      case 7:
      {
        // Calculate nodal positions
        OneDimensionalLegendreShape<7>::calculate_nodal_positions();
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<7> psi1(s[0]);
        OneDimensionalLegendreDShape<7> dpsi1ds(s[0]);
        // Loop over shapes and copy across
        for (unsigned i = 0; i < p_order(); i++)
        {
          psi(i) = psi1[i];
          dpsi(i, 0) = dpsi1ds[i];
        }
        break;
      }
      default:
        oomph_info
          << "\n ERROR: PRefineableQElement::dshape_local() exceeded maximum";
        oomph_info << "\n        polynomial order for shape functions."
                   << std::endl;
    }
  }

  //=======================================================================
  /// Second derivatives of shape functions for PRefineableQElement<DIM>
  ///  d2psids(i,0) = \f$ d^2 \psi_j / d s^2 \f$
  //=======================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<1, INITIAL_NNODE_1D>::d2shape_local(
    const Vector<double>& s, Shape& psi, DShape& dpsids, DShape& d2psids) const
  {
    std::ostringstream error_message;
    error_message
      << "\nd2shape_local currently not implemented for this element\n";
    throw OomphLibError(
      error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  //=================================================================
  /// Internal function to set up the hanging nodes on a particular
  /// edge of the element. (Not required in 1D.)
  //=================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<1, INITIAL_NNODE_1D>::binary_hang_helper(
    const int& value_id, const int& my_edge, std::ofstream& output_hangfile)
  {
  }

  //=======================================================================
  /// Rebuild the element from nodes found in its sons
  /// Adjusts its p-order to be the maximum of its sons' p-orders
  //=======================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<1, INITIAL_NNODE_1D>::rebuild_from_sons(
    Mesh*& mesh_pt)
  {
    // Get p-orders of sons
    unsigned n_sons = this->tree_pt()->nsons();
    Vector<unsigned> son_p_order(n_sons);
    unsigned max_son_p_order = 0;
    for (unsigned ison = 0; ison < n_sons; ison++)
    {
      PRefineableElement* el_pt = dynamic_cast<PRefineableElement*>(
        this->tree_pt()->son_pt(ison)->object_pt());
      son_p_order[ison] = el_pt->p_order();
      if (son_p_order[ison] > max_son_p_order)
        max_son_p_order = son_p_order[ison];
    }

    unsigned old_Nnode = this->nnode();
    unsigned old_P_order = this->p_order();
    // Set p-order of the element
    this->p_order() = max_son_p_order;

    // Change integration scheme
    delete this->integral_pt();
    switch (this->p_order())
    {
      case 2:
        this->set_integration_scheme(new GaussLobattoLegendre<1, 2>);
        break;
      case 3:
        this->set_integration_scheme(new GaussLobattoLegendre<1, 3>);
        break;
      case 4:
        this->set_integration_scheme(new GaussLobattoLegendre<1, 4>);
        break;
      case 5:
        this->set_integration_scheme(new GaussLobattoLegendre<1, 5>);
        break;
      case 6:
        this->set_integration_scheme(new GaussLobattoLegendre<1, 6>);
        break;
      case 7:
        this->set_integration_scheme(new GaussLobattoLegendre<1, 7>);
        break;
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       integration scheme.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }

    // Back up pointers to old nodes before they are lost
    Vector<Node*> old_node_pt(old_Nnode);
    for (unsigned n = 0; n < old_Nnode; n++)
    {
      old_node_pt[n] = this->node_pt(n);
    }

    // Allocate new space for Nodes (at the element level)
    this->set_n_node(this->p_order());

    // Copy vertex nodes and create new edge and internal nodes
    //---------------------------------------------------------

    // Copy vertex nodes
    this->node_pt(0) = old_node_pt[0];
    this->node_pt(this->p_order() - 1) = old_node_pt[old_P_order - 1];


    //=============================================================
    // Below this line is copied from RefineableQSpectralElement<2>

    // The timestepper should be the same for all nodes and node 0 should
    // never be deleted.
    if (this->node_pt(0) == 0)
    {
      throw OomphLibError("The vertex node (0) does not exist",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    TimeStepper* time_stepper_pt = this->node_pt(0)->time_stepper_pt();

    // Determine number of history values stored
    const unsigned ntstorage = time_stepper_pt->ntstorage();

    // Allocate storage for local coordinates
    Vector<double> s_fraction(1), s(1);

    // Determine the number of nodes in the element
    const unsigned n_node = this->nnode_1d();

    // Loop over the nodes in the element
    for (unsigned n = 0; n < n_node; n++)
    {
      // Get the fractional position of the node in the direction of s[0]
      s_fraction[0] = this->local_one_d_fraction_of_node(n, 0);

      // Determine the local coordinate in the father element
      s[0] = -1.0 + 2.0 * s_fraction[0];

      // If the node has not been built
      if (this->node_pt(n) == 0)
      {
        // Has the node been created by one of its neighbours?
        bool is_periodic = false;
        Node* created_node_pt =
          this->node_created_by_neighbour(s_fraction, is_periodic);

        // If it has, set the pointer
        if (created_node_pt != 0)
        {
          // If the node is periodic
          if (is_periodic)
          {
            throw OomphLibError("Cannot handle periodic nodes yet",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
          // Non-periodic case, just set the pointer
          else
          {
            this->node_pt(n) = created_node_pt;
          }
        }
        // Otherwise, we need to build it
        else
        {
          // First, find the son element in which the node should live

          // Find coordinate in the son
          Vector<double> s_in_son(1);
          using namespace BinaryTreeNames;
          int son = -10;
          // If s_fraction is between 0 and 0.5, we are in the left son
          if (s_fraction[0] < 0.5)
          {
            son = L;
            s_in_son[0] = -1.0 + 4.0 * s_fraction[0];
          }
          // Otherwise we are in the right son
          else
          {
            son = R;
            s_in_son[0] = -1.0 + 4.0 * (s_fraction[0] - 0.5);
          }

          // Get the pointer to the son element in which the new node
          // would live
          PRefineableQElement<1, INITIAL_NNODE_1D>* son_el_pt =
            dynamic_cast<PRefineableQElement<1, INITIAL_NNODE_1D>*>(
              this->tree_pt()->son_pt(son)->object_pt());

          // In 1D we should never be rebuilding an element's vertex nodes
          // (since they will never be deleted), so throw an error if we
          // appear to be doing so
#ifdef PARANOID
          if (n == 0 || n == n_node - 1)
          {
            std::string error_message =
              "I am trying to rebuild one of the (two) vertex nodes in\n";
            error_message +=
              "this 1D element. It should not have been possible to delete\n";
            error_message += "either of these!\n";

            throw OomphLibError(
              error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
          }
#endif

          // With this in mind we will always be creating normal "bulk" nodes
          this->node_pt(n) = this->construct_node(n, time_stepper_pt);

          // Now we set the position and values at the newly created node

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

            // Get position from son element -- this uses the macro element
            // representation if appropriate
            son_el_pt->get_x(t, s_in_son, x_prev);

            // Set the previous position of the new node
            this->node_pt(n)->x(t, 0) = x_prev[0];

            // Allocate storage for the previous values at the node
            // NOTE: the size of this vector is equal to the number of values
            // (e.g. 3 velocity components and 1 pressure, say)
            // associated with each node and NOT the number of history values
            // which the node stores!
            Vector<double> prev_values;

            // Get values from son element
            // Note: get_interpolated_values() sets Vector size itself.
            son_el_pt->get_interpolated_values(t, s_in_son, prev_values);

            // Determine the number of values at the new node
            const unsigned n_value = this->node_pt(n)->nvalue();

            // Loop over all values and set the previous values
            for (unsigned v = 0; v < n_value; v++)
            {
              this->node_pt(n)->set_value(t, v, prev_values[v]);
            }
          } // End of loop over history values

          // Add new node to mesh
          mesh_pt->add_node_pt(this->node_pt(n));

        } // End of case where we build the node ourselves

      } // End of if this node has not been built
    } // End of loop over nodes in element

    // Check if the element is an algebraic element
    // This is done on all nodes in the element after reconstruction
    // rather than as the nodes are built
    AlgebraicElementBase* alg_el_pt = dynamic_cast<AlgebraicElementBase*>(this);

    // If so, throw error
    if (alg_el_pt != 0)
    {
      std::string error_message =
        "Have not implemented rebuilding from sons for";
      error_message += "Algebraic p-refineable elements yet\n";

      throw OomphLibError(
        error_message,
        "PRefineableQElement<1,INITIAL_NNODE_1D>::rebuild_from_sons()",
        OOMPH_EXCEPTION_LOCATION);
    }
  }

  //=================================================================
  /// Check inter-element continuity of
  /// - nodal positions
  /// - (nodally) interpolated function values
  //====================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<1, INITIAL_NNODE_1D>::check_integrity(
    double& max_error)
  {
    RefineableQElement<1>::check_integrity(max_error);
  }

  /// /////////////////////////////////////////////////////////////
  //       2D PRefineableQElements
  /// /////////////////////////////////////////////////////////////

  /// Get local coordinates of node j in the element; vector sets its own size
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<2, INITIAL_NNODE_1D>::local_coordinate_of_node(
    const unsigned& n, Vector<double>& s) const
  {
    s.resize(2);
    unsigned Nnode_1d = this->nnode_1d();
    unsigned n0 = n % Nnode_1d;
    unsigned n1 = n / Nnode_1d;

    switch (Nnode_1d)
    {
      case 2:
        OneDimensionalLegendreShape<2>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<2>::nodal_position(n0);
        s[1] = OneDimensionalLegendreShape<2>::nodal_position(n1);
        break;
      case 3:
        OneDimensionalLegendreShape<3>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<3>::nodal_position(n0);
        s[1] = OneDimensionalLegendreShape<3>::nodal_position(n1);
        break;
      case 4:
        OneDimensionalLegendreShape<4>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<4>::nodal_position(n0);
        s[1] = OneDimensionalLegendreShape<4>::nodal_position(n1);
        break;
      case 5:
        OneDimensionalLegendreShape<5>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<5>::nodal_position(n0);
        s[1] = OneDimensionalLegendreShape<5>::nodal_position(n1);
        break;
      case 6:
        OneDimensionalLegendreShape<6>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<6>::nodal_position(n0);
        s[1] = OneDimensionalLegendreShape<6>::nodal_position(n1);
        break;
      case 7:
        OneDimensionalLegendreShape<7>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<7>::nodal_position(n0);
        s[1] = OneDimensionalLegendreShape<7>::nodal_position(n1);
        break;
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       shape functions.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        break;
    }
  }

  /// Get the local fractino of node j in the element
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<2, INITIAL_NNODE_1D>::local_fraction_of_node(
    const unsigned& n, Vector<double>& s_fraction)
  {
    this->local_coordinate_of_node(n, s_fraction);
    s_fraction[0] = 0.5 * (s_fraction[0] + 1.0);
    s_fraction[1] = 0.5 * (s_fraction[1] + 1.0);
  }

  template<unsigned INITIAL_NNODE_1D>
  double PRefineableQElement<2, INITIAL_NNODE_1D>::local_one_d_fraction_of_node(
    const unsigned& n1d, const unsigned& i)
  {
    switch (this->nnode_1d())
    {
      case 2:
        OneDimensionalLegendreShape<2>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<2>::nodal_position(n1d) + 1.0);
      case 3:
        OneDimensionalLegendreShape<3>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<3>::nodal_position(n1d) + 1.0);
      case 4:
        OneDimensionalLegendreShape<4>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<4>::nodal_position(n1d) + 1.0);
      case 5:
        OneDimensionalLegendreShape<5>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<5>::nodal_position(n1d) + 1.0);
      case 6:
        OneDimensionalLegendreShape<6>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<6>::nodal_position(n1d) + 1.0);
      case 7:
        OneDimensionalLegendreShape<7>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<7>::nodal_position(n1d) + 1.0);
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       shape functions.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        return 0.0;
    }
  }


  //==================================================================
  /// Return the node at the specified local coordinate
  //==================================================================
  template<unsigned INITIAL_NNODE_1D>
  Node* PRefineableQElement<2, INITIAL_NNODE_1D>::get_node_at_local_coordinate(
    const Vector<double>& s) const
  {
    unsigned Nnode_1d = this->nnode_1d();
    // Load the tolerance into a local variable
    double tol = FiniteElement::Node_location_tolerance;
    // There are two possible indices.
    Vector<int> index(2);

    // Loop over indices
    for (unsigned i = 0; i < 2; i++)
    {
      // Determine the index
      // -------------------

      bool is_found = false;

      // If we are at the lower limit, the index is zero
      if (std::fabs(s[i] + 1.0) < tol)
      {
        index[i] = 0;
        is_found = true;
      }
      // If we are at the upper limit, the index is the number of nodes minus 1
      else if (std::fabs(s[i] - 1.0) < tol)
      {
        index[i] = Nnode_1d - 1;
        is_found = true;
      }
      // Otherwise, we have to calculate the index in general
      else
      {
        // Compute Gauss-Lobatto-Legendre node positions
        Vector<double> z;
        Orthpoly::gll_nodes(Nnode_1d, z);
        // Loop over possible internal nodal positions
        for (unsigned n = 1; n < Nnode_1d - 1; n++)
        {
          if (std::fabs(z[n] - s[i]) < tol)
          {
            index[i] = n;
            is_found = true;
            break;
          }
        }
      }

      if (!is_found)
      {
        // No matching nodes
        return 0;
      }
    }
    // If we've got here we have a node, so let's return a pointer to it
    return this->node_pt(index[0] + Nnode_1d * index[1]);
  }

  //===================================================================
  /// If a neighbouring element has already created a node at
  /// a position corresponding to the local fractional position within the
  /// present element, s_fraction, return
  /// a pointer to that node. If not, return NULL (0). If the node is
  /// periodic the flag is_periodic will be true
  //===================================================================
  template<unsigned INITIAL_NNODE_1D>
  Node* PRefineableQElement<2, INITIAL_NNODE_1D>::node_created_by_neighbour(
    const Vector<double>& s_fraction, bool& is_periodic)
  {
    using namespace QuadTreeNames;

    // Calculate the edges on which the node lies
    Vector<int> edges;
    if (s_fraction[0] == 0.0)
    {
      edges.push_back(W);
    }
    if (s_fraction[0] == 1.0)
    {
      edges.push_back(E);
    }
    if (s_fraction[1] == 0.0)
    {
      edges.push_back(S);
    }
    if (s_fraction[1] == 1.0)
    {
      edges.push_back(N);
    }

    // Find the number of edges
    unsigned n_size = edges.size();
    // If there are no edges, then there is no neighbour, return 0
    if (n_size == 0)
    {
      return 0;
    }

    Vector<unsigned> translate_s(2);
    Vector<double> s_lo_neigh(2);
    Vector<double> s_hi_neigh(2);
    Vector<double> s(2);

    int neigh_edge, diff_level;
    bool in_neighbouring_tree;

    // Loop over the edges
    for (unsigned j = 0; j < n_size; j++)
    {
      // Find pointer to neighbouring element along edge
      QuadTree* neigh_pt;
      neigh_pt = quadtree_pt()->gteq_edge_neighbour(edges[j],
                                                    translate_s,
                                                    s_lo_neigh,
                                                    s_hi_neigh,
                                                    neigh_edge,
                                                    diff_level,
                                                    in_neighbouring_tree);

      // Neighbour exists
      if (neigh_pt != 0)
      {
        // Have any of its vertex nodes been created yet?
        // (Must look in incomplete neighbours because after the
        // pre-build they may contain pointers to the required nodes. e.g.
        // h-refinement of neighbouring linear and quadratic elements)
        bool a_vertex_node_is_built = false;
        QElement<2, INITIAL_NNODE_1D>* neigh_obj_pt =
          dynamic_cast<QElement<2, INITIAL_NNODE_1D>*>(neigh_pt->object_pt());
        if (neigh_obj_pt == 0)
        {
          throw OomphLibError("Not a quad element!",
                              "PRefineableQElement<2,INITIAL_NNODE_1D>::node_"
                              "created_by_neighbour()",
                              OOMPH_EXCEPTION_LOCATION);
        }
        for (unsigned vnode = 0; vnode < neigh_obj_pt->nvertex_node(); vnode++)
        {
          if (neigh_obj_pt->vertex_node_pt(vnode) != 0)
            a_vertex_node_is_built = true;
          break;
        }
        if (a_vertex_node_is_built)
        {
          // We now need to translate the nodal location
          // as defined in terms of the fractional coordinates of the present
          // element into those of its neighbour

          // Calculate the local coordinate in the neighbour
          // Note that we need to use the translation scheme in case
          // the local coordinates are swapped in the neighbour.
          for (unsigned i = 0; i < 2; i++)
          {
            s[i] = s_lo_neigh[i] +
                   s_fraction[translate_s[i]] * (s_hi_neigh[i] - s_lo_neigh[i]);
          }

          // Find the node in the neighbour
          Node* neighbour_node_pt =
            neigh_pt->object_pt()->get_node_at_local_coordinate(s);

          // If there is a node, return it
          if (neighbour_node_pt != 0)
          {
            // Now work out whether it's a periodic boundary
            // only possible if we have moved into a neighbouring tree
            if (in_neighbouring_tree)
            {
              // Return whether the neighbour is periodic
              is_periodic =
                quadtree_pt()->root_pt()->is_neighbour_periodic(edges[j]);
            }
            // Return the pointer to the neighbouring node
            return neighbour_node_pt;
          }
        }
      }
    }
    // Node not found, return null
    return 0;
  }

  //===================================================================
  /// If a neighbouring element's son has already created a node at
  /// a position corresponding to the local fractional position within the
  /// present element, s_fraction, return
  /// a pointer to that node. If not, return NULL (0). If the node is
  /// periodic the flag is_periodic will be true
  //===================================================================
  template<unsigned INITIAL_NNODE_1D>
  Node* PRefineableQElement<2, INITIAL_NNODE_1D>::
    node_created_by_son_of_neighbour(const Vector<double>& s_fraction,
                                     bool& is_periodic)
  {
    using namespace QuadTreeNames;

    // Calculate the edges on which the node lies
    Vector<int> edges;
    if (s_fraction[0] == 0.0)
    {
      edges.push_back(W);
    }
    if (s_fraction[0] == 1.0)
    {
      edges.push_back(E);
    }
    if (s_fraction[1] == 0.0)
    {
      edges.push_back(S);
    }
    if (s_fraction[1] == 1.0)
    {
      edges.push_back(N);
    }

    // Find the number of edges
    unsigned n_size = edges.size();
    // If there are no edges, then there is no neighbour, return 0
    if (n_size == 0)
    {
      return 0;
    }

    Vector<unsigned> translate_s(2);
    Vector<double> s_lo_neigh(2);
    Vector<double> s_hi_neigh(2);
    Vector<double> s(2);

    int neigh_edge, diff_level;
    bool in_neighbouring_tree;

    // Loop over the edges
    for (unsigned j = 0; j < n_size; j++)
    {
      // Find pointer to neighbouring element along edge
      QuadTree* neigh_pt;
      neigh_pt = quadtree_pt()->gteq_edge_neighbour(edges[j],
                                                    translate_s,
                                                    s_lo_neigh,
                                                    s_hi_neigh,
                                                    neigh_edge,
                                                    diff_level,
                                                    in_neighbouring_tree);

      // Neighbour exists
      if (neigh_pt != 0)
      {
        // Have its nodes been created yet?
        // (Must look in sons of unfinished neighbours too!!!)
        if (1)
        {
          // We now need to translate the nodal location
          // as defined in terms of the fractional coordinates of the present
          // element into those of its neighbour

          // Calculate the local coordinate in the neighbour
          // Note that we need to use the translation scheme in case
          // the local coordinates are swapped in the neighbour.
          for (unsigned i = 0; i < 2; i++)
          {
            s[i] = s_lo_neigh[i] +
                   s_fraction[translate_s[i]] * (s_hi_neigh[i] - s_lo_neigh[i]);
          }

          // Check if the element has sons
          if (neigh_pt->nsons() != 0)
          {
            // First, find the son element in which the node should live

            // Find coordinates in the sons
            Vector<double> s_in_son(2);
            int son = -10;
            // If negative on the west side
            if (s[0] < 0.0)
            {
              // On the south side
              if (s[1] < 0.0)
              {
                // It's the southwest son
                son = SW;
                s_in_son[0] = 1.0 + 2.0 * s[0];
                s_in_son[1] = 1.0 + 2.0 * s[1];
              }
              // On the north side
              else
              {
                // It's the northwest son
                son = NW;
                s_in_son[0] = 1.0 + 2.0 * s[0];
                s_in_son[1] = -1.0 + 2.0 * s[1];
              }
            }
            else
            {
              // On the south side
              if (s[1] < 0.0)
              {
                // It's the southeast son
                son = SE;
                s_in_son[0] = -1.0 + 2.0 * s[0];
                s_in_son[1] = 1.0 + 2.0 * s[1];
              }
              // On the north side
              else
              {
                // It's the northeast son
                son = NE;
                s_in_son[0] = -1.0 + 2.0 * s[0];
                s_in_son[1] = -1.0 + 2.0 * s[1];
              }
            }

            // Find the node in the neighbour's son
            Node* neighbour_son_node_pt =
              neigh_pt->son_pt(son)->object_pt()->get_node_at_local_coordinate(
                s_in_son);

            // If there is a node, return it
            if (neighbour_son_node_pt != 0)
            {
              // Now work out whether it's a periodic boundary
              // only possible if we have moved into a neighbouring tree
              if (in_neighbouring_tree)
              {
                // Return whether the neighbour is periodic
                is_periodic =
                  quadtree_pt()->root_pt()->is_neighbour_periodic(edges[j]);
              }
              // Return the pointer to the neighbouring node
              return neighbour_son_node_pt;
            }
          }
          else
          {
            // No sons to search in, so no node can be found
            return 0;
          }
        }
      }
    }
    // Node not found, return null
    return 0;
  }

  //==================================================================
  /// Set the correct p-order of the element based on that of its
  /// father. Then construct an integration scheme of the correct order.
  /// If an adopted father is specified, information from this is
  /// used instead of using the father found from the tree.
  //==================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<2, INITIAL_NNODE_1D>::initial_setup(
    Tree* const& adopted_father_pt, const unsigned& initial_p_order)
  {
    // Storage for pointer to my father (in quadtree impersonation)
    QuadTree* father_pt = 0;

    // Check if an adopted father has been specified
    if (adopted_father_pt != 0)
    {
      // Get pointer to my father (in quadtree impersonation)
      father_pt = dynamic_cast<QuadTree*>(adopted_father_pt);
    }
    // Check if element is in a tree
    else if (Tree_pt != 0)
    {
      // Get pointer to my father (in quadtree impersonation)
      father_pt = dynamic_cast<QuadTree*>(quadtree_pt()->father_pt());
    }
    // else
    // {
    //  throw OomphLibError(
    //         "Element not in a tree, and no adopted father has been
    //         specified!", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    // }

    // Check if we found a father
    if (father_pt != 0 || initial_p_order != 0)
    {
      if (father_pt != 0)
      {
        PRefineableQElement<2, INITIAL_NNODE_1D>* father_el_pt =
          dynamic_cast<PRefineableQElement<2, INITIAL_NNODE_1D>*>(
            father_pt->object_pt());
        if (father_el_pt != 0)
        {
          unsigned father_p_order = father_el_pt->p_order();
          // Set p-order to that of father
          P_order = father_p_order;
        }
      }
      else
      {
        P_order = initial_p_order;
      }

      // Now sort out the element...
      // (has p^2 nodes)
      unsigned new_n_node = P_order * P_order;

      // Allocate new space for Nodes (at the element level)
      this->set_n_node(new_n_node);

      // Set integration scheme
      delete this->integral_pt();
      switch (P_order)
      {
        case 2:
          this->set_integration_scheme(new GaussLobattoLegendre<2, 2>);
          break;
        case 3:
          this->set_integration_scheme(new GaussLobattoLegendre<2, 3>);
          break;
        case 4:
          this->set_integration_scheme(new GaussLobattoLegendre<2, 4>);
          break;
        case 5:
          this->set_integration_scheme(new GaussLobattoLegendre<2, 5>);
          break;
        case 6:
          this->set_integration_scheme(new GaussLobattoLegendre<2, 6>);
          break;
        case 7:
          this->set_integration_scheme(new GaussLobattoLegendre<2, 7>);
          break;
        default:
          std::ostringstream error_message;
          error_message << "\nERROR: Exceeded maximum polynomial order for";
          error_message << "\n       integration scheme.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }
  }

  //==================================================================
  /// Check the father element for any required nodes which
  /// already exist
  //==================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<2, INITIAL_NNODE_1D>::pre_build(
    Mesh*& mesh_pt, Vector<Node*>& new_node_pt)
  {
    using namespace QuadTreeNames;

    // Get the number of 1d nodes
    unsigned n_p = nnode_1d();

    // Check whether static father_bound needs to be created
    if (Father_bound[n_p].nrow() == 0)
    {
      setup_father_bounds();
    }

    // Pointer to my father (in quadtree impersonation)
    QuadTree* father_pt = dynamic_cast<QuadTree*>(quadtree_pt()->father_pt());

    // What type of son am I? Ask my quadtree representation...
    int son_type = Tree_pt->son_type();

    // Has somebody build me already? (If any nodes have not been built)
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

      // Return pointer to father element
      RefineableQElement<2>* father_el_pt =
        dynamic_cast<RefineableQElement<2>*>(father_pt->object_pt());

      // Timestepper should be the same for all nodes in father
      // element -- use it create timesteppers for new nodes
      TimeStepper* time_stepper_pt =
        father_el_pt->node_pt(0)->time_stepper_pt();

      // Number of history values (incl. present)
      unsigned ntstorage = time_stepper_pt->ntstorage();

      Vector<double> s_lo(2);
      Vector<double> s_hi(2);
      Vector<double> s(2);
      Vector<double> x(2);

      // Setup vertex coordinates in father element:
      //--------------------------------------------
      switch (son_type)
      {
        case SW:
          s_lo[0] = -1.0;
          s_hi[0] = 0.0;
          s_lo[1] = -1.0;
          s_hi[1] = 0.0;
          break;

        case SE:
          s_lo[0] = 0.0;
          s_hi[0] = 1.0;
          s_lo[1] = -1.0;
          s_hi[1] = 0.0;
          break;

        case NE:
          s_lo[0] = 0.0;
          s_hi[0] = 1.0;
          s_lo[1] = 0.0;
          s_hi[1] = 1.0;
          break;

        case NW:
          s_lo[0] = -1.0;
          s_hi[0] = 0.0;
          s_lo[1] = 0.0;
          s_hi[1] = 1.0;
          break;
      }

      // Pass macro element pointer on to sons and
      // set coordinates in macro element
      // hierher why can I see this?
      if (father_el_pt->macro_elem_pt() != 0)
      {
        set_macro_elem_pt(father_el_pt->macro_elem_pt());
        for (unsigned i = 0; i < 2; i++)
        {
          s_macro_ll(i) =
            father_el_pt->s_macro_ll(i) +
            0.5 * (s_lo[i] + 1.0) *
              (father_el_pt->s_macro_ur(i) - father_el_pt->s_macro_ll(i));
          s_macro_ur(i) =
            father_el_pt->s_macro_ll(i) +
            0.5 * (s_hi[i] + 1.0) *
              (father_el_pt->s_macro_ur(i) - father_el_pt->s_macro_ll(i));
        }
      }


      // If the father element hasn't been generated yet, we're stuck...
      if (father_el_pt->node_pt(0) == 0)
      {
        throw OomphLibError(
          "Trouble: father_el_pt->node_pt(0)==0\n Can't build son element!\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        unsigned jnod = 0;
        Vector<double> x_small(2);
        Vector<double> x_large(2);

        Vector<double> s_fraction(2);
        // Loop over nodes in element
        for (unsigned i0 = 0; i0 < n_p; i0++)
        {
          // Get the fractional position of the node in the direction of s[0]
          s_fraction[0] = this->local_one_d_fraction_of_node(i0, 0);
          // Local coordinate in father element
          s[0] = s_lo[0] + (s_hi[0] - s_lo[0]) * s_fraction[0];

          for (unsigned i1 = 0; i1 < n_p; i1++)
          {
            // Get the fractional position of the node in the direction of s[1]
            s_fraction[1] = this->local_one_d_fraction_of_node(i1, 1);
            // Local coordinate in father element
            s[1] = s_lo[1] + (s_hi[1] - s_lo[1]) * s_fraction[1];

            // Local node number
            jnod = i0 + n_p * i1;

            // Check whether the father's node is periodic if so, complain
            /* {
              Node* father_node_pt = father_el_pt->node_pt(jnod);
              if((father_node_pt->is_a_copy()) ||
              (father_node_pt->position_is_a_copy()))
              {
              throw OomphLibError(
              "Can't handle periodic nodes (yet).",
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
              }
              }*/

            // Initialise flag: So far, this node hasn't been built
            // or copied yet
            // bool node_done=false;

            // Get the pointer to the node in the father, returns NULL
            // if there is not node
            Node* created_node_pt =
              father_el_pt->get_node_at_local_coordinate(s);

            // Does this node already exist in father element?
            //------------------------------------------------
            if (created_node_pt != 0)
            {
              // Copy node across
              this->node_pt(jnod) = created_node_pt;

              // Make sure that we update the values at the node so that
              // they are consistent with the present representation.
              // This is only need for mixed interpolation where the value
              // at the father could now become active.

              // Loop over all history values
              for (unsigned t = 0; t < ntstorage; t++)
              {
                // Get values from father element
                // Note: get_interpolated_values() sets Vector size itself.
                Vector<double> prev_values;
                father_el_pt->get_interpolated_values(t, s, prev_values);
                // Find the minimum number of values
                //(either those stored at the node, or those returned by
                // the function)
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

              // Node has been created by copy
              // node_done=true;
            }
          } // End of vertical loop over nodes in element
        } // End of horizontal loop over nodes in element
      } // Sanity check: Father element has been generated

    } // End of nodes not built
    else
    {
      // If this is element updates by macro element node update then we need
      // to set the correct info in the nodes here since this won't get done
      // later in build() because we already have all our nodes created.
      MacroElementNodeUpdateElementBase* m_el_pt =
        dynamic_cast<MacroElementNodeUpdateElementBase*>(this);
      if (m_el_pt != 0)
      {
        // Get vector of geometric objects
        Vector<GeomObject*> geom_object_pt = m_el_pt->geom_object_pt();

        /// / Build update info by passing vector of geometric objects:
        /// / This sets the current element to be the update element
        /// / for all of the element's nodes -- this is reversed
        /// / if the element is ever un-refined in the father element's
        /// / rebuild_from_sons() function which overwrites this
        /// / assignment to avoid nasty segmentation faults that occur
        /// / when a node tries to update itself via an element that no
        /// / longer exists...
        m_el_pt->set_node_update_info(geom_object_pt);
      }
    }
  }

  //==================================================================
  /// p-refine the element inc times. (If inc<0 then p-unrefinement
  /// is performed.)
  //==================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<2, INITIAL_NNODE_1D>::p_refine(
    const int& inc, Mesh* const& mesh_pt, GeneralisedElement* const& clone_pt)
  {
    // Cast clone to correct type
    PRefineableQElement<2, INITIAL_NNODE_1D>* clone_el_pt =
      dynamic_cast<PRefineableQElement<2, INITIAL_NNODE_1D>*>(clone_pt);

    // If it is a MacroElement then we need to copy the geometric data too.
    MacroElementNodeUpdateElementBase* m_el_pt =
      dynamic_cast<MacroElementNodeUpdateElementBase*>(this);
    if (m_el_pt != 0)
    {
      MacroElementNodeUpdateElementBase* clone_m_el_pt =
        dynamic_cast<MacroElementNodeUpdateElementBase*>(clone_pt);

      // Store local copy of geom object vector, so it can be passed on
      // to son elements (and their nodes) during refinement
      unsigned ngeom_object = m_el_pt->ngeom_object();
      clone_m_el_pt->geom_object_pt().resize(ngeom_object);
      for (unsigned i = 0; i < ngeom_object; i++)
      {
        clone_m_el_pt->geom_object_pt()[i] = m_el_pt->geom_object_pt(i);
      }

      // clone_m_el_pt->set_node_update_info(m_el_pt->geom_object_pt());
    }

    // Check if we can use it
    if (clone_el_pt == 0)
    {
      throw OomphLibError(
        "Cloned copy must be a PRefineableQElement<2,INITIAL_NNODE_1D>!",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#ifdef PARANOID
    // Clone exists, so check that it is infact a copy of me
    else
    {
      // Flag to keep track of check
      bool clone_is_ok = true;

      // Does the clone have the correct p-order?
      clone_is_ok = clone_is_ok && (clone_el_pt->p_order() == this->p_order());

      if (!clone_is_ok)
      {
        std::ostringstream err_stream;
        err_stream << "Clone element has a different p-order from me,"
                   << std::endl
                   << "but it is supposed to be a copy!" << std::endl;

        throw OomphLibError(
          err_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Does the clone have the same number of nodes as me?
      clone_is_ok = clone_is_ok && (clone_el_pt->nnode() == this->nnode());

      if (!clone_is_ok)
      {
        std::ostringstream err_stream;
        err_stream << "Clone element has a different number of nodes from me,"
                   << std::endl
                   << "but it is supposed to be a copy!" << std::endl;

        throw OomphLibError(
          err_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Does the clone have the same nodes as me?
      for (unsigned n = 0; n < this->nnode(); n++)
      {
        clone_is_ok =
          clone_is_ok && (clone_el_pt->node_pt(n) == this->node_pt(n));
      }

      if (!clone_is_ok)
      {
        std::ostringstream err_stream;
        err_stream << "The nodes belonging to the clone element are different"
                   << std::endl
                   << "from mine, but it is supposed to be a copy!"
                   << std::endl;

        throw OomphLibError(
          err_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Is this a macro element?
      MacroElementNodeUpdateElementBase* m_el_pt =
        dynamic_cast<MacroElementNodeUpdateElementBase*>(this);
      if (m_el_pt != 0)
      {
        // Get macro element version of clone
        MacroElementNodeUpdateElementBase* clone_m_el_pt =
          dynamic_cast<MacroElementNodeUpdateElementBase*>(clone_el_pt);

        // Does the clone have the same geometric objects?
        Vector<GeomObject*> clone_geom_object_pt(
          clone_m_el_pt->geom_object_pt());
        Vector<GeomObject*> geom_object_pt(m_el_pt->geom_object_pt());

        clone_is_ok =
          clone_is_ok && (geom_object_pt.size() == clone_geom_object_pt.size());
        for (unsigned n = 0;
             n < std::min(geom_object_pt.size(), clone_geom_object_pt.size());
             n++)
        {
          clone_is_ok =
            clone_is_ok && (clone_geom_object_pt[n] == geom_object_pt[n]);
        }

        if (!clone_is_ok)
        {
          std::ostringstream err_stream;
          err_stream << "The clone element has different geometric objects"
                     << std::endl
                     << "from mine, but it is supposed to be a copy!"
                     << std::endl;

          throw OomphLibError(
            err_stream.str(),
            "PRefineableQElement<2,INITIAL_NNODE_1D>::p_refine()",
            OOMPH_EXCEPTION_LOCATION);
        }
      }

      // If we get to here then the clone has all the information we require
    }
#endif

    // Currently we can't handle the case of generalised coordinates
    // since we haven't established how they should be interpolated.
    // Buffer this case:
    if (clone_el_pt->node_pt(0)->nposition_type() != 1)
    {
      throw OomphLibError("Can't handle generalised nodal positions (yet).",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Timestepper should be the same for all nodes -- use it
    // to create timesteppers for new nodes
    TimeStepper* time_stepper_pt = this->node_pt(0)->time_stepper_pt();

    // Get number of history values (incl. present)
    unsigned ntstorage = time_stepper_pt->ntstorage();

    // Increment p-order of the element
    P_order += inc;

    // Change integration scheme
    delete this->integral_pt();
    switch (P_order)
    {
      case 2:
        this->set_integration_scheme(new GaussLobattoLegendre<2, 2>);
        break;
      case 3:
        this->set_integration_scheme(new GaussLobattoLegendre<2, 3>);
        break;
      case 4:
        this->set_integration_scheme(new GaussLobattoLegendre<2, 4>);
        break;
      case 5:
        this->set_integration_scheme(new GaussLobattoLegendre<2, 5>);
        break;
      case 6:
        this->set_integration_scheme(new GaussLobattoLegendre<2, 6>);
        break;
      case 7:
        this->set_integration_scheme(new GaussLobattoLegendre<2, 7>);
        break;
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       integration scheme.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }

    // Allocate new space for Nodes (at the element level)
    this->set_n_node(P_order * P_order);

    // Copy vertex nodes and create new edge and internal nodes
    //---------------------------------------------------------

    // Setup vertex coordinates in element:
    //-------------------------------------
    Vector<double> s_lo(2);
    Vector<double> s_hi(2);
    s_lo[0] = -1.0;
    s_hi[0] = 1.0;
    s_lo[1] = -1.0;
    s_hi[1] = 1.0;

    // Local coordinate in element
    Vector<double> s(2);

    unsigned jnod = 0;

    Vector<double> s_fraction(2);
    // Loop over nodes in element
    for (unsigned i0 = 0; i0 < P_order; i0++)
    {
      // Get the fractional position of the node in the direction of s[0]
      s_fraction[0] = this->local_one_d_fraction_of_node(i0, 0);
      // Local coordinate
      s[0] = s_lo[0] + (s_hi[0] - s_lo[0]) * s_fraction[0];

      for (unsigned i1 = 0; i1 < P_order; i1++)
      {
        // Get the fractional position of the node in the direction of s[1]
        s_fraction[1] = this->local_one_d_fraction_of_node(i1, 1);
        // Local coordinate
        s[1] = s_lo[1] + (s_hi[1] - s_lo[1]) * s_fraction[1];

        // Local node number
        jnod = i0 + P_order * i1;

        // Initialise flag: So far, this node hasn't been built
        // or copied yet
        bool node_done = false;

        // Get the pointer to the node in this element (or rather, its clone),
        // returns NULL if there is not node
        Node* created_node_pt = clone_el_pt->get_node_at_local_coordinate(s);

        // Does this node already exist in this element?
        //----------------------------------------------
        if (created_node_pt != 0)
        {
          // Copy node across
          this->node_pt(jnod) = created_node_pt;

          // Make sure that we update the values at the node so that
          // they are consistent with the present representation.
          // This is only need for mixed interpolation where the value
          // at the father could now become active.

          // Loop over all history values
          for (unsigned t = 0; t < ntstorage; t++)
          {
            // Get values from father element
            // Note: get_interpolated_values() sets Vector size itself.
            Vector<double> prev_values;
            clone_el_pt->get_interpolated_values(t, s, prev_values);
            // Find the minimum number of values
            //(either those stored at the node, or those returned by
            // the function)
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

          // Node has been created by copy
          node_done = true;
        }
        // Node does not exist in this element but might already
        //------------------------------------------------------
        // have been created by neighbouring elements
        //-------------------------------------------
        else
        {
          // Was the node created by one of its neighbours
          // Whether or not the node lies on an edge can be calculated
          // by from the fractional position
          bool is_periodic = false;
          created_node_pt = node_created_by_neighbour(s_fraction, is_periodic);

          // If the node was so created, assign the pointers
          if (created_node_pt != 0)
          {
            // If the node is periodic
            if (is_periodic)
            {
              // Now the node must be on a boundary, but we don't know which
              // one
              // The returned created_node_pt is actually the neighbouring
              // periodic node
              Node* neighbour_node_pt = created_node_pt;

              // Determine the edge on which the new node will live
              //(cannot be a vertex node)
              int my_bound = Tree::OMEGA;
              if (s_fraction[0] == 0.0) my_bound = QuadTreeNames::W;
              else if (s_fraction[0] == 1.0)
                my_bound = QuadTreeNames::E;
              else if (s_fraction[1] == 0.0)
                my_bound = QuadTreeNames::S;
              else if (s_fraction[1] == 1.0)
                my_bound = QuadTreeNames::N;

              // Storage for the set of Mesh boundaries on which the
              // appropriate edge lives.
              // [New nodes should always be mid-edge nodes and therefore
              // only live on one boundary but just to play it safe...]
              std::set<unsigned> boundaries;
              // Only get the boundaries if we are at the edge of
              // an element. Nodes in the centre of an element cannot be
              // on Mesh boundaries
              if (my_bound != Tree::OMEGA)
              {
                clone_el_pt->get_boundaries(my_bound, boundaries);
              }

#ifdef PARANOID
              // Case where a new node lives on more than one boundary
              // seems fishy enough to flag
              if (boundaries.size() > 1)
              {
                throw OomphLibError(
                  "boundaries.size()!=1 seems a bit strange..\n",
                  OOMPH_CURRENT_FUNCTION,
                  OOMPH_EXCEPTION_LOCATION);
              }

              // Case when there are no boundaries, we are in big trouble
              if (boundaries.size() == 0)
              {
                std::ostringstream error_stream;
                error_stream << "Periodic node is not on a boundary...\n"
                             << "Coordinates: " << created_node_pt->x(0) << " "
                             << created_node_pt->x(1) << "\n";
                throw OomphLibError(error_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif

              // Create node and set the pointer to it from the element
              created_node_pt =
                this->construct_boundary_node(jnod, time_stepper_pt);
              // Make the node periodic from the neighbour
              created_node_pt->make_periodic(neighbour_node_pt);

              // Loop over # of history values
              for (unsigned t = 0; t < ntstorage; t++)
              {
                // Get position from father element -- this uses the macro
                // element representation if appropriate. If the node
                // turns out to be a hanging node later on, then
                // its position gets adjusted in line with its
                // hanging node interpolation.
                Vector<double> x_prev(2);
                clone_el_pt->get_x(t, s, x_prev);
                // Set previous positions of the new node
                for (unsigned i = 0; i < 2; i++)
                {
                  created_node_pt->x(t, i) = x_prev[i];
                }
              }

              // Check if we need to add nodes to the mesh
              if (mesh_pt != 0)
              {
                // Next, we Update the boundary lookup schemes
                // Loop over the boundaries stored in the set
                for (std::set<unsigned>::iterator it = boundaries.begin();
                     it != boundaries.end();
                     ++it)
                {
                  // Add the node to the boundary
                  mesh_pt->add_boundary_node(*it, created_node_pt);

                  // If we have set an intrinsic coordinate on this
                  // mesh boundary then it must also be interpolated on
                  // the new node
                  // Now interpolate the intrinsic boundary coordinate
                  if (mesh_pt->boundary_coordinate_exists(*it) == true)
                  {
                    Vector<double> zeta(1);
                    clone_el_pt->interpolated_zeta_on_edge(
                      *it, my_bound, s, zeta);

                    created_node_pt->set_coordinates_on_boundary(*it, zeta);
                  }
                }

                // Make sure that we add the node to the mesh
                mesh_pt->add_node_pt(created_node_pt);
              }
            } // End of periodic case
            // Otherwise the node is not periodic, so just set the
            // pointer to the neighbours node
            else
            {
              this->node_pt(jnod) = created_node_pt;
            }
            // Node has been created
            node_done = true;
          }
          // Node does not exist in neighbour element but might already
          //-----------------------------------------------------------
          // have been created by a son of a neighbouring element
          //-----------------------------------------------------
          else
          {
            // Was the node created by one of its neighbours' sons
            // Whether or not the node lies on an edge can be calculated
            // by from the fractional position
            bool is_periodic = false;
            created_node_pt =
              node_created_by_son_of_neighbour(s_fraction, is_periodic);

            // If the node was so created, assign the pointers
            if (created_node_pt != 0)
            {
              // If the node is periodic
              if (is_periodic)
              {
                // Now the node must be on a boundary, but we don't know which
                // one
                // The returned created_node_pt is actually the neighbouring
                // periodic node
                Node* neighbour_node_pt = created_node_pt;

                // Determine the edge on which the new node will live
                //(cannot be a vertex node)
                int my_bound = Tree::OMEGA;
                if (s_fraction[0] == 0.0) my_bound = QuadTreeNames::W;
                else if (s_fraction[0] == 1.0)
                  my_bound = QuadTreeNames::E;
                else if (s_fraction[1] == 0.0)
                  my_bound = QuadTreeNames::S;
                else if (s_fraction[1] == 1.0)
                  my_bound = QuadTreeNames::N;

                // Storage for the set of Mesh boundaries on which the
                // appropriate edge lives.
                // [New nodes should always be mid-edge nodes and therefore
                // only live on one boundary but just to play it safe...]
                std::set<unsigned> boundaries;
                // Only get the boundaries if we are at the edge of
                // an element. Nodes in the centre of an element cannot be
                // on Mesh boundaries
                if (my_bound != Tree::OMEGA)
                {
                  clone_el_pt->get_boundaries(my_bound, boundaries);
                }

#ifdef PARANOID
                // Case where a new node lives on more than one boundary
                // seems fishy enough to flag
                if (boundaries.size() > 1)
                {
                  throw OomphLibError(
                    "boundaries.size()!=1 seems a bit strange..\n",
                    OOMPH_CURRENT_FUNCTION,
                    OOMPH_EXCEPTION_LOCATION);
                }

                // Case when there are no boundaries, we are in big trouble
                if (boundaries.size() == 0)
                {
                  std::ostringstream error_stream;
                  error_stream << "Periodic node is not on a boundary...\n"
                               << "Coordinates: " << created_node_pt->x(0)
                               << " " << created_node_pt->x(1) << "\n";
                  throw OomphLibError(error_stream.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
#endif

                // Create node and set the pointer to it from the element
                created_node_pt =
                  this->construct_boundary_node(jnod, time_stepper_pt);
                // Make the node periodic from the neighbour
                created_node_pt->make_periodic(neighbour_node_pt);

                // Loop over # of history values
                for (unsigned t = 0; t < ntstorage; t++)
                {
                  // Get position from father element -- this uses the macro
                  // element representation if appropriate. If the node
                  // turns out to be a hanging node later on, then
                  // its position gets adjusted in line with its
                  // hanging node interpolation.
                  Vector<double> x_prev(2);
                  clone_el_pt->get_x(t, s, x_prev);
                  // Set previous positions of the new node
                  for (unsigned i = 0; i < 2; i++)
                  {
                    created_node_pt->x(t, i) = x_prev[i];
                  }
                }

                // Check if we need to add nodes to the mesh
                if (mesh_pt != 0)
                {
                  // Next, we Update the boundary lookup schemes
                  // Loop over the boundaries stored in the set
                  for (std::set<unsigned>::iterator it = boundaries.begin();
                       it != boundaries.end();
                       ++it)
                  {
                    // Add the node to the boundary
                    mesh_pt->add_boundary_node(*it, created_node_pt);

                    // If we have set an intrinsic coordinate on this
                    // mesh boundary then it must also be interpolated on
                    // the new node
                    // Now interpolate the intrinsic boundary coordinate
                    if (mesh_pt->boundary_coordinate_exists(*it) == true)
                    {
                      Vector<double> zeta(1);
                      clone_el_pt->interpolated_zeta_on_edge(
                        *it, my_bound, s, zeta);

                      created_node_pt->set_coordinates_on_boundary(*it, zeta);
                    }
                  }

                  // Make sure that we add the node to the mesh
                  mesh_pt->add_node_pt(created_node_pt);
                }
              } // End of periodic case
              // Otherwise the node is not periodic, so just set the
              // pointer to the neighbours node
              else
              {
                this->node_pt(jnod) = created_node_pt;
              }
              // Node has been created
              node_done = true;
            } // Node does not exist in son of neighbouring element
          } // Node does not exist in neighbouring element
        } // Node does not exist in this element

        // Node has not been built anywhere ---> build it here
        if (!node_done)
        {
          // Firstly, we need to determine whether or not a node lies
          // on the boundary before building it, because
          // we actually assign a different type of node on boundaries.

          // Determine the edge on which the new node will live
          //(cannot be a vertex node)
          int my_bound = Tree::OMEGA;
          if (s_fraction[0] == 0.0) my_bound = QuadTreeNames::W;
          else if (s_fraction[0] == 1.0)
            my_bound = QuadTreeNames::E;
          else if (s_fraction[1] == 0.0)
            my_bound = QuadTreeNames::S;
          else if (s_fraction[1] == 1.0)
            my_bound = QuadTreeNames::N;

          // Storage for the set of Mesh boundaries on which the
          // appropriate edge lives.
          // [New nodes should always be mid-edge nodes and therefore
          // only live on one boundary but just to play it safe...]
          std::set<unsigned> boundaries;
          // Only get the boundaries if we are at the edge of
          // an element. Nodes in the centre of an element cannot be
          // on Mesh boundaries
          if (my_bound != Tree::OMEGA)
          {
            clone_el_pt->get_boundaries(my_bound, boundaries);
          }

#ifdef PARANOID
          // Case where a new node lives on more than one boundary
          // seems fishy enough to flag
          if (boundaries.size() > 1)
          {
            throw OomphLibError("boundaries.size()!=1 seems a bit strange..\n",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif

          // If the node lives on a mesh boundary,
          // then we need to create a boundary node
          if (boundaries.size() > 0)
          {
            // Create node and set the pointer to it from the element
            created_node_pt =
              this->construct_boundary_node(jnod, time_stepper_pt);

            // Now we need to work out whether to pin the values at
            // the new node based on the boundary conditions applied at
            // its Mesh boundary

            // Get the boundary conditions from the father
            Vector<int> bound_cons(this->ncont_interpolated_values());
            clone_el_pt->get_bcs(my_bound, bound_cons);

            // Loop over the values and pin, if necessary
            unsigned n_value = created_node_pt->nvalue();
            for (unsigned k = 0; k < n_value; k++)
            {
              if (bound_cons[k])
              {
                created_node_pt->pin(k);
              }
            }

            // Solid node? If so, deal with the positional boundary
            // conditions:
            SolidNode* solid_node_pt =
              dynamic_cast<SolidNode*>(created_node_pt);
            if (solid_node_pt != 0)
            {
              // Get the positional boundary conditions from the father:
              unsigned n_dim = created_node_pt->ndim();
              Vector<int> solid_bound_cons(n_dim);
              RefineableSolidQElement<2>* clone_solid_el_pt =
                dynamic_cast<RefineableSolidQElement<2>*>(clone_el_pt);
#ifdef PARANOID
              if (clone_solid_el_pt == 0)
              {
                std::string error_message =
                  "We have a SolidNode outside a refineable SolidElement\n";
                error_message +=
                  "during mesh refinement -- this doesn't make sense";

                throw OomphLibError(error_message,
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif
              clone_solid_el_pt->get_solid_bcs(my_bound, solid_bound_cons);

              // Loop over the positions and pin, if necessary
              for (unsigned k = 0; k < n_dim; k++)
              {
                if (solid_bound_cons[k])
                {
                  solid_node_pt->pin_position(k);
                }
              }
            } // End of if solid_node_pt


            // Check if we need to add nodes to the mesh
            if (mesh_pt != 0)
            {
              // Next, we Update the boundary lookup schemes
              // Loop over the boundaries stored in the set
              for (std::set<unsigned>::iterator it = boundaries.begin();
                   it != boundaries.end();
                   ++it)
              {
                // Add the node to the boundary
                mesh_pt->add_boundary_node(*it, created_node_pt);

                // If we have set an intrinsic coordinate on this
                // mesh boundary then it must also be interpolated on
                // the new node
                // Now interpolate the intrinsic boundary coordinate
                if (mesh_pt->boundary_coordinate_exists(*it) == true)
                {
                  Vector<double> zeta(1);
                  clone_el_pt->interpolated_zeta_on_edge(
                    *it, my_bound, s, zeta);

                  created_node_pt->set_coordinates_on_boundary(*it, zeta);
                }
              }
            }
          }
          // Otherwise the node is not on a Mesh boundary and
          // we create a normal "bulk" node
          else
          {
            // Create node and set the pointer to it from the element
            created_node_pt = this->construct_node(jnod, time_stepper_pt);
          }

          // Now we set the position and values at the newly created node

          // In the first instance use macro element or FE representation
          // to create past and present nodal positions.
          // (THIS STEP SHOULD NOT BE SKIPPED FOR ALGEBRAIC
          // ELEMENTS AS NOT ALL OF THEM NECESSARILY IMPLEMENT
          // NONTRIVIAL NODE UPDATE FUNCTIONS. CALLING
          // THE NODE UPDATE FOR SUCH ELEMENTS/NODES WILL LEAVE
          // THEIR NODAL POSITIONS WHERE THEY WERE (THIS IS APPROPRIATE
          // ONCE THEY HAVE BEEN GIVEN POSITIONS) BUT WILL
          // NOT ASSIGN SENSIBLE INITIAL POSITONS!

          // Loop over # of history values
          for (unsigned t = 0; t < ntstorage; t++)
          {
            // Get position from father element -- this uses the macro
            // element representation if appropriate. If the node
            // turns out to be a hanging node later on, then
            // its position gets adjusted in line with its
            // hanging node interpolation.
            Vector<double> x_prev(2);
            clone_el_pt->get_x(t, s, x_prev);

            // Set previous positions of the new node
            for (unsigned i = 0; i < 2; i++)
            {
              created_node_pt->x(t, i) = x_prev[i];
            }
          }

          // Loop over all history values
          for (unsigned t = 0; t < ntstorage; t++)
          {
            // Get values from father element
            // Note: get_interpolated_values() sets Vector size itself.
            Vector<double> prev_values;
            clone_el_pt->get_interpolated_values(t, s, prev_values);
            // Initialise the values at the new node
            unsigned n_value = created_node_pt->nvalue();
            for (unsigned k = 0; k < n_value; k++)
            {
              created_node_pt->set_value(t, k, prev_values[k]);
            }
          }

          // Add new node to mesh (if requested)
          if (mesh_pt != 0)
          {
            mesh_pt->add_node_pt(created_node_pt);
          }

          AlgebraicElementBase* alg_el_pt =
            dynamic_cast<AlgebraicElementBase*>(this);

          // If we do have an algebraic element
          if (alg_el_pt != 0)
          {
            std::string error_message = "Have not implemented p-refinement for";
            error_message += "Algebraic p-refineable elements yet\n";

            throw OomphLibError(
              error_message,
              "PRefineableQElement<2,INITIAL_NNODE_1D>::p_refine()",
              OOMPH_EXCEPTION_LOCATION);
          }

        } // End of case when we build the node ourselves

        // Check if the element is an algebraic element
        AlgebraicElementBase* alg_el_pt =
          dynamic_cast<AlgebraicElementBase*>(this);

        // If the element is an algebraic element, setup
        // node position (past and present) from algebraic node update
        // function. This over-writes previous assingments that
        // were made based on the macro-element/FE representation.
        // NOTE: YES, THIS NEEDS TO BE CALLED REPEATEDLY IF THE
        // NODE IS MEMBER OF MULTIPLE ELEMENTS: THEY ALL ASSIGN
        // THE SAME NODAL POSITIONS BUT WE NEED TO ADD THE REMESH
        // INFO FOR *ALL* ROOT ELEMENTS!
        if (alg_el_pt != 0)
        {
          // Build algebraic node update info for new node
          // This sets up the node update data for all node update
          // functions that are shared by all nodes in the father
          // element
          alg_el_pt->setup_algebraic_node_update(
            this->node_pt(jnod), s, clone_el_pt);
        }

      } // End of vertical loop over nodes in element

    } // End of horizontal loop over nodes in element


    // If the element is a MacroElementNodeUpdateElement, set
    // the update parameters for the current element's nodes --
    // all this needs is the vector of (pointers to the)
    // geometric objects that affect the MacroElement-based
    // node update -- this needs to be done to set the node
    // update info for newly created nodes
    MacroElementNodeUpdateElementBase* clone_m_el_pt =
      dynamic_cast<MacroElementNodeUpdateElementBase*>(clone_el_pt);
    if (clone_m_el_pt != 0)
    {
      // Get vector of geometric objects from father (construct vector
      // via copy operation)
      Vector<GeomObject*> geom_object_pt(clone_m_el_pt->geom_object_pt());

      // Cast current element to MacroElementNodeUpdateElement:
      MacroElementNodeUpdateElementBase* m_el_pt =
        dynamic_cast<MacroElementNodeUpdateElementBase*>(this);

#ifdef PARANOID
      if (m_el_pt == 0)
      {
        std::string error_message =
          "Failed to cast to MacroElementNodeUpdateElementBase*\n";
        error_message +=
          "Strange -- if my clone is a MacroElementNodeUpdateElement\n";
        error_message += "then I should be too....\n";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // Build update info by passing vector of geometric objects:
      // This sets the current element to be the update element
      // for all of the element's nodes -- this is reversed
      // if the element is ever un-refined in the father element's
      // rebuild_from_sons() function which overwrites this
      // assignment to avoid nasty segmentation faults that occur
      // when a node tries to update itself via an element that no
      // longer exists...
      m_el_pt->set_node_update_info(geom_object_pt);

      /// / Now loop over nodes in element
      // unsigned n_node = this->nnode();
      // for (unsigned j=0;j<n_node;j++)
      // {
      //  // Get local coordinate in element (Vector sets its own size)
      //  Vector<double> s_in_node_update_element;
      //  this->local_coordinate_of_node(j,s_in_node_update_element);
      //
      //  // Pass the lot to the node
      //  static_cast<MacroElementNodeUpdateNode*>(this->node_pt(j))->
      //   set_node_update_info(this,s_in_node_update_element,m_el_pt->geom_object_pt());
      // }

      /// /BENFLAG:
      // std::cout << "Checking that all the nodes have this as their update
      // element..." << std::endl;
      /// /std::cout << "this = " << this << std::endl;
      // for(unsigned j=0; j<this->nnode(); j++)
      // {
      //  //std::cout << this->node_pt(j) << ":   [" << this->node_pt(j)->x(0)
      //  << "," << this->node_pt(j)->x(1) << "]   update element: " <<
      //  dynamic_cast<MacroElementNodeUpdateNode*>(this->node_pt(j))->node_update_element_pt()
      //  << std::endl; MacroElementNodeUpdateNode* mac_nod_pt =
      //  dynamic_cast<MacroElementNodeUpdateNode*>(this->node_pt(j));
      //  if(mac_nod_pt->node_update_element_pt()!=this)
      //   {
      //    std::cout << "Something's not right! Update element is wrong..." <<
      //    std::endl;
      //   }
      //  FiniteElement* up_el_pt =
      //  dynamic_cast<FiniteElement*>(mac_nod_pt->node_update_element_pt());
      //  bool not_good = true;
      //  for(unsigned l=0; l<up_el_pt->nnode(); l++)
      //   {
      //    if(up_el_pt->node_pt(l)==mac_nod_pt)
      //     {
      //      not_good = false;
      //      break;
      //     }
      //   }
      //  if(not_good==true)
      //   {
      //    std::cout << "Macro node doesn't belong to its update element!" <<
      //    std::endl;
      //   }
      // }


      // Loop over all nodes in element again, to re-set the positions
      // This must be done using the new element's macro-element
      // representation, rather than the old version which may be
      // of a different p-order!
      for (unsigned i0 = 0; i0 < P_order; i0++)
      {
        // Get the fractional position of the node in the direction of s[0]
        s_fraction[0] = this->local_one_d_fraction_of_node(i0, 0);
        // Local coordinate
        s[0] = s_lo[0] + (s_hi[0] - s_lo[0]) * s_fraction[0];

        for (unsigned i1 = 0; i1 < P_order; i1++)
        {
          // Get the fractional position of the node in the direction of s[1]
          s_fraction[1] = this->local_one_d_fraction_of_node(i1, 1);
          // Local coordinate
          s[1] = s_lo[1] + (s_hi[1] - s_lo[1]) * s_fraction[1];

          // Local node number
          jnod = i0 + P_order * i1;

          // Loop over # of history values
          for (unsigned t = 0; t < ntstorage; t++)
          {
            // Get position from father element -- this uses the macro
            // element representation if appropriate. If the node
            // turns out to be a hanging node later on, then
            // its position gets adjusted in line with its
            // hanging node interpolation.
            Vector<double> x_prev(2);
            this->get_x(t, s, x_prev);

            // Set previous positions of the new node
            for (unsigned i = 0; i < 2; i++)
            {
              this->node_pt(jnod)->x(t, i) = x_prev[i];
            }
          }
        }
      }
    }

    // Not necessary to delete the old nodes since all original nodes are in the
    // current mesh and so will be pruned as part of the mesh adaption process.

    // Do any further-build required
    this->further_build();
  }

  //=======================================================================
  /// Shape functions for PRefineableQElement<DIM>
  //=======================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<2, INITIAL_NNODE_1D>::shape(const Vector<double>& s,
                                                       Shape& psi) const
  {
    switch (p_order())
    {
      case 2:
      {
        // Call the OneDimensional Shape functions
        OneDimensionalLegendreShape<2>::calculate_nodal_positions();
        OneDimensionalLegendreShape<2> psi1(s[0]), psi2(s[1]);

        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < 2; i++)
        {
          for (unsigned j = 0; j < 2; j++)
          {
            psi(2 * i + j) = psi2[i] * psi1[j];
          }
        }
        break;
      }
      case 3:
      {
        // Call the OneDimensional Shape functions
        OneDimensionalLegendreShape<3>::calculate_nodal_positions();
        OneDimensionalLegendreShape<3> psi1(s[0]), psi2(s[1]);

        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < 3; i++)
        {
          for (unsigned j = 0; j < 3; j++)
          {
            psi(3 * i + j) = psi2[i] * psi1[j];
          }
        }
        break;
      }
      case 4:
      {
        // Call the OneDimensional Shape functions
        OneDimensionalLegendreShape<4>::calculate_nodal_positions();
        OneDimensionalLegendreShape<4> psi1(s[0]), psi2(s[1]);

        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < 4; i++)
        {
          for (unsigned j = 0; j < 4; j++)
          {
            psi(4 * i + j) = psi2[i] * psi1[j];
          }
        }
        break;
      }
      case 5:
      {
        // Call the OneDimensional Shape functions
        OneDimensionalLegendreShape<5>::calculate_nodal_positions();
        OneDimensionalLegendreShape<5> psi1(s[0]), psi2(s[1]);

        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < 5; i++)
        {
          for (unsigned j = 0; j < 5; j++)
          {
            psi(5 * i + j) = psi2[i] * psi1[j];
          }
        }
        break;
      }
      case 6:
      {
        // Call the OneDimensional Shape functions
        OneDimensionalLegendreShape<6>::calculate_nodal_positions();
        OneDimensionalLegendreShape<6> psi1(s[0]), psi2(s[1]);

        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < 6; i++)
        {
          for (unsigned j = 0; j < 6; j++)
          {
            psi(6 * i + j) = psi2[i] * psi1[j];
          }
        }
        break;
      }
      case 7:
      {
        // Call the OneDimensional Shape functions
        OneDimensionalLegendreShape<7>::calculate_nodal_positions();
        OneDimensionalLegendreShape<7> psi1(s[0]), psi2(s[1]);

        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < 7; i++)
        {
          for (unsigned j = 0; j < 7; j++)
          {
            psi(7 * i + j) = psi2[i] * psi1[j];
          }
        }
        break;
      }
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       polynomial order for shape functions.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
  }

  //=======================================================================
  /// Derivatives of shape functions for PRefineableQElement<DIM>
  //=======================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<2, INITIAL_NNODE_1D>::dshape_local(
    const Vector<double>& s, Shape& psi, DShape& dpsids) const
  {
    switch (p_order())
    {
      case 2:
      {
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<2>::calculate_nodal_positions();
        OneDimensionalLegendreShape<2> psi1(s[0]), psi2(s[1]);
        OneDimensionalLegendreDShape<2> dpsi1ds(s[0]), dpsi2ds(s[1]);

        // Index for the shape functions
        unsigned index = 0;
        // Loop over shape functions in element
        for (unsigned i = 0; i < 2; i++)
        {
          for (unsigned j = 0; j < 2; j++)
          {
            // Assign the values
            dpsids(index, 0) = psi2[i] * dpsi1ds[j];
            dpsids(index, 1) = dpsi2ds[i] * psi1[j];
            psi[index] = psi2[i] * psi1[j];
            // Increment the index
            ++index;
          }
        }
        break;
      }
      case 3:
      {
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<3>::calculate_nodal_positions();
        OneDimensionalLegendreShape<3> psi1(s[0]), psi2(s[1]);
        OneDimensionalLegendreDShape<3> dpsi1ds(s[0]), dpsi2ds(s[1]);

        // Index for the shape functions
        unsigned index = 0;
        // Loop over shape functions in element
        for (unsigned i = 0; i < 3; i++)
        {
          for (unsigned j = 0; j < 3; j++)
          {
            // Assign the values
            dpsids(index, 0) = psi2[i] * dpsi1ds[j];
            dpsids(index, 1) = dpsi2ds[i] * psi1[j];
            psi[index] = psi2[i] * psi1[j];
            // Increment the index
            ++index;
          }
        }
        break;
      }
      case 4:
      {
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<4>::calculate_nodal_positions();
        OneDimensionalLegendreShape<4> psi1(s[0]), psi2(s[1]);
        OneDimensionalLegendreDShape<4> dpsi1ds(s[0]), dpsi2ds(s[1]);

        // Index for the shape functions
        unsigned index = 0;
        // Loop over shape functions in element
        for (unsigned i = 0; i < 4; i++)
        {
          for (unsigned j = 0; j < 4; j++)
          {
            // Assign the values
            dpsids(index, 0) = psi2[i] * dpsi1ds[j];
            dpsids(index, 1) = dpsi2ds[i] * psi1[j];
            psi[index] = psi2[i] * psi1[j];
            // Increment the index
            ++index;
          }
        }
        break;
      }
      case 5:
      {
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<5>::calculate_nodal_positions();
        OneDimensionalLegendreShape<5> psi1(s[0]), psi2(s[1]);
        OneDimensionalLegendreDShape<5> dpsi1ds(s[0]), dpsi2ds(s[1]);

        // Index for the shape functions
        unsigned index = 0;
        // Loop over shape functions in element
        for (unsigned i = 0; i < 5; i++)
        {
          for (unsigned j = 0; j < 5; j++)
          {
            // Assign the values
            dpsids(index, 0) = psi2[i] * dpsi1ds[j];
            dpsids(index, 1) = dpsi2ds[i] * psi1[j];
            psi[index] = psi2[i] * psi1[j];
            // Increment the index
            ++index;
          }
        }
        break;
      }
      case 6:
      {
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<6>::calculate_nodal_positions();
        OneDimensionalLegendreShape<6> psi1(s[0]), psi2(s[1]);
        OneDimensionalLegendreDShape<6> dpsi1ds(s[0]), dpsi2ds(s[1]);

        // Index for the shape functions
        unsigned index = 0;
        // Loop over shape functions in element
        for (unsigned i = 0; i < 6; i++)
        {
          for (unsigned j = 0; j < 6; j++)
          {
            // Assign the values
            dpsids(index, 0) = psi2[i] * dpsi1ds[j];
            dpsids(index, 1) = dpsi2ds[i] * psi1[j];
            psi[index] = psi2[i] * psi1[j];
            // Increment the index
            ++index;
          }
        }
        break;
      }
      case 7:
      {
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<7>::calculate_nodal_positions();
        OneDimensionalLegendreShape<7> psi1(s[0]), psi2(s[1]);
        OneDimensionalLegendreDShape<7> dpsi1ds(s[0]), dpsi2ds(s[1]);

        // Index for the shape functions
        unsigned index = 0;
        // Loop over shape functions in element
        for (unsigned i = 0; i < 7; i++)
        {
          for (unsigned j = 0; j < 7; j++)
          {
            // Assign the values
            dpsids(index, 0) = psi2[i] * dpsi1ds[j];
            dpsids(index, 1) = dpsi2ds[i] * psi1[j];
            psi[index] = psi2[i] * psi1[j];
            // Increment the index
            ++index;
          }
        }
        break;
      }
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       polynomial order for shape functions.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
  }

  //=======================================================================
  /// Second derivatives of shape functions for PRefineableQElement<DIM>
  ///  d2psids(i,0) = \f$ d^2 \psi_j / d s^2 \f$
  //=======================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<2, INITIAL_NNODE_1D>::d2shape_local(
    const Vector<double>& s, Shape& psi, DShape& dpsids, DShape& d2psids) const
  {
    std::ostringstream error_message;
    error_message
      << "\nd2shape_local currently not implemented for this element\n";
    throw OomphLibError(
      error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  //=======================================================================
  /// Rebuild the element from nodes found in its sons
  /// Adjusts its p-order to be the maximum of its sons' p-orders
  //=======================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<2, INITIAL_NNODE_1D>::rebuild_from_sons(
    Mesh*& mesh_pt)
  {
    using namespace QuadTreeNames;

    // Get p-orders of sons
    unsigned n_sons = this->tree_pt()->nsons();
    Vector<unsigned> son_p_order(n_sons);
    unsigned max_son_p_order = 0;
    for (unsigned ison = 0; ison < n_sons; ison++)
    {
      PRefineableElement* el_pt = dynamic_cast<PRefineableElement*>(
        this->tree_pt()->son_pt(ison)->object_pt());
      son_p_order[ison] = el_pt->p_order();
      if (son_p_order[ison] > max_son_p_order)
        max_son_p_order = son_p_order[ison];
    }

    unsigned old_Nnode = this->nnode();
    unsigned old_P_order = this->p_order();
    // Set p-order of the element
    this->p_order() = max_son_p_order;

    // Change integration scheme
    delete this->integral_pt();
    switch (this->p_order())
    {
      case 2:
        this->set_integration_scheme(new GaussLobattoLegendre<2, 2>);
        break;
      case 3:
        this->set_integration_scheme(new GaussLobattoLegendre<2, 3>);
        break;
      case 4:
        this->set_integration_scheme(new GaussLobattoLegendre<2, 4>);
        break;
      case 5:
        this->set_integration_scheme(new GaussLobattoLegendre<2, 5>);
        break;
      case 6:
        this->set_integration_scheme(new GaussLobattoLegendre<2, 6>);
        break;
      case 7:
        this->set_integration_scheme(new GaussLobattoLegendre<2, 7>);
        break;
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       integration scheme.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }

    // Back up pointers to old nodes before they are lost
    Vector<Node*> old_node_pt(old_Nnode);
    for (unsigned n = 0; n < old_Nnode; n++)
    {
      old_node_pt[n] = this->node_pt(n);
    }

    // Allocate new space for Nodes (at the element level)
    this->set_n_node(this->p_order() * this->p_order());

    // Copy vertex nodes which were populated in the pre-build
    this->node_pt(0) = old_node_pt[0];
    this->node_pt(this->p_order() - 1) = old_node_pt[old_P_order - 1];
    this->node_pt(this->p_order() * (this->p_order() - 1)) =
      old_node_pt[(old_P_order) * (old_P_order - 1)];
    this->node_pt(this->p_order() * this->p_order() - 1) =
      old_node_pt[(old_P_order) * (old_P_order)-1];

    // Copy midpoint nodes from sons if new p-order is odd
    if (this->p_order() % 2 == 1)
    {
      // Work out which is midpoint node
      unsigned midpoint = (this->p_order() - 1) / 2;

      // Bottom edge
      this->node_pt(midpoint) = dynamic_cast<RefineableQElement<2>*>(
                                  quadtree_pt()->son_pt(SW)->object_pt())
                                  ->vertex_node_pt(1);
      // Left edge
      this->node_pt(midpoint * this->p_order()) =
        dynamic_cast<RefineableQElement<2>*>(
          quadtree_pt()->son_pt(SW)->object_pt())
          ->vertex_node_pt(2);
      // Top edge
      this->node_pt((this->p_order() - 1) * this->p_order() + midpoint) =
        dynamic_cast<RefineableQElement<2>*>(
          quadtree_pt()->son_pt(NE)->object_pt())
          ->vertex_node_pt(2);
      // Right edge
      this->node_pt((midpoint + 1) * this->p_order() - 1) =
        dynamic_cast<RefineableQElement<2>*>(
          quadtree_pt()->son_pt(NE)->object_pt())
          ->vertex_node_pt(1);
    }


    // The timestepper should be the same for all nodes and node 0 should
    // never be deleted.
    if (this->node_pt(0) == 0)
    {
      throw OomphLibError("The Corner node (0) does not exist",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    TimeStepper* time_stepper_pt = this->node_pt(0)->time_stepper_pt();
    unsigned ntstorage = time_stepper_pt->ntstorage();

    unsigned jnod = 0;
    Vector<double> s_fraction(2), s(2);
    // Loop over the nodes in the element
    unsigned n_p = this->nnode_1d();
    for (unsigned i0 = 0; i0 < n_p; i0++)
    {
      // Get the fractional position of the node
      s_fraction[0] = this->local_one_d_fraction_of_node(i0, 0);
      // Local coordinate
      s[0] = -1.0 + 2.0 * s_fraction[0];

      for (unsigned i1 = 0; i1 < n_p; i1++)
      {
        // Get the fractional position of the node in the direction of s[1]
        s_fraction[1] = this->local_one_d_fraction_of_node(i1, 1);
        // Local coordinate in father element
        s[1] = -1.0 + 2.0 * s_fraction[1];

        // Set the local node number
        jnod = i0 + n_p * i1;

        // Initialise flag: So far, this node hasn't been built
        // or copied yet
        bool node_done = false;

        // Get the pointer to the node in this element, returns NULL
        // if there is not node
        Node* created_node_pt = this->get_node_at_local_coordinate(s);

        // Does this node already exist in this element?
        //----------------------------------------------
        if (created_node_pt != 0)
        {
          // Copy node across
          this->node_pt(jnod) = created_node_pt;

          // Node has been created by copy
          node_done = true;
        }
        // Node does not exist in this element but might already
        //------------------------------------------------------
        // have been created by neighbouring elements
        //-------------------------------------------
        else
        {
          // Was the node created by one of its neighbours
          // Whether or not the node lies on an edge can be calculated
          // by from the fractional position
          bool is_periodic = false;
          created_node_pt = node_created_by_neighbour(s_fraction, is_periodic);

          // If the node was so created, assign the pointers
          if (created_node_pt != 0)
          {
            // If the node is periodic
            if (is_periodic)
            {
              throw OomphLibError("Cannot handle periodic nodes yet",
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
            // Non-periodic case, just set the pointer
            else
            {
              this->node_pt(jnod) = created_node_pt;
            }
            // Node has been created
            node_done = true;
          }
        } // Node does not exist in this element

        // Node has not been built anywhere ---> build it here
        if (!node_done)
        {
          // First, find the son element in which the node should live

          // Find coordinates in the sons
          Vector<double> s_in_son(2);
          using namespace QuadTreeNames;
          int son = -10;
          // If negative on the west side
          if (s_fraction[0] < 0.5)
          {
            // On the south side
            if (s_fraction[1] < 0.5)
            {
              // It's the southwest son
              son = SW;
              s_in_son[0] = -1.0 + 4.0 * s_fraction[0];
              s_in_son[1] = -1.0 + 4.0 * s_fraction[1];
            }
            // On the north side
            else
            {
              // It's the northwest son
              son = NW;
              s_in_son[0] = -1.0 + 4.0 * s_fraction[0];
              s_in_son[1] = -1.0 + 4.0 * (s_fraction[1] - 0.5);
            }
          }
          else
          {
            // On the south side
            if (s_fraction[1] < 0.5)
            {
              // It's the southeast son
              son = SE;
              s_in_son[0] = -1.0 + 4.0 * (s_fraction[0] - 0.5);
              s_in_son[1] = -1.0 + 4.0 * s_fraction[1];
            }
            // On the north side
            else
            {
              // It's the northeast son
              son = NE;
              s_in_son[0] = -1.0 + 4.0 * (s_fraction[0] - 0.5);
              s_in_son[1] = -1.0 + 4.0 * (s_fraction[1] - 0.5);
            }
          }

          // Get the pointer to the son element in which the new node
          // would live
          PRefineableQElement<2, INITIAL_NNODE_1D>* son_el_pt =
            dynamic_cast<PRefineableQElement<2, INITIAL_NNODE_1D>*>(
              this->tree_pt()->son_pt(son)->object_pt());

          // If we are rebuilding, then worry about the boundary conditions
          // Find the boundary of the node
          // Initially none
          int boundary = Tree::OMEGA;
          // If we are on the western boundary
          if (i0 == 0)
          {
            boundary = W;
          }
          // If we are on the eastern boundary
          else if (i0 == n_p - 1)
          {
            boundary = E;
          }

          // If we are on the southern boundary
          if (i1 == 0)
          {
            // If we already have already set the boundary, we're on a corner
            switch (boundary)
            {
              case W:
                boundary = SW;
                break;
              case E:
                boundary = SE;
                break;
              // Boundary not set
              default:
                boundary = S;
                break;
            }
          }
          // If we are the northern bounadry
          else if (i1 == n_p - 1)
          {
            // If we already have a boundary
            switch (boundary)
            {
              case W:
                boundary = NW;
                break;
              case E:
                boundary = NE;
                break;
              default:
                boundary = N;
                break;
            }
          }

          // set of boundaries that this edge in the son lives on
          std::set<unsigned> boundaries;

          // Now get the boundary conditions from the son
          // The boundaries will be common to the son because there can be
          // no rotations here
          if (boundary != Tree::OMEGA)
          {
            son_el_pt->get_boundaries(boundary, boundaries);
          }

          // If the node lives on a boundary:
          // Construct a boundary node,
          // Get boundary conditions and
          // update all lookup schemes
          if (boundaries.size() > 0)
          {
            // Construct the new node
            created_node_pt =
              this->construct_boundary_node(jnod, time_stepper_pt);

            // Get the boundary conditions from the son
            Vector<int> bound_cons(this->ncont_interpolated_values());
            son_el_pt->get_bcs(boundary, bound_cons);

            // Loop over the values and pin if necessary
            unsigned nval = created_node_pt->nvalue();
            for (unsigned k = 0; k < nval; k++)
            {
              if (bound_cons[k])
              {
                created_node_pt->pin(k);
              }
            }

            // Solid node? If so, deal with the positional boundary
            // conditions:
            SolidNode* solid_node_pt =
              dynamic_cast<SolidNode*>(created_node_pt);
            if (solid_node_pt != 0)
            {
              // Get the positional boundary conditions from the father:
              unsigned n_dim = created_node_pt->ndim();
              Vector<int> solid_bound_cons(n_dim);
              RefineableSolidQElement<2>* son_solid_el_pt =
                dynamic_cast<RefineableSolidQElement<2>*>(son_el_pt);
#ifdef PARANOID
              if (son_solid_el_pt == 0)
              {
                std::string error_message =
                  "We have a SolidNode outside a refineable SolidElement\n";
                error_message +=
                  "during mesh refinement -- this doesn't make sense\n";

                throw OomphLibError(error_message,
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif
              son_solid_el_pt->get_solid_bcs(boundary, solid_bound_cons);

              // Loop over the positions and pin, if necessary
              for (unsigned k = 0; k < n_dim; k++)
              {
                if (solid_bound_cons[k])
                {
                  solid_node_pt->pin_position(k);
                }
              }
            } // End of if solid_node_pt


            // Next we update the boundary look-up schemes
            // Loop over the boundaries stored in the set
            for (std::set<unsigned>::iterator it = boundaries.begin();
                 it != boundaries.end();
                 ++it)
            {
              // Add the node to the boundary
              mesh_pt->add_boundary_node(*it, created_node_pt);

              // If we have set an intrinsic coordinate on this
              // mesh boundary then it must also be interpolated on
              // the new node
              // Now interpolate the intrinsic boundary coordinate
              if (mesh_pt->boundary_coordinate_exists(*it) == true)
              {
                Vector<double> zeta(1);
                son_el_pt->interpolated_zeta_on_edge(
                  *it, boundary, s_in_son, zeta);

                created_node_pt->set_coordinates_on_boundary(*it, zeta);
              }
            }
          }
          // Otherwise the node is not on a Mesh boundary
          // and we create a normal "bulk" node
          else
          {
            // Construct the new node
            created_node_pt = this->construct_node(jnod, time_stepper_pt);
          }

          // Now we set the position and values at the newly created node

          // In the first instance use macro element or FE representation
          // to create past and present nodal positions.
          // (THIS STEP SHOULD NOT BE SKIPPED FOR ALGEBRAIC
          // ELEMENTS AS NOT ALL OF THEM NECESSARILY IMPLEMENT
          // NONTRIVIAL NODE UPDATE FUNCTIONS. CALLING
          // THE NODE UPDATE FOR SUCH ELEMENTS/NODES WILL LEAVE
          // THEIR NODAL POSITIONS WHERE THEY WERE (THIS IS APPROPRIATE
          // ONCE THEY HAVE BEEN GIVEN POSITIONS) BUT WILL
          // NOT ASSIGN SENSIBLE INITIAL POSITONS!

          // Loop over # of history values
          // Loop over # of history values
          for (unsigned t = 0; t < ntstorage; t++)
          {
            using namespace QuadTreeNames;
            // Get the position from the son
            Vector<double> x_prev(2);

            // Now let's fill in the value
            son_el_pt->get_x(t, s_in_son, x_prev);
            for (unsigned i = 0; i < 2; i++)
            {
              created_node_pt->x(t, i) = x_prev[i];
            }
          }

          // Now set up the values
          // Loop over all history values
          for (unsigned t = 0; t < ntstorage; t++)
          {
            // Get values from father element
            // Note: get_interpolated_values() sets Vector size itself.
            Vector<double> prev_values;
            son_el_pt->get_interpolated_values(t, s_in_son, prev_values);

            // Initialise the values at the new node
            for (unsigned k = 0; k < created_node_pt->nvalue(); k++)
            {
              created_node_pt->set_value(t, k, prev_values[k]);
            }
          }

          // Add the node to the mesh
          mesh_pt->add_node_pt(created_node_pt);

          // Check if the element is an algebraic element
          AlgebraicElementBase* alg_el_pt =
            dynamic_cast<AlgebraicElementBase*>(this);

          // If we do have an algebraic element
          if (alg_el_pt != 0)
          {
            std::string error_message =
              "Have not implemented rebuilding from sons for";
            error_message += "Algebraic p-refineable elements yet\n";

            throw OomphLibError(
              error_message,
              "PRefineableQElement<2,INITIAL_NNODE_1D>::rebuild_from_sons()",
              OOMPH_EXCEPTION_LOCATION);
          }

        } // End of the case when we build the node ourselves
      }
    } // End of loop over all nodes in element


    // If the element is a MacroElementNodeUpdateElement, set the update
    // parameters for the current element's nodes. These need to be reset
    // (as in MacroElementNodeUpdateElement<ELEMENT>::rebuild_from_sons())
    // because the nodes in this element have changed
    MacroElementNodeUpdateElementBase* m_el_pt =
      dynamic_cast<MacroElementNodeUpdateElementBase*>(this);
    if (m_el_pt != 0)
    {
      // Loop over the nodes
      for (unsigned j = 0; j < this->nnode(); j++)
      {
        // Get local coordinate in element (Vector sets its own size)
        Vector<double> s_in_node_update_element;
        this->local_coordinate_of_node(j, s_in_node_update_element);

        // Get vector of geometric objects
        Vector<GeomObject*> geom_object_pt(m_el_pt->geom_object_pt());

        // Pass the lot to the node
        static_cast<MacroElementNodeUpdateNode*>(this->node_pt(j))
          ->set_node_update_info(
            this, s_in_node_update_element, geom_object_pt);
      }
    }

    // MacroElementNodeUpdateElementBase* m_el_pt=dynamic_cast<
    // MacroElementNodeUpdateElementBase*>(this);
    // if(m_el_pt!=0)
    // {
    // Loop over all nodes in element again, to re-set the positions
    // This must be done using the new element's macro-element
    // representation, rather than the old version which may be
    // of a different p-order!
    for (unsigned i0 = 0; i0 < n_p; i0++)
    {
      // Get the fractional position of the node
      s_fraction[0] = this->local_one_d_fraction_of_node(i0, 0);
      // Local coordinate
      s[0] = -1.0 + 2.0 * s_fraction[0];

      for (unsigned i1 = 0; i1 < n_p; i1++)
      {
        // Get the fractional position of the node in the direction of s[1]
        s_fraction[1] = this->local_one_d_fraction_of_node(i1, 1);
        // Local coordinate in father element
        s[1] = -1.0 + 2.0 * s_fraction[1];

        // Set the local node number
        jnod = i0 + n_p * i1;

        // Loop over # of history values
        for (unsigned t = 0; t < ntstorage; t++)
        {
          // Get position from father element -- this uses the macro
          // element representation if appropriate. If the node
          // turns out to be a hanging node later on, then
          // its position gets adjusted in line with its
          // hanging node interpolation.
          Vector<double> x_prev(2);
          this->get_x(t, s, x_prev);

          // Set previous positions of the new node
          for (unsigned i = 0; i < 2; i++)
          {
            this->node_pt(jnod)->x(t, i) = x_prev[i];
          }
        }
      }
    }
    // }
  }

  //=================================================================
  /// Check inter-element continuity of
  /// - nodal positions
  /// - (nodally) interpolated function values
  /// Overloaded to not check differences in the value. Mortaring
  /// doesn't enforce strong continuity between elements.
  //====================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<2, INITIAL_NNODE_1D>::check_integrity(
    double& max_error)
  {
    // Overloaded to *not* check for continuity in value of interpolated
    // variables. This is necessary because mortaring does not ensure continuity
    // across element boundaries. It therefore makes no sense to test for this.

    // Dummy set max_error to 0
    max_error = 0.0;

    // With macro-elements, (strong) continuity in position is nolonger
    // guaranteed either, so we don't check for this either. In fact, we do
    // nothing at all.
    if (this->macro_elem_pt() != 0)
    {
      // We have a macro element, so do nothing!
      return;
    }

    using namespace QuadTreeNames;

    // Number of nodes along edge
    unsigned n_p = nnode_1d();

    // Number of timesteps (incl. present) for which continuity is
    // to be checked.
    unsigned n_time = 1;

    // Initialise errors
    max_error = 0.0;
    Vector<double> max_error_x(2, 0.0);
    double max_error_val = 0.0;

    Vector<int> edges(4);
    edges[0] = S;
    edges[1] = N;
    edges[2] = W;
    edges[3] = E;

    // Loop over the edges
    for (unsigned edge_counter = 0; edge_counter < 4; edge_counter++)
    {
      Vector<unsigned> translate_s(2);
      Vector<double> s(2), s_lo_neigh(2), s_hi_neigh(2), s_fraction(2);
      int neigh_edge, diff_level;
      bool in_neighbouring_tree;

      // Find pointer to neighbour in this direction
      QuadTree* neigh_pt;
      neigh_pt = quadtree_pt()->gteq_edge_neighbour(edges[edge_counter],
                                                    translate_s,
                                                    s_lo_neigh,
                                                    s_hi_neigh,
                                                    neigh_edge,
                                                    diff_level,
                                                    in_neighbouring_tree);

      // Neighbour exists and has existing nodes
      if ((neigh_pt != 0) && (neigh_pt->object_pt()->nodes_built()))
      {
        // Need to exclude periodic nodes from this check
        // There are only periodic nodes if we are in a neighbouring tree
        bool is_periodic = false;
        if (in_neighbouring_tree)
        {
          // Is it periodic
          is_periodic = this->tree_pt()->root_pt()->is_neighbour_periodic(
            edges[edge_counter]);
        }

        // We also need to exclude edges which may have hanging nodes
        // because mortaring does not guarantee (strong) continuity
        // in position or in value at nonconforming element boundaries
        bool exclude_this_edge = false;
        if (diff_level != 0)
        {
          // h-type nonconformity (dependent)
          exclude_this_edge = true;
        }
        else if (neigh_pt->nsons() != 0)
        {
          // h-type nonconformity (master)
          exclude_this_edge = true;
        }
        else
        {
          unsigned my_p_order = this->p_order();
          unsigned neigh_p_order =
            dynamic_cast<PRefineableQElement*>(neigh_pt->object_pt())
              ->p_order();
          if (my_p_order != neigh_p_order)
          {
            // p-type nonconformity
            exclude_this_edge = true;
          }
        }

        // With macro-elements, (strong) continuity in position is nolonger
        // guaranteed either, so we don't check for this either. In fact, we do
        // nothing at all.
        if (dynamic_cast<FiniteElement*>(neigh_pt->object_pt())
              ->macro_elem_pt() != 0)
        {
          // We have a macro element, so do nothing!
          break;
        }

        // Check conforming edges
        if (!exclude_this_edge)
        {
          // Loop over nodes along the edge
          for (unsigned i0 = 0; i0 < n_p; i0++)
          {
            // Storage for pointer to the local node
            Node* local_node_pt = 0;

            switch (edge_counter)
            {
              case 0:
                // Local fraction of node
                s_fraction[0] = this->local_one_d_fraction_of_node(i0, 0);
                s_fraction[1] = 0.0;
                // Get pointer to local node
                local_node_pt = this->node_pt(i0);
                break;

              case 1:
                // Local fraction of node
                s_fraction[0] = this->local_one_d_fraction_of_node(i0, 0);
                s_fraction[1] = 1.0;
                // Get pointer to local node
                local_node_pt = this->node_pt(i0 + n_p * (n_p - 1));
                break;

              case 2:
                // Local fraction of node
                s_fraction[0] = 0.0;
                s_fraction[1] = this->local_one_d_fraction_of_node(i0, 1);
                // Get pointer to local node
                local_node_pt = this->node_pt(n_p * i0);
                break;

              case 3:
                // Local fraction of node
                s_fraction[0] = 1.0;
                s_fraction[1] = this->local_one_d_fraction_of_node(i0, 1);
                // Get pointer to local node
                local_node_pt = this->node_pt(n_p - 1 + n_p * i0);
                break;
            }

            // Calculate the local coordinate and the local coordinate as viewed
            // from the neighbour
            Vector<double> s_in_neighb(2);
            for (unsigned i = 0; i < 2; i++)
            {
              // Local coordinate in this element
              s[i] = -1.0 + 2.0 * s_fraction[i];
              // Local coordinate in the neighbour
              s_in_neighb[i] =
                s_lo_neigh[i] +
                s_fraction[translate_s[i]] * (s_hi_neigh[i] - s_lo_neigh[i]);
            }

            // Loop over timesteps
            for (unsigned t = 0; t < n_time; t++)
            {
              // Get the nodal position from neighbour element
              Vector<double> x_in_neighb(2);
              neigh_pt->object_pt()->interpolated_x(
                t, s_in_neighb, x_in_neighb);

              // Check error only if the node is NOT periodic
              if (is_periodic == false)
              {
                for (int i = 0; i < 2; i++)
                {
                  // Find the spatial error
                  double err =
                    std::fabs(local_node_pt->x(t, i) - x_in_neighb[i]);

                  // If it's bigger than our tolerance, say so
                  if (err > 1e-9)
                  {
                    oomph_info << "errx " << err << " " << t << " "
                               << local_node_pt->x(t, i) << " "
                               << x_in_neighb[i] << std::endl;

                    oomph_info << "at " << local_node_pt->x(0) << " "
                               << local_node_pt->x(1) << std::endl;
                  }

                  // If it's bigger than the previous max error, it is the
                  // new max error!
                  if (err > max_error_x[i])
                  {
                    max_error_x[i] = err;
                  }
                }
              }

              // Get the values from neighbour element. Note: # of values
              // gets set by routine (because in general we don't know
              // how many interpolated values a certain element has
              Vector<double> values_in_neighb;
              neigh_pt->object_pt()->get_interpolated_values(
                t, s_in_neighb, values_in_neighb);

              // Get the values in current element.
              Vector<double> values;
              this->get_interpolated_values(t, s, values);

              // Now figure out how many continuously interpolated values there
              // are
              unsigned num_val =
                neigh_pt->object_pt()->ncont_interpolated_values();

              // Check error
              for (unsigned ival = 0; ival < num_val; ival++)
              {
                double err = std::fabs(values[ival] - values_in_neighb[ival]);

                if (err > 1.0e-10)
                {
                  oomph_info << local_node_pt->x(0) << " "
                             << local_node_pt->x(1) << " \n# "
                             << "erru (S)" << err << " " << ival << " "
                             << this->get_node_number(local_node_pt) << " "
                             << values[ival] << " " << values_in_neighb[ival]
                             << std::endl;
                }

                if (err > max_error_val)
                {
                  max_error_val = err;
                }
              }
            }
          }
        }
      }
    }

    max_error = max_error_x[0];
    if (max_error_x[1] > max_error) max_error = max_error_x[1];
    if (max_error_val > max_error) max_error = max_error_val;

    if (max_error > 1e-9)
    {
      oomph_info << "\n#------------------------------------ \n#Max error ";
      oomph_info << max_error_x[0] << " " << max_error_x[1] << " "
                 << max_error_val << std::endl;
      oomph_info << "#------------------------------------ \n " << std::endl;
    }
  }

  //=================================================================
  /// Internal function to set up the hanging nodes on a particular
  /// edge of the element.
  /// Implements the mortarting method to enforce continuity weakly
  /// across non-conforming element boundaries \f$\Gamma\f$ using an
  /// integral matching condition
  /// \f[ \int_\Gamma (u_{\mbox{S}} - u_{\mbox{M}}) \psi \mbox{d} s = 0 \f]
  /// for all polynomials \f$\psi\f$ on \f$\Gamma\f$ of degree at most
  /// p-2 (where p is the spectral-order of the dependent element) and a
  /// vertex matching condition
  /// \f[ (u_{\mbox{S}} - u_{\mbox{M}})\big\vert_{\partial\Gamma} = 0.\f]
  ///
  /// The algorithm works as follows:
  ///  - First the element determines if its edge my_edge is on the
  ///    master or dependent side of the non-conformity.
  /// At h-type non-conformities
  ///    we choose long edges to be masters, and at p-type nonconformities the
  ///    edge with lower p-order is the master.
  ///  - Mortaring is performed by the dependent side.
  ///  - If a vertex node of the mortar is shared between master and dependent
  ///    element then the vertex matching condition is enforced automatically.
  ///    Otherwise it must be imposed by constraining its value to that at on
  ///    the master side.
  ///  - The integral matching condition is discretised and the mortar test
  ///    functions \f$ \psi \f$ are chosen to be derivatives of Legendre
  ///    polynomials of degree p-1.
  ///  - The mortar mass matrix M is constructed. Its entries are the
  ///    mortar test functions evaluated at the dependent nodal positions,
  ///    so it is diagonal.
  ///  - The local projection matrix is constructed for the master element by
  ///    applying the discretised integral matching condition along the mortar
  ///    using the appropriate quadrature order.
  ///  - The global projection matrix is then assembled by subtracting
  ///    contributions from the mortar vertex nodes.
  ///  - The mortar system \f$ M\xi^s = P\hat{\xi^m} \f$ is constructed,
  ///    where \f$ \xi^m \f$ and \f$ \xi^s \f$ are the nodal values at
  ///    the master and dependent nodes respectively.
  ///  - The conformity matrix \f$ C = M^{-1}P \f$ is computed. This is
  ///    straightforward since the mass matrix is diagonal.
  ///  - Finally, the master nodes and weights for each dependent node
  ///    are read from
  ///    the conformity matrix and stored in the dependents' hanging schemes.
  ///
  /// The positions of the dependent nodes are set to be consistent with their
  /// hanging schemes.
  //=================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<2, INITIAL_NNODE_1D>::quad_hang_helper(
    const int& value_id, const int& my_edge, std::ofstream& output_hangfile)
  {
    using namespace QuadTreeNames;

    Vector<unsigned> translate_s(2);
    Vector<double> s_lo_neigh(2);
    Vector<double> s_hi_neigh(2);
    int neigh_edge, diff_level;
    bool in_neighbouring_tree;

    // Find pointer to neighbour in this direction
    QuadTree* neigh_pt;
    neigh_pt = this->quadtree_pt()->gteq_edge_neighbour(my_edge,
                                                        translate_s,
                                                        s_lo_neigh,
                                                        s_hi_neigh,
                                                        neigh_edge,
                                                        diff_level,
                                                        in_neighbouring_tree);

    // Work out master/dependent edges
    //----------------------------

    // Set up booleans
    // bool h_type_master = false;
    bool h_type_dependent = false;
    // bool p_type_master = false;
    bool p_type_dependent = false;

    // Neighbour exists and all nodes have been created
    if (neigh_pt != 0)
    {
      // Check if neighbour is bigger than me
      if (diff_level != 0)
      {
        // Dependent at h-type non-conformity
        h_type_dependent = true;
      }
      // Check if neighbour is the same size as me
      else if (neigh_pt->nsons() == 0)
      {
        // Neighbour is same size as me
        // Find p-orders of each element
        unsigned my_p_order =
          dynamic_cast<PRefineableQElement<2, INITIAL_NNODE_1D>*>(this)
            ->p_order();
        unsigned neigh_p_order =
          dynamic_cast<PRefineableQElement<2, INITIAL_NNODE_1D>*>(
            neigh_pt->object_pt())
            ->p_order();

        // Check for p-type non-conformity
        if (neigh_p_order == my_p_order)
        {
          // At a conforming interface
        }
        else if (neigh_p_order < my_p_order)
        {
          // Dependent at p-type non-conformity
          p_type_dependent = true;
        }
        else
        {
          // Master at p-type non-conformity
          // p_type_master = true;
        }
      }
      // Neighbour must be smaller than me
      else
      {
        // Master at h-type non-conformity
        // h_type_master = true;
      }
    }
    else
    {
      // Edge is on a boundary
    }

    // Now do the mortaring
    //---------------------
    if (h_type_dependent || p_type_dependent)
    {
      // Compute the active coordinate index along the this side of mortar
      unsigned active_coord_index;
      if (my_edge == N || my_edge == S) active_coord_index = 0;
      else if (my_edge == E || my_edge == W)
        active_coord_index = 1;
      else
      {
        throw OomphLibError("Cannot transform coordinates",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Get pointer to neighbouring master element (in p-refineable form)
      PRefineableQElement<2, INITIAL_NNODE_1D>* neigh_obj_pt;
      neigh_obj_pt = dynamic_cast<PRefineableQElement<2, INITIAL_NNODE_1D>*>(
        neigh_pt->object_pt());

      // Create vector of master and dependent nodes
      //----------------------------------------
      Vector<Node*> master_node_pt, dependent_node_pt;
      Vector<unsigned> master_node_number, dependent_node_number;
      Vector<Vector<double>> dependent_node_s_fraction;

      // Number of nodes in one dimension
      const unsigned my_n_p = this->ninterpolating_node_1d(value_id);
      const unsigned neigh_n_p = neigh_obj_pt->ninterpolating_node_1d(value_id);

      // Test for the periodic node case
      // Are we crossing a periodic boundary
      bool is_periodic = false;
      if (in_neighbouring_tree)
      {
        is_periodic =
          this->tree_pt()->root_pt()->is_neighbour_periodic(my_edge);
      }

      // If it is periodic we actually need to get the node in
      // the neighbour of the neighbour (which will be a parent of
      // the present element) so that the "fixed" coordinate is
      // correctly calculated.
      // The idea is to replace the neigh_pt and associated data
      // with those of the neighbour of the neighbour
      if (is_periodic)
      {
        throw OomphLibError(
          "Cannot do mortaring with periodic hanging nodes yet! (Fix me!)",
          "PRefineableQElement<2,INITIAL_NNODE_1D>::quad_hang_helper()",
          OOMPH_EXCEPTION_LOCATION);
      } // End of special treatment for periodic hanging nodes

      // Storage for pointers to the nodes and their numbers along the master
      // edge
      unsigned neighbour_node_number = 0;
      Node* neighbour_node_pt = 0;

      // Loop over nodes along the edge
      for (unsigned i0 = 0; i0 < neigh_n_p; i0++)
      {
        // Find the neighbour's node
        switch (neigh_edge)
        {
          case N:
            neighbour_node_number = i0 + neigh_n_p * (neigh_n_p - 1);
            neighbour_node_pt = neigh_obj_pt->interpolating_node_pt(
              neighbour_node_number, value_id);
            break;

          case S:
            neighbour_node_number = i0;
            neighbour_node_pt = neigh_obj_pt->interpolating_node_pt(
              neighbour_node_number, value_id);
            break;

          case E:
            neighbour_node_number = neigh_n_p - 1 + neigh_n_p * i0;
            neighbour_node_pt = neigh_obj_pt->interpolating_node_pt(
              neighbour_node_number, value_id);
            break;

          case W:
            neighbour_node_number = neigh_n_p * i0;
            neighbour_node_pt = neigh_obj_pt->interpolating_node_pt(
              neighbour_node_number, value_id);
            break;

          default:
            throw OomphLibError("my_edge not N, S, W, E\n",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        // Set node as master node
        master_node_number.push_back(neighbour_node_number);
        master_node_pt.push_back(neighbour_node_pt);
      }

      // Storage for pointers to the local nodes and their numbers along my edge
      unsigned local_node_number = 0;
      Node* local_node_pt = 0;

      // Loop over the nodes along my edge
      for (unsigned i0 = 0; i0 < my_n_p; i0++)
      {
        // Storage for the fractional position of the node in the element
        Vector<double> s_fraction(2);

        // Find the local node and the fractional position of the node
        // which depends on the edge, of course
        switch (my_edge)
        {
          case N:
            s_fraction[0] =
              local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
            s_fraction[1] = 1.0;
            local_node_number = i0 + my_n_p * (my_n_p - 1);
            local_node_pt =
              this->interpolating_node_pt(local_node_number, value_id);
            break;

          case S:
            s_fraction[0] =
              local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
            s_fraction[1] = 0.0;
            local_node_number = i0;
            local_node_pt =
              this->interpolating_node_pt(local_node_number, value_id);
            break;

          case E:
            s_fraction[0] = 1.0;
            s_fraction[1] =
              local_one_d_fraction_of_interpolating_node(i0, 1, value_id);
            local_node_number = my_n_p - 1 + my_n_p * i0;
            local_node_pt =
              this->interpolating_node_pt(local_node_number, value_id);
            break;

          case W:
            s_fraction[0] = 0.0;
            s_fraction[1] =
              local_one_d_fraction_of_interpolating_node(i0, 1, value_id);
            local_node_number = my_n_p * i0;
            local_node_pt =
              this->interpolating_node_pt(local_node_number, value_id);
            break;

          default:
            throw OomphLibError("my_edge not N, S, W, E\n",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        // Add node to vector of dependent element nodes
        dependent_node_number.push_back(local_node_number);
        dependent_node_pt.push_back(local_node_pt);

        // Store node's local fraction
        dependent_node_s_fraction.push_back(s_fraction);
      }

      // Store the number of dependent and master nodes for use later
      const unsigned n_dependent_nodes = dependent_node_pt.size();
      const unsigned n_master_nodes = master_node_pt.size();
      const unsigned dependent_element_nnode_1d = my_n_p;
      const unsigned master_element_nnode_1d = neigh_n_p;

      // Storage for master shapes
      Shape master_shapes(neigh_obj_pt->ninterpolating_node(value_id));

      // Get master and dependent nodal positions
      //-------------------------------------
      Vector<double> dependent_nodal_position;
      Vector<double> dependent_weight(dependent_element_nnode_1d);
      Orthpoly::gll_nodes(
        dependent_element_nnode_1d, dependent_nodal_position, dependent_weight);
      Vector<double> master_nodal_position;
      Vector<double> master_weight(master_element_nnode_1d);
      Orthpoly::gll_nodes(
        master_element_nnode_1d, master_nodal_position, master_weight);

      // Apply the vertex matching condition
      //------------------------------------
      // Vertiex matching is ensured automatically in cases where there is a
      // node at each end of the mortar that is shared between the master and
      // dependent elements. Where this is not the case, the vertex matching
      // condition must be enforced by constraining the value of the unknown at
      // the node on the dependent side to be the same as the value at that
      // location in the master.

      // Store positions of the mortar vertex/non-vertex nodes in the dependent
      // element
      const unsigned n_mortar_vertices = 2;
      Vector<unsigned> vertex_pos(n_mortar_vertices);
      vertex_pos[0] = 0;
      vertex_pos[1] = this->ninterpolating_node_1d(value_id) - 1;
      Vector<unsigned> non_vertex_pos(my_n_p - n_mortar_vertices);
      for (unsigned i = 0; i < my_n_p - n_mortar_vertices; i++)
      {
        non_vertex_pos[i] = i + 1;
      }

      // Check if the mortar vertices are shared
      for (unsigned v = 0; v < n_mortar_vertices; v++)
      {
        // Search master node storage for the node
        bool node_is_shared = false;
        for (unsigned i = 0; i < master_node_pt.size(); i++)
        {
          if (dependent_node_pt[vertex_pos[v]] == master_node_pt[i])
          {
            node_is_shared = true;
            break;
          }
        }

        // If the node is not shared then we must constrain its value by setting
        // up a hanging scheme
        if (!node_is_shared)
        {
          // Calculate weights. These are just the master shapes evaluated at
          // this dependent node's position

          // Work out this node's location in the master
          Vector<double> s_in_neigh(2);
          for (unsigned i = 0; i < 2; i++)
          {
            s_in_neigh[i] =
              s_lo_neigh[i] +
              dependent_node_s_fraction[vertex_pos[v]][translate_s[i]] *
                (s_hi_neigh[i] - s_lo_neigh[i]);
          }

          // Get master shapes at dependent nodal positions
          neigh_obj_pt->interpolating_basis(
            s_in_neigh, master_shapes, value_id);

          // Remove any existing hanging node info
          // (This may be a bit wasteful, but guarantees correctness)
          dependent_node_pt[vertex_pos[v]]->set_nonhanging();

          // Set up hanging scheme for this node
          HangInfo* explicit_hang_pt = new HangInfo(n_master_nodes);
          for (unsigned m = 0; m < n_master_nodes; m++)
          {
            explicit_hang_pt->set_master_node_pt(
              m, master_node_pt[m], master_shapes[master_node_number[m]]);
          }

          // Make node hang
          dependent_node_pt[vertex_pos[v]]->set_hanging_pt(explicit_hang_pt,
                                                           -1);
        }
      }

      // Check that there are actually dependent nodes for which we still need
      // to construct a hanging scheme. If not then there is nothing more to do.
      if (n_dependent_nodes - n_mortar_vertices > 0)
      {
        // Assemble mass matrix for mortar
        //--------------------------------
        Vector<double> psi(n_dependent_nodes - n_mortar_vertices);
        Vector<double> diag_M(n_dependent_nodes - n_mortar_vertices);
        Vector<Vector<double>> shared_node_M(n_mortar_vertices);
        for (unsigned i = 0; i < shared_node_M.size(); i++)
        {
          shared_node_M[i].resize(n_dependent_nodes - n_mortar_vertices);
        }

        // Fill in part corresponding to dependent nodal positions (unknown)
        for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
        {
          // Use L'Hosptal's rule:
          psi[i] =
            pow(-1.0, int((dependent_element_nnode_1d - 1) - i - 1)) *
            -Orthpoly::ddlegendre(dependent_element_nnode_1d - 1,
                                  dependent_nodal_position[non_vertex_pos[i]]);
          // Put in contribution
          diag_M[i] = psi[i] * dependent_weight[non_vertex_pos[i]];
        }

        // Fill in part corresponding to dependent element vertices (known)
        for (unsigned v = 0; v < shared_node_M.size(); v++)
        {
          for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
          {
            // Check if denominator is zero
            if (std::fabs(dependent_nodal_position[non_vertex_pos[i]] -
                          dependent_nodal_position[vertex_pos[v]]) >= 1.0e-8)
            {
              // We're ok
              psi[i] =
                pow(-1.0, int((dependent_element_nnode_1d - 1) - i - 1)) *
                Orthpoly::dlegendre(dependent_element_nnode_1d - 1,
                                    dependent_nodal_position[vertex_pos[v]]) /
                (dependent_nodal_position[non_vertex_pos[i]] -
                 dependent_nodal_position[vertex_pos[v]]);
            }
            // Check if numerator is zero
            else if (std::fabs(Orthpoly::dlegendre(
                       dependent_element_nnode_1d - 1,
                       dependent_nodal_position[vertex_pos[v]])) < 1.0e-8)
            {
              // We can use l'hopital's rule
              psi[i] =
                pow(-1.0, int((dependent_element_nnode_1d - 1) - i - 1)) *
                -Orthpoly::ddlegendre(
                  dependent_element_nnode_1d - 1,
                  dependent_nodal_position[non_vertex_pos[i]]);
            }
            else
            {
              // We can't use l'hopital's rule
              throw OomphLibError(
                "Cannot use l'Hopital's rule. Dividing by zero is not allowed!",
                "PRefineableQElement<2,INITIAL_NNODE_1D>::quad_hang_helper()",
                OOMPH_EXCEPTION_LOCATION);
            }
            // Put in contribution
            shared_node_M[v][i] = psi[i] * dependent_weight[vertex_pos[v]];
          }
        }

        // Assemble local projection matrix for mortar
        //--------------------------------------------

        // Have only one local projection matrix because there is just one
        // master
        Vector<Vector<double>> P(n_dependent_nodes - n_mortar_vertices);
        for (unsigned i = 0; i < P.size(); i++)
        {
          P[i].resize(n_master_nodes, 0.0);
        }

        // Storage for local coordinate
        Vector<double> s(2);

        // Sum contributions from master element shapes (quadrature).
        // The order of quadrature must be high enough to evaluate a polynomial
        // of order N_s + N_m - 3 exactly, where N_s = n_dependent_nodes, N_m =
        // n_master_nodes.
        // (Use pointers for the quadrature knots and weights so that
        // data is not unnecessarily copied)
        // unsigned quadrature_order =
        // std::max(dependent_element_nnode_1d,master_element_nnode_1d);
        Vector<double>*quadrature_knot, *quadrature_weight;
        if (dependent_element_nnode_1d >= master_element_nnode_1d)
        {
          // Use the same quadrature order as the dependent element (me)
          quadrature_knot = &dependent_nodal_position;
          quadrature_weight = &dependent_weight;
        }
        else
        {
          // Use the same quadrature order as the master element (neighbour)
          quadrature_knot = &master_nodal_position;
          quadrature_weight = &master_weight;
        }

        // Quadrature loop
        for (unsigned q = 0; q < (*quadrature_weight).size(); q++)
        {
          // Evaluate mortar test functions at the quadrature knots in the
          // dependent
          // s[active_coord_index] = (*quadrature_knot)[q];
          Vector<double> s_on_mortar(1);
          s_on_mortar[0] = (*quadrature_knot)[q];

          // Get psi
          for (unsigned k = 0; k < n_dependent_nodes - n_mortar_vertices; k++)
          {
            // Check if denominator is zero
            if (std::fabs(dependent_nodal_position[non_vertex_pos[k]] -
                          s_on_mortar[0]) >= 1.0e-08)
            {
              // We're ok
              psi[k] =
                pow(-1.0, int((dependent_element_nnode_1d - 1) - k - 1)) *
                Orthpoly::dlegendre(dependent_element_nnode_1d - 1,
                                    s_on_mortar[0]) /
                (dependent_nodal_position[non_vertex_pos[k]] - s_on_mortar[0]);
            }
            // Check if numerator is zero
            else if (std::fabs(Orthpoly::dlegendre(
                       dependent_element_nnode_1d - 1, s_on_mortar[0])) <
                     1.0e-8)
            {
              // We can use l'Hopital's rule
              psi[k] =
                pow(-1.0, int((dependent_element_nnode_1d - 1) - k - 1)) *
                -Orthpoly::ddlegendre(dependent_element_nnode_1d - 1,
                                      s_on_mortar[0]);
            }
            else
            {
              // We can't use l'hopital's rule
              throw OomphLibError(
                "Cannot use l'Hopital's rule. Dividing by zero is not allowed!",
                "PRefineableQElement<2,INITIAL_NNODE_1D>::quad_hang_helper()",
                OOMPH_EXCEPTION_LOCATION);
            }
          }

          // Convert coordinate on mortar to local fraction in dependent element
          Vector<double> s_fraction(2);
          for (unsigned i = 0; i < 2; i++)
          {
            s_fraction[i] = (i == active_coord_index) ?
                              0.5 * (s_on_mortar[0] + 1.0) :
                              dependent_node_s_fraction[vertex_pos[0]][i];
          }

          // Project active coordinate into master element
          Vector<double> s_in_neigh(2);
          for (unsigned i = 0; i < 2; i++)
          {
            s_in_neigh[i] = s_lo_neigh[i] + s_fraction[translate_s[i]] *
                                              (s_hi_neigh[i] - s_lo_neigh[i]);
          }

          // Evaluate master shapes at projections of local quadrature knots
          neigh_obj_pt->interpolating_basis(
            s_in_neigh, master_shapes, value_id);

          // Populate local projection matrix
          for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
          {
            for (unsigned j = 0; j < n_master_nodes; j++)
            {
              P[i][j] += master_shapes[master_node_number[j]] * psi[i] *
                         (*quadrature_weight)[q];
            }
          }
        }

        // Assemble global projection matrices for mortar
        //-----------------------------------------------
        // Need to subtract contributions from the "known unknowns"
        // corresponding to the nodes at the vertices of the mortar

        // Assemble contributions from mortar vertex nodes
        for (unsigned v = 0; v < n_mortar_vertices; v++)
        {
          // Convert coordinate on mortar to local fraction in dependent element
          Vector<double> s_fraction(2);
          for (unsigned i = 0; i < 2; i++)
          {
            s_fraction[i] =
              (i == active_coord_index) ?
                0.5 * (dependent_nodal_position[vertex_pos[v]] + 1.0) :
                dependent_node_s_fraction[vertex_pos[0]][i];
          }

          // Project active coordinate into master element
          Vector<double> s_in_neigh(2);
          for (unsigned i = 0; i < 2; i++)
          {
            s_in_neigh[i] = s_lo_neigh[i] + s_fraction[translate_s[i]] *
                                              (s_hi_neigh[i] - s_lo_neigh[i]);
          }

          // Get master shapes at dependent nodal positions
          neigh_obj_pt->interpolating_basis(
            s_in_neigh, master_shapes, value_id);

          for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
          {
            for (unsigned k = 0; k < n_master_nodes; k++)
            {
              P[i][k] -=
                master_shapes[master_node_number[k]] * shared_node_M[v][i];
            }
          }
        }

        // Solve mortar system
        //--------------------
        for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
        {
          for (unsigned j = 0; j < n_master_nodes; j++)
          {
            P[i][j] /= diag_M[i];
          }
        }

        // Create structures to hold the hanging info
        //-------------------------------------------
        Vector<HangInfo*> hang_info_pt(n_dependent_nodes);
        for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
        {
          hang_info_pt[i] = new HangInfo(n_master_nodes);
        }

        // Copy information to hanging nodes
        //----------------------------------
        for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
        {
          for (unsigned j = 0; j < n_master_nodes; j++)
          {
            hang_info_pt[i]->set_master_node_pt(j, master_node_pt[j], P[i][j]);
          }
        }

        // Set pointers to hanging info
        //-----------------------------
        for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
        {
          // Check that the node shoule actually hang.
          // That is, if the polynomial orders of the elements at a p-type
          // non-conormity are both odd then the middle node on the edge is
          // shared but a hanging scheme has been constructed for it.
          bool node_is_really_shared = false;
          for (unsigned m = 0; m < hang_info_pt[i]->nmaster(); m++)
          {
            // We can simply check if the hanging scheme lists itself as a
            // master
            if (hang_info_pt[i]->master_node_pt(m) ==
                dependent_node_pt[non_vertex_pos[i]])
            {
              node_is_really_shared = true;

#ifdef PARANOID
              // Also check the corresponding weight: it should be 1
              if (std::fabs(hang_info_pt[i]->master_weight(m) - 1.0) > 1.0e-06)
              {
                throw OomphLibError(
                  "Something fishy here -- with shared node at a mortar vertex",
                  "PRefineableQElemen<2,INITIAL_NNODE_1D>t::quad_hang_helper()",
                  OOMPH_EXCEPTION_LOCATION);
              }
#endif
            }
          }

          // Now we can make the node hang, if it isn't a shared node
          if (!node_is_really_shared)
          {
            dependent_node_pt[non_vertex_pos[i]]->set_hanging_pt(
              hang_info_pt[i], -1);
          }
        }

      } // End of case where there are still dependent nodes

#ifdef PARANOID
      // Check all dependent nodes, hanging or otherwise
      for (unsigned i = 0; i < n_dependent_nodes; i++)
      {
        // Check that weights sum to 1 for those that hang
        if (dependent_node_pt[i]->is_hanging())
        {
          // Check that weights sum to 1
          double weight_sum = 0.0;
          for (unsigned m = 0;
               m < dependent_node_pt[i]->hanging_pt()->nmaster();
               m++)
          {
            weight_sum += dependent_node_pt[i]->hanging_pt()->master_weight(m);
          }

          // Warn if not
          if (std::fabs(weight_sum - 1.0) > 1.0e-08)
          {
            oomph_info << "Sum of master weights: " << weight_sum << std::endl;
            OomphLibWarning(
              "Weights in hanging scheme do not sum to 1",
              "PRefineableQElement<2,INITIAL_NNODE_1D>::quad_hang_helper()",
              OOMPH_EXCEPTION_LOCATION);
          }
        }
        else
        {
          // Check that this node is shared with the master element if it
          // isn't hanging
          bool is_master = false;
          for (unsigned n = 0; n < n_master_nodes; n++)
          {
            if (dependent_node_pt[i] == master_node_pt[n])
            {
              // Node is a master
              is_master = true;
              break;
            }
          }

          if (!is_master)
          {
            // Throw error
            std::ostringstream error_string;
            error_string << "This node in the dependent element is neither"
                         << std::endl
                         << "hanging or shared with a master element."
                         << std::endl;

            throw OomphLibError(
              error_string.str(),
              "PRefineableQElement<2,INITIAL_NNODE_1D>::quad_hang_helper()",
              OOMPH_EXCEPTION_LOCATION);
          }
        }
      }
#endif

      // Finally, Loop over all dependent nodes and fine-tune their positions
      //-----------------------------------------------------------------
      // Here we simply set the node's positions to be consistent
      // with the hanging scheme. This is not strictly necessary
      // because it is done in the mesh adaptation before the node
      // becomes non-hanging later on. We make no attempt to ensure
      // (strong) continuity in the position across the mortar.
      for (unsigned i = 0; i < n_dependent_nodes; i++)
      {
        // Only fine-tune hanging nodes
        if (dependent_node_pt[i]->is_hanging())
        {
          // If we are doing the position, then
          if (value_id == -1)
          {
            // Get the local coordinate of this dependent node
            Vector<double> s_local(2);
            this->local_coordinate_of_node(dependent_node_number[i], s_local);

            // Get the position from interpolation in this element via
            // the hanging scheme
            Vector<double> x_in_neighb(2);
            this->interpolated_x(s_local, x_in_neighb);

            // Fine adjust the coordinates (macro map will pick up boundary
            // accurately but will lead to different element edges)
            dependent_node_pt[i]->x(0) = x_in_neighb[0];
            dependent_node_pt[i]->x(1) = x_in_neighb[1];
          }
        }
      }
    } // End of case where this interface is to be mortared
  }

  /// /////////////////////////////////////////////////////////////
  //       3D PRefineableQElements
  /// /////////////////////////////////////////////////////////////

  /// Get local coordinates of node j in the element; vector sets its own size
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<3, INITIAL_NNODE_1D>::local_coordinate_of_node(
    const unsigned& n, Vector<double>& s) const
  {
    s.resize(3);
    unsigned Nnode_1d = this->nnode_1d();
    unsigned n0 = n % Nnode_1d;
    unsigned n1 = unsigned(double(n) / double(Nnode_1d)) % Nnode_1d;
    unsigned n2 = unsigned(double(n) / double(Nnode_1d * Nnode_1d));

    switch (Nnode_1d)
    {
      case 2:
        OneDimensionalLegendreShape<2>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<2>::nodal_position(n0);
        s[1] = OneDimensionalLegendreShape<2>::nodal_position(n1);
        s[2] = OneDimensionalLegendreShape<2>::nodal_position(n2);
        break;
      case 3:
        OneDimensionalLegendreShape<3>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<3>::nodal_position(n0);
        s[1] = OneDimensionalLegendreShape<3>::nodal_position(n1);
        s[2] = OneDimensionalLegendreShape<3>::nodal_position(n2);
        break;
      case 4:
        OneDimensionalLegendreShape<4>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<4>::nodal_position(n0);
        s[1] = OneDimensionalLegendreShape<4>::nodal_position(n1);
        s[2] = OneDimensionalLegendreShape<4>::nodal_position(n2);
        break;
      case 5:
        OneDimensionalLegendreShape<5>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<5>::nodal_position(n0);
        s[1] = OneDimensionalLegendreShape<5>::nodal_position(n1);
        s[2] = OneDimensionalLegendreShape<5>::nodal_position(n2);
        break;
      case 6:
        OneDimensionalLegendreShape<6>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<6>::nodal_position(n0);
        s[1] = OneDimensionalLegendreShape<6>::nodal_position(n1);
        s[2] = OneDimensionalLegendreShape<6>::nodal_position(n2);
        break;
      case 7:
        OneDimensionalLegendreShape<7>::calculate_nodal_positions();
        s[0] = OneDimensionalLegendreShape<7>::nodal_position(n0);
        s[1] = OneDimensionalLegendreShape<7>::nodal_position(n1);
        s[2] = OneDimensionalLegendreShape<7>::nodal_position(n2);
        break;
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       shape functions.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        break;
    }
  }

  /// Get the local fractino of node j in the element
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<3, INITIAL_NNODE_1D>::local_fraction_of_node(
    const unsigned& n, Vector<double>& s_fraction)
  {
    this->local_coordinate_of_node(n, s_fraction);
    s_fraction[0] = 0.5 * (s_fraction[0] + 1.0);
    s_fraction[1] = 0.5 * (s_fraction[1] + 1.0);
    s_fraction[2] = 0.5 * (s_fraction[2] + 1.0);
  }

  template<unsigned INITIAL_NNODE_1D>
  double PRefineableQElement<3, INITIAL_NNODE_1D>::local_one_d_fraction_of_node(
    const unsigned& n1d, const unsigned& i)
  {
    switch (this->nnode_1d())
    {
      case 2:
        OneDimensionalLegendreShape<2>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<2>::nodal_position(n1d) + 1.0);
      case 3:
        OneDimensionalLegendreShape<3>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<3>::nodal_position(n1d) + 1.0);
      case 4:
        OneDimensionalLegendreShape<4>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<4>::nodal_position(n1d) + 1.0);
      case 5:
        OneDimensionalLegendreShape<5>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<5>::nodal_position(n1d) + 1.0);
      case 6:
        OneDimensionalLegendreShape<6>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<6>::nodal_position(n1d) + 1.0);
      case 7:
        OneDimensionalLegendreShape<7>::calculate_nodal_positions();
        return 0.5 *
               (OneDimensionalLegendreShape<7>::nodal_position(n1d) + 1.0);
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       shape functions.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        return 0.0;
    }
  }

  //==================================================================
  /// Return the node at the specified local coordinate
  //==================================================================
  template<unsigned INITIAL_NNODE_1D>
  Node* PRefineableQElement<3, INITIAL_NNODE_1D>::get_node_at_local_coordinate(
    const Vector<double>& s) const
  {
    unsigned Nnode_1d = this->nnode_1d();
    // Load the tolerance into a local variable
    double tol = FiniteElement::Node_location_tolerance;
    // There are two possible indices.
    Vector<int> index(3);

    // Loop over indices
    for (unsigned i = 0; i < 3; i++)
    {
      // Determine the index
      // -------------------

      bool is_found = false;

      // If we are at the lower limit, the index is zero
      if (std::fabs(s[i] + 1.0) < tol)
      {
        index[i] = 0;
        is_found = true;
      }
      // If we are at the upper limit, the index is the number of nodes minus 1
      else if (std::fabs(s[i] - 1.0) < tol)
      {
        index[i] = Nnode_1d - 1;
        is_found = true;
      }
      // Otherwise, we have to calculate the index in general
      else
      {
        // Compute Gauss-Lobatto-Legendre node positions
        Vector<double> z;
        Orthpoly::gll_nodes(Nnode_1d, z);
        // Loop over possible internal nodal positions
        for (unsigned n = 1; n < Nnode_1d - 1; n++)
        {
          if (std::fabs(z[n] - s[i]) < tol)
          {
            index[i] = n;
            is_found = true;
            break;
          }
        }
      }

      if (!is_found)
      {
        // No matching nodes
        return 0;
      }
    }
    // If we've got here we have a node, so let's return a pointer to it
    return this->node_pt(index[0] + Nnode_1d * index[1] +
                         Nnode_1d * Nnode_1d * index[2]);
  }

  //===================================================================
  /// If a neighbouring element has already created a node at
  /// a position corresponding to the local fractional position within the
  /// present element, s_fraction, return
  /// a pointer to that node. If not, return NULL (0).
  //===================================================================
  template<unsigned INITIAL_NNODE_1D>
  Node* PRefineableQElement<3, INITIAL_NNODE_1D>::node_created_by_neighbour(
    const Vector<double>& s_fraction, bool& is_periodic)
  {
    using namespace OcTreeNames;

    // Calculate the faces/edges on which the node lies
    Vector<int> faces;
    Vector<int> edges;

    if (s_fraction[0] == 0.0)
    {
      faces.push_back(L);
      if (s_fraction[1] == 0.0)
      {
        edges.push_back(LD);
      }
      if (s_fraction[2] == 0.0)
      {
        edges.push_back(LB);
      }
      if (s_fraction[1] == 1.0)
      {
        edges.push_back(LU);
      }
      if (s_fraction[2] == 1.0)
      {
        edges.push_back(LF);
      }
    }

    if (s_fraction[0] == 1.0)
    {
      faces.push_back(R);
      if (s_fraction[1] == 0.0)
      {
        edges.push_back(RD);
      }
      if (s_fraction[2] == 0.0)
      {
        edges.push_back(RB);
      }
      if (s_fraction[1] == 1.0)
      {
        edges.push_back(RU);
      }
      if (s_fraction[2] == 1.0)
      {
        edges.push_back(RF);
      }
    }

    if (s_fraction[1] == 0.0)
    {
      faces.push_back(D);
      if (s_fraction[2] == 0.0)
      {
        edges.push_back(DB);
      }
      if (s_fraction[2] == 1.0)
      {
        edges.push_back(DF);
      }
    }

    if (s_fraction[1] == 1.0)
    {
      faces.push_back(U);
      if (s_fraction[2] == 0.0)
      {
        edges.push_back(UB);
      }
      if (s_fraction[2] == 1.0)
      {
        edges.push_back(UF);
      }
    }

    if (s_fraction[2] == 0.0)
    {
      faces.push_back(B);
    }

    if (s_fraction[2] == 1.0)
    {
      faces.push_back(F);
    }

    // Find the number of faces
    unsigned n_face = faces.size();

    // Find the number of edges
    unsigned n_edge = edges.size();

    Vector<unsigned> translate_s(3);
    Vector<double> s_lo_neigh(3);
    Vector<double> s_hi_neigh(3);
    Vector<double> s(3);

    int neigh_face, diff_level;
    bool in_neighbouring_tree;

    // Loop over the faces on which the node lies
    //------------------------------------------
    for (unsigned j = 0; j < n_face; j++)
    {
      // Find pointer to neighbouring element along face
      OcTree* neigh_pt;
      neigh_pt = octree_pt()->gteq_face_neighbour(faces[j],
                                                  translate_s,
                                                  s_lo_neigh,
                                                  s_hi_neigh,
                                                  neigh_face,
                                                  diff_level,
                                                  in_neighbouring_tree);

      // Neighbour exists
      if (neigh_pt != 0)
      {
        // Have any of its vertex nodes been created yet?
        // (Must look in incomplete neighbours because after the
        // pre-build they may contain pointers to the required nodes. e.g.
        // h-refinement of neighbouring linear and quadratic elements)
        bool a_vertex_node_is_built = false;
        QElement<3, INITIAL_NNODE_1D>* neigh_obj_pt =
          dynamic_cast<QElement<3, INITIAL_NNODE_1D>*>(neigh_pt->object_pt());
        if (neigh_obj_pt == 0)
        {
          throw OomphLibError("Not a quad element!",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
        for (unsigned vnode = 0; vnode < neigh_obj_pt->nvertex_node(); vnode++)
        {
          if (neigh_obj_pt->vertex_node_pt(vnode) != 0)
            a_vertex_node_is_built = true;
          break;
        }
        if (a_vertex_node_is_built)
        {
          // We now need to translate the nodal location, defined in terms
          // of the fractional coordinates of the present element into
          // those of its neighbour. For this we use the information returned
          // to use from the octree function.

          // Calculate the local coordinate in the neighbour
          // Note that we need to use the translation scheme in case
          // the local coordinates are swapped in the neighbour.
          for (unsigned i = 0; i < 3; i++)
          {
            s[i] = s_lo_neigh[i] +
                   s_fraction[translate_s[i]] * (s_hi_neigh[i] - s_lo_neigh[i]);
          }

          // Find the node in the neighbour
          Node* neighbour_node_pt =
            neigh_pt->object_pt()->get_node_at_local_coordinate(s);

          // If there is a node, return it
          if (neighbour_node_pt != 0)
          {
            // Now work out whether it's a periodic boundary
            // only possible if we have moved into a neighbouring tree
            if (in_neighbouring_tree)
            {
              // Return whether the neighbour is periodic
              is_periodic =
                octree_pt()->root_pt()->is_neighbour_periodic(faces[j]);
            }

            return neighbour_node_pt;
          }
        }
      }
    } // End of loop over faces


    // Loop over the edges on which the node lies
    //------------------------------------------
    for (unsigned j = 0; j < n_edge; j++)
    {
      // Even if we restrict ourselves to true edge neighbours (i.e.
      // elements that are not also face neighbours) there may be multiple
      // edge neighbours across the edges between multiple root octrees.
      // When making the first call to OcTree::gteq_true_edge_neighbour(...)
      // we simply return the first one of these multiple edge neighbours
      // (if there are any at all, of course) and also return the total number
      // of true edge neighbours. If the node in question already exists
      // on the first edge neighbour we're done. If it doesn't it may exist
      // on other edge neighbours so we repeat the process over all
      // other edge neighbours (bailing out if a node is found, of course).

      // Initially return the zero-th true edge neighbour
      unsigned i_root_edge_neighbour = 0;

      // Initialise the total number of true edge neighbours
      unsigned nroot_edge_neighbour = 0;

      // Keep searching until we've found the node or until we've checked
      // all available edge neighbours
      bool keep_searching = true;
      while (keep_searching)
      {
        // Find pointer to neighbouring element along edge
        OcTree* neigh_pt;
        neigh_pt = octree_pt()->gteq_true_edge_neighbour(edges[j],
                                                         i_root_edge_neighbour,
                                                         nroot_edge_neighbour,
                                                         translate_s,
                                                         s_lo_neigh,
                                                         s_hi_neigh,
                                                         neigh_face,
                                                         diff_level);

        // Neighbour exists
        if (neigh_pt != 0)
        {
          // Have any of its vertex nodes been created yet?
          // (Must look in incomplete neighbours because after the
          // pre-build they may contain pointers to the required nodes. e.g.
          // h-refinement of neighbouring linear and quadratic elements)
          bool a_vertex_node_is_built = false;
          QElement<3, INITIAL_NNODE_1D>* neigh_obj_pt =
            dynamic_cast<QElement<3, INITIAL_NNODE_1D>*>(neigh_pt->object_pt());
          if (neigh_obj_pt == 0)
          {
            throw OomphLibError("Not a quad element!",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
          for (unsigned vnode = 0; vnode < neigh_obj_pt->nvertex_node();
               vnode++)
          {
            if (neigh_obj_pt->vertex_node_pt(vnode) != 0)
              a_vertex_node_is_built = true;
            break;
          }
          if (a_vertex_node_is_built)
          {
            // We now need to translate the nodal location, defined in terms
            // of the fractional coordinates of the present element into
            // those of its neighbour. For this we use the information returned
            // to use from the octree function.

            // Calculate the local coordinate in the neighbour
            // Note that we need to use the translation scheme in case
            // the local coordinates are swapped in the neighbour.
            for (unsigned i = 0; i < 3; i++)
            {
              s[i] = s_lo_neigh[i] + s_fraction[translate_s[i]] *
                                       (s_hi_neigh[i] - s_lo_neigh[i]);
            }

            // Find the node in the neighbour
            Node* neighbour_node_pt =
              neigh_pt->object_pt()->get_node_at_local_coordinate(s);

            // If there is a node, return it
            if (neighbour_node_pt != 0)
            {
              // Get the faces on which the edge lies
              Vector<int> faces_attached_to_edge =
                OcTree::faces_of_common_edge(edges[j]);

              // Get the number of entries in the vector
              unsigned n_faces_attached_to_edge = faces_attached_to_edge.size();

              // Loop over the faces
              for (unsigned i_face = 0; i_face < n_faces_attached_to_edge;
                   i_face++)
              {
                // Is the node periodic in the face direction?
                is_periodic = octree_pt()->root_pt()->is_neighbour_periodic(
                  faces_attached_to_edge[i_face]);

                // Check if the edge is periodic in the i_face-th face direction
                if (is_periodic)
                {
                  // We're done!
                  break;
                }
              } // for (unsigned
                // i_face=0;i_face<n_faces_attached_to_edge;i_face++)

              return neighbour_node_pt;
            }
          }
        }

        // Keep searching, but only if there are further edge neighbours
        // Try next root edge neighbour
        i_root_edge_neighbour++;

        // Have we exhausted the search
        if (i_root_edge_neighbour >= nroot_edge_neighbour)
        {
          keep_searching = false;
        }

      } // End of while keep searching over all true edge neighbours

    } // End of loop over edges

    // Node not found, return null
    return 0;
  }

  //===================================================================
  /// If a neighbouring element's son has already created a node at
  /// a position corresponding to the local fractional position within the
  /// present element, s_fraction, return
  /// a pointer to that node. If not, return NULL (0). If the node is
  /// periodic the flag is_periodic will be true
  //===================================================================
  template<unsigned INITIAL_NNODE_1D>
  Node* PRefineableQElement<3, INITIAL_NNODE_1D>::
    node_created_by_son_of_neighbour(const Vector<double>& s_fraction,
                                     bool& is_periodic)
  {
    using namespace OcTreeNames;

    // Calculate the faces/edges on which the node lies
    Vector<int> faces;
    Vector<int> edges;

    if (s_fraction[0] == 0.0)
    {
      faces.push_back(L);
      if (s_fraction[1] == 0.0)
      {
        edges.push_back(LD);
      }
      if (s_fraction[2] == 0.0)
      {
        edges.push_back(LB);
      }
      if (s_fraction[1] == 1.0)
      {
        edges.push_back(LU);
      }
      if (s_fraction[2] == 1.0)
      {
        edges.push_back(LF);
      }
    }

    if (s_fraction[0] == 1.0)
    {
      faces.push_back(R);
      if (s_fraction[1] == 0.0)
      {
        edges.push_back(RD);
      }
      if (s_fraction[2] == 0.0)
      {
        edges.push_back(RB);
      }
      if (s_fraction[1] == 1.0)
      {
        edges.push_back(RU);
      }
      if (s_fraction[2] == 1.0)
      {
        edges.push_back(RF);
      }
    }

    if (s_fraction[1] == 0.0)
    {
      faces.push_back(D);
      if (s_fraction[2] == 0.0)
      {
        edges.push_back(DB);
      }
      if (s_fraction[2] == 1.0)
      {
        edges.push_back(DF);
      }
    }

    if (s_fraction[1] == 1.0)
    {
      faces.push_back(U);
      if (s_fraction[2] == 0.0)
      {
        edges.push_back(UB);
      }
      if (s_fraction[2] == 1.0)
      {
        edges.push_back(UF);
      }
    }

    if (s_fraction[2] == 0.0)
    {
      faces.push_back(B);
    }

    if (s_fraction[2] == 1.0)
    {
      faces.push_back(F);
    }

    // Find the number of faces
    unsigned n_face = faces.size();

    // Find the number of edges
    unsigned n_edge = edges.size();

    Vector<unsigned> translate_s(3);
    Vector<double> s_lo_neigh(3);
    Vector<double> s_hi_neigh(3);
    Vector<double> s(3);

    int neigh_face, diff_level;
    bool in_neighbouring_tree;

    // Loop over the faces on which the node lies
    //------------------------------------------
    for (unsigned j = 0; j < n_face; j++)
    {
      // Find pointer to neighbouring element along face
      OcTree* neigh_pt;
      neigh_pt = octree_pt()->gteq_face_neighbour(faces[j],
                                                  translate_s,
                                                  s_lo_neigh,
                                                  s_hi_neigh,
                                                  neigh_face,
                                                  diff_level,
                                                  in_neighbouring_tree);

      // Neighbour exists
      if (neigh_pt != 0)
      {
        // Have its nodes been created yet?
        // (Must look in sons of unfinished neighbours too!!!)
        if (true)
        {
          // We now need to translate the nodal location, defined in terms
          // of the fractional coordinates of the present element into
          // those of its neighbour. For this we use the information returned
          // to use from the octree function.

          // Calculate the local coordinate in the neighbour
          // Note that we need to use the translation scheme in case
          // the local coordinates are swapped in the neighbour.
          for (unsigned i = 0; i < 3; i++)
          {
            s[i] = s_lo_neigh[i] +
                   s_fraction[translate_s[i]] * (s_hi_neigh[i] - s_lo_neigh[i]);
          }

          // Check if the element has sons
          if (neigh_pt->nsons() != 0)
          {
            // First, find the son element in which the node should live

            // Find coordinates in the sons
            Vector<double> s_in_son(3);
            int son = -10;
            // On the left
            if (s[0] < 0.0)
            {
              // On the bottom
              if (s[1] < 0.0)
              {
                // On the back
                if (s[2] < 0.0)
                {
                  // It's the LDB son
                  son = OcTreeNames::LDB;
                  s_in_son[0] = 1.0 + 2.0 * s[0];
                  s_in_son[1] = 1.0 + 2.0 * s[1];
                  s_in_son[2] = 1.0 + 2.0 * s[2];
                }
                // On the front
                else
                {
                  // It's the LDF son
                  son = OcTreeNames::LDF;
                  s_in_son[0] = 1.0 + 2.0 * s[0];
                  s_in_son[1] = 1.0 + 2.0 * s[1];
                  s_in_son[2] = -1.0 + 2.0 * s[2];
                }
              }
              // On the top
              else
              {
                // On the back
                if (s[2] < 0.0)
                {
                  // It's the LUB son
                  son = OcTreeNames::LUB;
                  s_in_son[0] = 1.0 + 2.0 * s[0];
                  s_in_son[1] = -1.0 + 2.0 * s[1];
                  s_in_son[2] = 1.0 + 2.0 * s[2];
                }
                // On the front
                else
                {
                  // It's the LUF son
                  son = OcTreeNames::LUF;
                  s_in_son[0] = 1.0 + 2.0 * s[0];
                  s_in_son[1] = -1.0 + 2.0 * s[1];
                  s_in_son[2] = -1.0 + 2.0 * s[2];
                }
              }
            }
            // On the right
            else
            {
              // On the bottom
              if (s[1] < 0.0)
              {
                // On the back
                if (s[2] < 0.0)
                {
                  // It's the RDB son
                  son = OcTreeNames::RDB;
                  s_in_son[0] = -1.0 + 2.0 * s[0];
                  s_in_son[1] = 1.0 + 2.0 * s[1];
                  s_in_son[2] = 1.0 + 2.0 * s[2];
                }
                // On the front
                else
                {
                  // It's the RDF son
                  son = OcTreeNames::RDF;
                  s_in_son[0] = -1.0 + 2.0 * s[0];
                  s_in_son[1] = 1.0 + 2.0 * s[1];
                  s_in_son[2] = -1.0 + 2.0 * s[2];
                }
              }
              // On the top
              else
              {
                // On the back
                if (s[2] < 0.0)
                {
                  // It's the RUB son
                  son = OcTreeNames::RUB;
                  s_in_son[0] = -1.0 + 2.0 * s[0];
                  s_in_son[1] = -1.0 + 2.0 * s[1];
                  s_in_son[2] = 1.0 + 2.0 * s[2];
                }
                // On the front
                else
                {
                  // It's the RUF son
                  son = OcTreeNames::RUF;
                  s_in_son[0] = -1.0 + 2.0 * s[0];
                  s_in_son[1] = -1.0 + 2.0 * s[1];
                  s_in_son[2] = -1.0 + 2.0 * s[2];
                }
              }
            }

            // Find the node in the neighbour's son
            Node* neighbour_son_node_pt =
              neigh_pt->son_pt(son)->object_pt()->get_node_at_local_coordinate(
                s_in_son);

            // If there is a node, return it
            if (neighbour_son_node_pt != 0)
            {
              // Now work out whether it's a periodic boundary
              // only possible if we have moved into a neighbouring tree
              if (in_neighbouring_tree)
              {
                // Return whether the neighbour is periodic
                is_periodic =
                  octree_pt()->root_pt()->is_neighbour_periodic(faces[j]);
              }

              // Return the pointer to the neighbouring node
              return neighbour_son_node_pt;
            }
          }
        }
      }
    } // End of loop over faces


    // Loop over the edges on which the node lies
    //------------------------------------------
    for (unsigned j = 0; j < n_edge; j++)
    {
      // Even if we restrict ourselves to true edge neighbours (i.e.
      // elements that are not also face neighbours) there may be multiple
      // edge neighbours across the edges between multiple root octrees.
      // When making the first call to OcTree::gteq_true_edge_neighbour(...)
      // we simply return the first one of these multiple edge neighbours
      // (if there are any at all, of course) and also return the total number
      // of true edge neighbours. If the node in question already exists
      // on the first edge neighbour we're done. If it doesn't it may exist
      // on other edge neighbours so we repeat the process over all
      // other edge neighbours (bailing out if a node is found, of course).

      // Initially return the zero-th true edge neighbour
      unsigned i_root_edge_neighbour = 0;

      // Initialise the total number of true edge neighbours
      unsigned nroot_edge_neighbour = 0;

      // Keep searching until we've found the node or until we've checked
      // all available edge neighbours
      bool keep_searching = true;
      while (keep_searching)
      {
        // Find pointer to neighbouring element along edge
        OcTree* neigh_pt;
        neigh_pt = octree_pt()->gteq_true_edge_neighbour(edges[j],
                                                         i_root_edge_neighbour,
                                                         nroot_edge_neighbour,
                                                         translate_s,
                                                         s_lo_neigh,
                                                         s_hi_neigh,
                                                         neigh_face,
                                                         diff_level);

        // Neighbour exists
        if (neigh_pt != 0)
        {
          // Have its nodes been created yet?
          // (Must look in sons of unfinished neighbours too!!!)
          if (true)
          {
            // We now need to translate the nodal location, defined in terms
            // of the fractional coordinates of the present element into
            // those of its neighbour. For this we use the information returned
            // to use from the octree function.

            // Calculate the local coordinate in the neighbour
            // Note that we need to use the translation scheme in case
            // the local coordinates are swapped in the neighbour.
            for (unsigned i = 0; i < 3; i++)
            {
              s[i] = s_lo_neigh[i] + s_fraction[translate_s[i]] *
                                       (s_hi_neigh[i] - s_lo_neigh[i]);
            }

            // Check if the element has sons
            if (neigh_pt->nsons() != 0)
            {
              // First, find the son element in which the node should live

              // Find coordinates in the sons
              Vector<double> s_in_son(3);
              int son = -10;
              // On the left
              if (s[0] < 0.0)
              {
                // On the bottom
                if (s[1] < 0.0)
                {
                  // On the back
                  if (s[2] < 0.0)
                  {
                    // It's the LDB son
                    son = OcTreeNames::LDB;
                    s_in_son[0] = 1.0 + 2.0 * s[0];
                    s_in_son[1] = 1.0 + 2.0 * s[1];
                    s_in_son[2] = 1.0 + 2.0 * s[2];
                  }
                  // On the front
                  else
                  {
                    // It's the LDF son
                    son = OcTreeNames::LDF;
                    s_in_son[0] = 1.0 + 2.0 * s[0];
                    s_in_son[1] = 1.0 + 2.0 * s[1];
                    s_in_son[2] = -1.0 + 2.0 * s[2];
                  }
                }
                // On the top
                else
                {
                  // On the back
                  if (s[2] < 0.0)
                  {
                    // It's the LUB son
                    son = OcTreeNames::LUB;
                    s_in_son[0] = 1.0 + 2.0 * s[0];
                    s_in_son[1] = -1.0 + 2.0 * s[1];
                    s_in_son[2] = 1.0 + 2.0 * s[2];
                  }
                  // On the front
                  else
                  {
                    // It's the LUF son
                    son = OcTreeNames::LUF;
                    s_in_son[0] = 1.0 + 2.0 * s[0];
                    s_in_son[1] = -1.0 + 2.0 * s[1];
                    s_in_son[2] = -1.0 + 2.0 * s[2];
                  }
                }
              }
              // On the right
              else
              {
                // On the bottom
                if (s[1] < 0.0)
                {
                  // On the back
                  if (s[2] < 0.0)
                  {
                    // It's the RDB son
                    son = OcTreeNames::RDB;
                    s_in_son[0] = -1.0 + 2.0 * s[0];
                    s_in_son[1] = 1.0 + 2.0 * s[1];
                    s_in_son[2] = 1.0 + 2.0 * s[2];
                  }
                  // On the front
                  else
                  {
                    // It's the RDF son
                    son = OcTreeNames::RDF;
                    s_in_son[0] = -1.0 + 2.0 * s[0];
                    s_in_son[1] = 1.0 + 2.0 * s[1];
                    s_in_son[2] = -1.0 + 2.0 * s[2];
                  }
                }
                // On the top
                else
                {
                  // On the back
                  if (s[2] < 0.0)
                  {
                    // It's the RUB son
                    son = OcTreeNames::RUB;
                    s_in_son[0] = -1.0 + 2.0 * s[0];
                    s_in_son[1] = -1.0 + 2.0 * s[1];
                    s_in_son[2] = 1.0 + 2.0 * s[2];
                  }
                  // On the front
                  else
                  {
                    // It's the RUF son
                    son = OcTreeNames::RUF;
                    s_in_son[0] = -1.0 + 2.0 * s[0];
                    s_in_son[1] = -1.0 + 2.0 * s[1];
                    s_in_son[2] = -1.0 + 2.0 * s[2];
                  }
                }
              }

              // Find the node in the neighbour's son
              Node* neighbour_son_node_pt =
                neigh_pt->son_pt(son)
                  ->object_pt()
                  ->get_node_at_local_coordinate(s_in_son);

              // If there is a node, return it
              if (neighbour_son_node_pt != 0)
              {
                // Get the faces on which the edge lies
                Vector<int> faces_attached_to_edge =
                  OcTree::faces_of_common_edge(edges[j]);

                // Get the number of entries in the vector
                unsigned n_faces_attached_to_edge =
                  faces_attached_to_edge.size();

                // Loop over the faces
                for (unsigned i_face = 0; i_face < n_faces_attached_to_edge;
                     i_face++)
                {
                  // Is the node periodic in the face direction?
                  is_periodic = octree_pt()->root_pt()->is_neighbour_periodic(
                    faces_attached_to_edge[i_face]);

                  // Check if the edge is periodic in the i_face-th face
                  // direction
                  if (is_periodic)
                  {
                    // We're done!
                    break;
                  }
                } // for (unsigned
                  // i_face=0;i_face<n_faces_attached_to_edge;i_face++)

                // Return the pointer to the neighbouring node
                return neighbour_son_node_pt;
              }
            }
          }
        }

        // Keep searching, but only if there are further edge neighbours
        // Try next root edge neighbour
        i_root_edge_neighbour++;

        // Have we exhausted the search
        if (i_root_edge_neighbour >= nroot_edge_neighbour)
        {
          keep_searching = false;
        }

      } // End of while keep searching over all true edge neighbours

    } // End of loop over edges

    // Node not found, return null
    return 0;
  }

  //==================================================================
  /// Set the correct p-order of the element based on that of its
  /// father. Then construct an integration scheme of the correct order.
  /// If an adopted father is specified, information from this is
  /// used instead of using the father found from the tree.
  //==================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<3, INITIAL_NNODE_1D>::initial_setup(
    Tree* const& adopted_father_pt, const unsigned& initial_p_order)
  {
    // Storage for pointer to my father (in octree impersonation)
    OcTree* father_pt = 0;

    // Check if an adopted father has been specified
    if (adopted_father_pt != 0)
    {
      // Get pointer to my father (in octree impersonation)
      father_pt = dynamic_cast<OcTree*>(adopted_father_pt);
    }
    // Check if element is in a tree
    else if (Tree_pt != 0)
    {
      // Get pointer to my father (in octree impersonation)
      father_pt = dynamic_cast<OcTree*>(octree_pt()->father_pt());
    }
    else
    {
      // throw OomphLibError(
      //       "Element not in a tree, and no adopted father has been
      //       specified!",
      //       "PRefineableQElement<2,INITIAL_NNODE_1D>::initial_setup()",
      //       OOMPH_EXCEPTION_LOCATION);
    }

    // Check if element has father
    if (father_pt != 0 || initial_p_order != 0)
    {
      if (father_pt != 0)
      {
        PRefineableQElement<3, INITIAL_NNODE_1D>* father_el_pt =
          dynamic_cast<PRefineableQElement<3, INITIAL_NNODE_1D>*>(
            father_pt->object_pt());
        if (father_el_pt != 0)
        {
          unsigned father_p_order = father_el_pt->p_order();
          // Set p-order to that of father
          P_order = father_p_order;
        }
      }
      else
      {
        P_order = initial_p_order;
      }

      // Now sort out the element...
      // (has p^3 nodes)
      unsigned new_n_node = P_order * P_order * P_order;

      // Allocate new space for Nodes (at the element level)
      this->set_n_node(new_n_node);

      // Set integration scheme
      delete this->integral_pt();
      switch (P_order)
      {
        case 2:
          this->set_integration_scheme(new GaussLobattoLegendre<3, 2>);
          break;
        case 3:
          this->set_integration_scheme(new GaussLobattoLegendre<3, 3>);
          break;
        case 4:
          this->set_integration_scheme(new GaussLobattoLegendre<3, 4>);
          break;
        case 5:
          this->set_integration_scheme(new GaussLobattoLegendre<3, 5>);
          break;
        case 6:
          this->set_integration_scheme(new GaussLobattoLegendre<3, 6>);
          break;
        case 7:
          this->set_integration_scheme(new GaussLobattoLegendre<3, 7>);
          break;
        default:
          std::ostringstream error_message;
          error_message << "\nERROR: Exceeded maximum polynomial order for";
          error_message << "\n       integration scheme.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }
  }

  //==================================================================
  /// Check the father element for any required nodes which
  /// already exist
  //==================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<3, INITIAL_NNODE_1D>::pre_build(
    Mesh*& mesh_pt, Vector<Node*>& new_node_pt)
  {
    using namespace OcTreeNames;

    // Get the number of 1d nodes
    unsigned n_p = nnode_1d();

    // Check whether static father_bound needs to be created
    if (Father_bound[n_p].nrow() == 0)
    {
      setup_father_bounds();
    }

    // Pointer to my father (in quadtree impersonation)
    OcTree* father_pt = dynamic_cast<OcTree*>(octree_pt()->father_pt());

    // What type of son am I? Ask my quadtree representation...
    int son_type = Tree_pt->son_type();

    // Has somebody build me already? (If any nodes have not been built)
    if (!nodes_built())
    {
#ifdef PARANOID
      if (father_pt == 0)
      {
        std::string error_message =
          "Something fishy here: I have no father and yet \n";
        error_message += "I have no nodes. Who has created me then?!\n";

        throw OomphLibError(
          error_message,
          "PRefineableQElement<3,INITIAL_NNODE_1D>::pre_build()",
          OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Return pointer to father element
      RefineableQElement<3>* father_el_pt =
        dynamic_cast<RefineableQElement<3>*>(father_pt->object_pt());

      // Timestepper should be the same for all nodes in father
      // element -- use it create timesteppers for new nodes
      TimeStepper* time_stepper_pt =
        father_el_pt->node_pt(0)->time_stepper_pt();

      // Number of history values (incl. present)
      unsigned ntstorage = time_stepper_pt->ntstorage();

      /// / Pass pointer to time object:
      // this->time_pt()=father_el_pt->time_pt();

      Vector<double> s_lo(3);
      Vector<double> s_hi(3);
      Vector<double> s(3);
      Vector<double> x(3);

      // Setup vertex coordinates in father element:
      //--------------------------------------------
      switch (son_type)
      {
        case LDB:
          s_lo[0] = -1.0;
          s_hi[0] = 0.0;
          s_lo[1] = -1.0;
          s_hi[1] = 0.0;
          s_lo[2] = -1.0;
          s_hi[2] = 0.0;
          break;

        case LDF:
          s_lo[0] = -1.0;
          s_hi[0] = 0.0;
          s_lo[1] = -1.0;
          s_hi[1] = 0.0;
          s_lo[2] = 0.0;
          s_hi[2] = 1.0;
          break;

        case LUB:
          s_lo[0] = -1.0;
          s_hi[0] = 0.0;
          s_lo[1] = 0.0;
          s_hi[1] = 1.0;
          s_lo[2] = -1.0;
          s_hi[2] = 0.0;
          break;

        case LUF:
          s_lo[0] = -1.0;
          s_hi[0] = 0.0;
          s_lo[1] = 0.0;
          s_hi[1] = 1.0;
          s_lo[2] = 0.0;
          s_hi[2] = 1.0;
          break;

        case RDB:
          s_lo[0] = 0.0;
          s_hi[0] = 1.0;
          s_lo[1] = -1.0;
          s_hi[1] = 0.0;
          s_lo[2] = -1.0;
          s_hi[2] = 0.0;
          break;

        case RDF:
          s_lo[0] = 0.0;
          s_hi[0] = 1.0;
          s_lo[1] = -1.0;
          s_hi[1] = 0.0;
          s_lo[2] = 0.0;
          s_hi[2] = 1.0;
          break;

        case RUB:
          s_lo[0] = 0.0;
          s_hi[0] = 1.0;
          s_lo[1] = 0.0;
          s_hi[1] = 1.0;
          s_lo[2] = -1.0;
          s_hi[2] = 0.0;
          break;

        case RUF:
          s_lo[0] = 0.0;
          s_hi[0] = 1.0;
          s_lo[1] = 0.0;
          s_hi[1] = 1.0;
          s_lo[2] = 0.0;
          s_hi[2] = 1.0;
          break;
      }

      /// / Pass macro element pointer on to sons and
      /// / set coordinates in macro element
      /// / hierher why can I see this?
      // if(father_el_pt->macro_elem_pt()!=0)
      // {
      //  set_macro_elem_pt(father_el_pt->macro_elem_pt());
      //  for(unsigned i=0;i<2;i++)
      //   {
      //    s_macro_ll(i)=      father_el_pt->s_macro_ll(i)+
      //     0.5*(s_lo[i]+1.0)*(father_el_pt->s_macro_ur(i)-
      //                        father_el_pt->s_macro_ll(i));
      //    s_macro_ur(i)=      father_el_pt->s_macro_ll(i)+
      //     0.5*(s_hi[i]+1.0)*(father_el_pt->s_macro_ur(i)-
      //                        father_el_pt->s_macro_ll(i));
      //   }
      // }


      // If the father element hasn't been generated yet, we're stuck...
      if (father_el_pt->node_pt(0) == 0)
      {
        throw OomphLibError(
          "Trouble: father_el_pt->node_pt(0)==0\n Can't build son element!\n",
          "PRefineableQElement<3,INITIAL_NNODE_1D>::pre_build()",
          OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        unsigned jnod = 0;

        Vector<double> s_fraction(3);
        // Loop over nodes in element
        for (unsigned i0 = 0; i0 < n_p; i0++)
        {
          // Get the fractional position of the node in the direction of s[0]
          s_fraction[0] = this->local_one_d_fraction_of_node(i0, 0);
          // Local coordinate in father element
          s[0] = s_lo[0] + (s_hi[0] - s_lo[0]) * s_fraction[0];

          for (unsigned i1 = 0; i1 < n_p; i1++)
          {
            // Get the fractional position of the node in the direction of s[1]
            s_fraction[1] = this->local_one_d_fraction_of_node(i1, 1);
            // Local coordinate in father element
            s[1] = s_lo[1] + (s_hi[1] - s_lo[1]) * s_fraction[1];

            for (unsigned i2 = 0; i2 < n_p; i2++)
            {
              // Get the fractional position of the node in the direction of
              // s[2]
              s_fraction[2] = this->local_one_d_fraction_of_node(i2, 2);
              // Local coordinate in father element
              s[2] = s_lo[2] + (s_hi[2] - s_lo[2]) * s_fraction[2];

              // Local node number
              jnod = i0 + n_p * i1 + n_p * n_p * i2;

              // Get the pointer to the node in the father, returns NULL
              // if there is not node
              Node* created_node_pt =
                father_el_pt->get_node_at_local_coordinate(s);

              // Does this node already exist in father element?
              //------------------------------------------------
              if (created_node_pt != 0)
              {
                // Copy node across
                this->node_pt(jnod) = created_node_pt;

                // Make sure that we update the values at the node so that
                // they are consistent with the present representation.
                // This is only need for mixed interpolation where the value
                // at the father could now become active.

                // Loop over all history values
                for (unsigned t = 0; t < ntstorage; t++)
                {
                  // Get values from father element
                  // Note: get_interpolated_values() sets Vector size itself.
                  Vector<double> prev_values;
                  father_el_pt->get_interpolated_values(t, s, prev_values);
                  // Find the minimum number of values
                  //(either those stored at the node, or those returned by
                  // the function)
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
              }
            } // End of loop i2
          } // End of loop i1
        } // End of loop i0
      } // Sanity check: Father element has been generated

    } // End of nodes not built
  }

  //==================================================================
  /// p-refine the element inc times. (If inc<0 then p-unrefinement
  /// is performed.)
  //==================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<3, INITIAL_NNODE_1D>::p_refine(
    const int& inc, Mesh* const& mesh_pt, GeneralisedElement* const& clone_pt)
  {
    // Number of dimensions
    unsigned n_dim = 3;

    // Cast clone to correct type
    PRefineableQElement<3, INITIAL_NNODE_1D>* clone_el_pt =
      dynamic_cast<PRefineableQElement<3, INITIAL_NNODE_1D>*>(clone_pt);

    // Check if we can use it
    if (clone_el_pt == 0)
    {
      throw OomphLibError(
        "Cloned copy must be a PRefineableQElement<3,INITIAL_NNODE_1D>!",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#ifdef PARANOID
    // Clone exists, so check that it is infact a copy of me
    else
    {
      // Flag to keep track of check
      bool clone_is_ok = true;

      // Does the clone have the correct p-order?
      clone_is_ok = clone_is_ok && (clone_el_pt->p_order() == this->p_order());

      if (!clone_is_ok)
      {
        std::ostringstream err_stream;
        err_stream << "Clone element has a different p-order from me,"
                   << std::endl
                   << "but it is supposed to be a copy!" << std::endl;

        throw OomphLibError(
          err_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Does the clone have the same number of nodes as me?
      clone_is_ok = clone_is_ok && (clone_el_pt->nnode() == this->nnode());

      if (!clone_is_ok)
      {
        std::ostringstream err_stream;
        err_stream << "Clone element has a different number of nodes from me,"
                   << std::endl
                   << "but it is supposed to be a copy!" << std::endl;

        throw OomphLibError(
          err_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Does the clone have the same nodes as me?
      for (unsigned n = 0; n < this->nnode(); n++)
      {
        clone_is_ok =
          clone_is_ok && (clone_el_pt->node_pt(n) == this->node_pt(n));
      }

      if (!clone_is_ok)
      {
        std::ostringstream err_stream;
        err_stream << "The nodes belonging to the clone element are different"
                   << std::endl
                   << "from mine, but it is supposed to be a copy!"
                   << std::endl;

        throw OomphLibError(
          err_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // If we get to here then the clone has all the information we require
    }
#endif

    // Currently we can't handle the case of generalised coordinates
    // since we haven't established how they should be interpolated.
    // Buffer this case:
    if (clone_el_pt->node_pt(0)->nposition_type() != 1)
    {
      throw OomphLibError("Can't handle generalised nodal positions (yet).",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Timestepper should be the same for all nodes -- use it
    // to create timesteppers for new nodes
    TimeStepper* time_stepper_pt = this->node_pt(0)->time_stepper_pt();

    // Get number of history values (incl. present)
    unsigned ntstorage = time_stepper_pt->ntstorage();

    // Increment p-order of the element
    this->P_order += inc;

    // Change integration scheme
    delete this->integral_pt();
    switch (this->P_order)
    {
      case 2:
        this->set_integration_scheme(new GaussLobattoLegendre<3, 2>);
        break;
      case 3:
        this->set_integration_scheme(new GaussLobattoLegendre<3, 3>);
        break;
      case 4:
        this->set_integration_scheme(new GaussLobattoLegendre<3, 4>);
        break;
      case 5:
        this->set_integration_scheme(new GaussLobattoLegendre<3, 5>);
        break;
      case 6:
        this->set_integration_scheme(new GaussLobattoLegendre<3, 6>);
        break;
      case 7:
        this->set_integration_scheme(new GaussLobattoLegendre<3, 7>);
        break;
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       integration scheme.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }

    // Allocate new space for Nodes (at the element level)
    this->set_n_node(this->P_order * this->P_order * this->P_order);

    // Copy vertex nodes and create new edge and internal nodes
    //---------------------------------------------------------

    // Setup vertex coordinates in element:
    //-------------------------------------
    Vector<double> s_lo(n_dim, 0.0);
    Vector<double> s_hi(n_dim, 0.0);
    s_lo[0] = -1.0;
    s_hi[0] = 1.0;
    s_lo[1] = -1.0;
    s_hi[1] = 1.0;
    s_lo[2] = -1.0;
    s_hi[2] = 1.0;

    // Local coordinate in element
    Vector<double> s(n_dim, 0.0);

    unsigned jnod = 0;

    Vector<double> s_fraction(n_dim, 0.0);
    // Loop over nodes in element
    for (unsigned i0 = 0; i0 < this->P_order; i0++)
    {
      // Get the fractional position of the node in the direction of s[0]
      s_fraction[0] = this->local_one_d_fraction_of_node(i0, 0);
      // Local coordinate
      s[0] = s_lo[0] + (s_hi[0] - s_lo[0]) * s_fraction[0];

      for (unsigned i1 = 0; i1 < this->P_order; i1++)
      {
        // Get the fractional position of the node in the direction of s[1]
        s_fraction[1] = this->local_one_d_fraction_of_node(i1, 1);
        // Local coordinate
        s[1] = s_lo[1] + (s_hi[1] - s_lo[1]) * s_fraction[1];

        for (unsigned i2 = 0; i2 < this->P_order; i2++)
        {
          // Get the fractional position of the node in the direction of s[2]
          s_fraction[2] = this->local_one_d_fraction_of_node(i2, 2);
          // Local coordinate
          s[2] = s_lo[2] + (s_hi[2] - s_lo[2]) * s_fraction[2];

          // Local node number
          jnod = i0 + this->P_order * i1 + this->P_order * this->P_order * i2;

          // Initialise flag: So far, this node hasn't been built
          // or copied yet
          bool node_done = false;

          // Get the pointer to the node in this element (or rather, its clone),
          // returns NULL if there is not node
          Node* created_node_pt = clone_el_pt->get_node_at_local_coordinate(s);
          // created_node_pt = 0;

          // Does this node already exist in this element?
          //----------------------------------------------
          if (created_node_pt != 0)
          {
            // Copy node across
            this->node_pt(jnod) = created_node_pt;

            // Make sure that we update the values at the node so that
            // they are consistent with the present representation.
            // This is only need for mixed interpolation where the value
            // at the father could now become active.

            // Loop over all history values
            for (unsigned t = 0; t < ntstorage; t++)
            {
              // Get values from father element
              // Note: get_interpolated_values() sets Vector size itself.
              Vector<double> prev_values;
              clone_el_pt->get_interpolated_values(t, s, prev_values);
              // Find the minimum number of values
              //(either those stored at the node, or those returned by
              // the function)
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

            // Node has been created by copy
            node_done = true;
          }
          // Node does not exist in this element but might already
          //------------------------------------------------------
          // have been created by neighbouring elements
          //-------------------------------------------
          else
          {
            // Was the node created by one of its neighbours
            // Whether or not the node lies on an edge can be calculated
            // by from the fractional position
            bool is_periodic = false;
            created_node_pt =
              this->node_created_by_neighbour(s_fraction, is_periodic);

            // If the node was so created, assign the pointers
            if (created_node_pt != 0)
            {
              // If the node is periodic
              if (is_periodic)
              {
                // Now the node must be on a boundary, but we don't know which
                // one
                // The returned created_node_pt is actually the neighbouring
                // periodic node
                Node* neighbour_node_pt = created_node_pt;

                // Determine the edge on which the new node will live
                //(cannot be a vertex node)
                int my_bound = Tree::OMEGA;
                // If we are on the left face
                if (i0 == 0)
                {
                  my_bound = OcTreeNames::L;
                }
                // If we are on the right face
                else if (i0 == this->P_order - 1)
                {
                  my_bound = OcTreeNames::R;
                }

                // If we are on the bottom face
                if (i1 == 0)
                {
                  // If we already have already set the boundary, we're on an
                  // edge
                  switch (my_bound)
                  {
                    case OcTreeNames::L:
                      my_bound = OcTreeNames::LD;
                      break;
                    case OcTreeNames::R:
                      my_bound = OcTreeNames::RD;
                      break;
                    // Boundary not set
                    default:
                      my_bound = OcTreeNames::D;
                      break;
                  }
                }
                // If we are the top face
                else if (i1 == this->P_order - 1)
                {
                  // If we already have a boundary
                  switch (my_bound)
                  {
                    case OcTreeNames::L:
                      my_bound = OcTreeNames::LU;
                      break;
                    case OcTreeNames::R:
                      my_bound = OcTreeNames::RU;
                      break;
                    default:
                      my_bound = OcTreeNames::U;
                      break;
                  }
                }

                // If we are on the back face
                if (i2 == 0)
                {
                  // If we already have already set the boundary, we're on an
                  // edge
                  switch (my_bound)
                  {
                    case OcTreeNames::L:
                      my_bound = OcTreeNames::LB;
                      break;
                    case OcTreeNames::LD:
                      my_bound = OcTreeNames::LDB;
                      break;
                    case OcTreeNames::LU:
                      my_bound = OcTreeNames::LUB;
                      break;
                    case OcTreeNames::R:
                      my_bound = OcTreeNames::RB;
                      break;
                    case OcTreeNames::RD:
                      my_bound = OcTreeNames::RDB;
                      break;
                    case OcTreeNames::RU:
                      my_bound = OcTreeNames::RUB;
                      break;
                    case OcTreeNames::D:
                      my_bound = OcTreeNames::DB;
                      break;
                    case OcTreeNames::U:
                      my_bound = OcTreeNames::UB;
                      break;
                    // Boundary not set
                    default:
                      my_bound = OcTreeNames::B;
                      break;
                  }
                }
                // If we are the front face
                else if (i2 == this->P_order - 1)
                {
                  // If we already have a boundary
                  switch (my_bound)
                  {
                    case OcTreeNames::L:
                      my_bound = OcTreeNames::LF;
                      break;
                    case OcTreeNames::LD:
                      my_bound = OcTreeNames::LDF;
                      break;
                    case OcTreeNames::LU:
                      my_bound = OcTreeNames::LUF;
                      break;
                    case OcTreeNames::R:
                      my_bound = OcTreeNames::RF;
                      break;
                    case OcTreeNames::RD:
                      my_bound = OcTreeNames::RDF;
                      break;
                    case OcTreeNames::RU:
                      my_bound = OcTreeNames::RUF;
                      break;
                    case OcTreeNames::D:
                      my_bound = OcTreeNames::DF;
                      break;
                    case OcTreeNames::U:
                      my_bound = OcTreeNames::UF;
                      break;
                    default:
                      my_bound = OcTreeNames::F;
                      break;
                  }
                }

                // Storage for the set of Mesh boundaries on which the
                // appropriate edge lives.
                // [New nodes should always be mid-edge nodes and therefore
                // only live on one boundary but just to play it safe...]
                std::set<unsigned> boundaries;
                // Only get the boundaries if we are at the edge of
                // an element. Nodes in the centre of an element cannot be
                // on Mesh boundaries
                if (my_bound != Tree::OMEGA)
                {
                  clone_el_pt->get_boundaries(my_bound, boundaries);
                }

#ifdef PARANOID
                // Case where a new node lives on more than one boundary
                // seems fishy enough to flag
                if (boundaries.size() > 2)
                {
                  throw OomphLibError(
                    "boundaries.size()>2 seems a bit strange..\n",
                    OOMPH_CURRENT_FUNCTION,
                    OOMPH_EXCEPTION_LOCATION);
                }

                // Case when there are no boundaries, we are in big trouble
                if (boundaries.size() == 0)
                {
                  std::ostringstream error_stream;
                  error_stream << "Periodic node is not on a boundary...\n"
                               << "Coordinates: " << created_node_pt->x(0)
                               << " " << created_node_pt->x(1) << "\n";
                  throw OomphLibError(error_stream.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
#endif

                // Create node and set the pointer to it from the element
                created_node_pt =
                  this->construct_boundary_node(jnod, time_stepper_pt);
                // Make the node periodic from the neighbour
                created_node_pt->make_periodic(neighbour_node_pt);

                // Loop over # of history values
                for (unsigned t = 0; t < ntstorage; t++)
                {
                  // Get position from father element -- this uses the macro
                  // element representation if appropriate. If the node
                  // turns out to be a hanging node later on, then
                  // its position gets adjusted in line with its
                  // hanging node interpolation.
                  Vector<double> x_prev(n_dim, 0.0);
                  clone_el_pt->get_x(t, s, x_prev);
                  // Set previous positions of the new node
                  for (unsigned i = 0; i < n_dim; i++)
                  {
                    created_node_pt->x(t, i) = x_prev[i];
                  }
                }

                // Check if we need to add nodes to the mesh
                if (mesh_pt != 0)
                {
                  // Next, we Update the boundary lookup schemes
                  // Loop over the boundaries stored in the set
                  for (std::set<unsigned>::iterator it = boundaries.begin();
                       it != boundaries.end();
                       ++it)
                  {
                    // Add the node to the boundary
                    mesh_pt->add_boundary_node(*it, created_node_pt);

                    // If we have set an intrinsic coordinate on this
                    // mesh boundary then it must also be interpolated on
                    // the new node
                    // Now interpolate the intrinsic boundary coordinate
                    if (mesh_pt->boundary_coordinate_exists(*it) == true)
                    {
                      Vector<double> zeta(2, 0.0);
                      clone_el_pt->interpolated_zeta_on_face(
                        *it, my_bound, s, zeta);

                      created_node_pt->set_coordinates_on_boundary(*it, zeta);
                    }
                  }

                  // Make sure that we add the node to the mesh
                  mesh_pt->add_node_pt(created_node_pt);
                }
              } // End of periodic case
              // Otherwise the node is not periodic, so just set the
              // pointer to the neighbours node
              else
              {
                this->node_pt(jnod) = created_node_pt;
              }
              // Node has been created
              node_done = true;
            }
            // Node does not exist in neighbour element but might already
            //-----------------------------------------------------------
            // have been created by a son of a neighbouring element
            //-----------------------------------------------------
            else
            {
              // Was the node created by one of its neighbours' sons
              // Whether or not the node lies on an edge can be calculated
              // by from the fractional position
              bool is_periodic = false;
              created_node_pt =
                this->node_created_by_son_of_neighbour(s_fraction, is_periodic);

              // If the node was so created, assign the pointers
              if (created_node_pt != 0)
              {
                // If the node is periodic
                if (is_periodic)
                {
                  // Now the node must be on a boundary, but we don't know which
                  // one
                  // The returned created_node_pt is actually the neighbouring
                  // periodic node
                  Node* neighbour_node_pt = created_node_pt;

                  // Determine the edge on which the new node will live
                  //(cannot be a vertex node)
                  int my_bound = Tree::OMEGA;
                  // If we are on the left face
                  if (i0 == 0)
                  {
                    my_bound = OcTreeNames::L;
                  }
                  // If we are on the right face
                  else if (i0 == this->P_order - 1)
                  {
                    my_bound = OcTreeNames::R;
                  }

                  // If we are on the bottom face
                  if (i1 == 0)
                  {
                    // If we already have already set the boundary, we're on an
                    // edge
                    switch (my_bound)
                    {
                      case OcTreeNames::L:
                        my_bound = OcTreeNames::LD;
                        break;
                      case OcTreeNames::R:
                        my_bound = OcTreeNames::RD;
                        break;
                      // Boundary not set
                      default:
                        my_bound = OcTreeNames::D;
                        break;
                    }
                  }
                  // If we are the top face
                  else if (i1 == this->P_order - 1)
                  {
                    // If we already have a boundary
                    switch (my_bound)
                    {
                      case OcTreeNames::L:
                        my_bound = OcTreeNames::LU;
                        break;
                      case OcTreeNames::R:
                        my_bound = OcTreeNames::RU;
                        break;
                      default:
                        my_bound = OcTreeNames::U;
                        break;
                    }
                  }

                  // If we are on the back face
                  if (i2 == 0)
                  {
                    // If we already have already set the boundary, we're on an
                    // edge
                    switch (my_bound)
                    {
                      case OcTreeNames::L:
                        my_bound = OcTreeNames::LB;
                        break;
                      case OcTreeNames::LD:
                        my_bound = OcTreeNames::LDB;
                        break;
                      case OcTreeNames::LU:
                        my_bound = OcTreeNames::LUB;
                        break;
                      case OcTreeNames::R:
                        my_bound = OcTreeNames::RB;
                        break;
                      case OcTreeNames::RD:
                        my_bound = OcTreeNames::RDB;
                        break;
                      case OcTreeNames::RU:
                        my_bound = OcTreeNames::RUB;
                        break;
                      case OcTreeNames::D:
                        my_bound = OcTreeNames::DB;
                        break;
                      case OcTreeNames::U:
                        my_bound = OcTreeNames::UB;
                        break;
                      // Boundary not set
                      default:
                        my_bound = OcTreeNames::B;
                        break;
                    }
                  }
                  // If we are the front face
                  else if (i2 == this->P_order - 1)
                  {
                    // If we already have a boundary
                    switch (my_bound)
                    {
                      case OcTreeNames::L:
                        my_bound = OcTreeNames::LF;
                        break;
                      case OcTreeNames::LD:
                        my_bound = OcTreeNames::LDF;
                        break;
                      case OcTreeNames::LU:
                        my_bound = OcTreeNames::LUF;
                        break;
                      case OcTreeNames::R:
                        my_bound = OcTreeNames::RF;
                        break;
                      case OcTreeNames::RD:
                        my_bound = OcTreeNames::RDF;
                        break;
                      case OcTreeNames::RU:
                        my_bound = OcTreeNames::RUF;
                        break;
                      case OcTreeNames::D:
                        my_bound = OcTreeNames::DF;
                        break;
                      case OcTreeNames::U:
                        my_bound = OcTreeNames::UF;
                        break;
                      default:
                        my_bound = OcTreeNames::F;
                        break;
                    }
                  }

                  // Storage for the set of Mesh boundaries on which the
                  // appropriate edge lives.
                  // [New nodes should always be mid-edge nodes and therefore
                  // only live on one boundary but just to play it safe...]
                  std::set<unsigned> boundaries;
                  // Only get the boundaries if we are at the edge of
                  // an element. Nodes in the centre of an element cannot be
                  // on Mesh boundaries
                  if (my_bound != Tree::OMEGA)
                  {
                    clone_el_pt->get_boundaries(my_bound, boundaries);
                  }

#ifdef PARANOID
                  // Case where a new node lives on more than one boundary
                  // seems fishy enough to flag
                  if (boundaries.size() > 2)
                  {
                    throw OomphLibError(
                      "boundaries.size()>2 seems a bit strange..\n",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
                  }

                  // Case when there are no boundaries, we are in big trouble
                  if (boundaries.size() == 0)
                  {
                    std::ostringstream error_stream;
                    error_stream << "Periodic node is not on a boundary...\n"
                                 << "Coordinates: " << created_node_pt->x(0)
                                 << " " << created_node_pt->x(1) << "\n";
                    throw OomphLibError(error_stream.str(),
                                        OOMPH_CURRENT_FUNCTION,
                                        OOMPH_EXCEPTION_LOCATION);
                  }
#endif

                  // Create node and set the pointer to it from the element
                  created_node_pt =
                    this->construct_boundary_node(jnod, time_stepper_pt);
                  // Make the node periodic from the neighbour
                  created_node_pt->make_periodic(neighbour_node_pt);

                  // Loop over # of history values
                  for (unsigned t = 0; t < ntstorage; t++)
                  {
                    // Get position from father element -- this uses the macro
                    // element representation if appropriate. If the node
                    // turns out to be a hanging node later on, then
                    // its position gets adjusted in line with its
                    // hanging node interpolation.
                    Vector<double> x_prev(n_dim);
                    clone_el_pt->get_x(t, s, x_prev);
                    // Set previous positions of the new node
                    for (unsigned i = 0; i < n_dim; i++)
                    {
                      created_node_pt->x(t, i) = x_prev[i];
                    }
                  }

                  // Check if we need to add nodes to the mesh
                  if (mesh_pt != 0)
                  {
                    // Next, we Update the boundary lookup schemes
                    // Loop over the boundaries stored in the set
                    for (std::set<unsigned>::iterator it = boundaries.begin();
                         it != boundaries.end();
                         ++it)
                    {
                      // Add the node to the boundary
                      mesh_pt->add_boundary_node(*it, created_node_pt);

                      // If we have set an intrinsic coordinate on this
                      // mesh boundary then it must also be interpolated on
                      // the new node
                      // Now interpolate the intrinsic boundary coordinate
                      if (mesh_pt->boundary_coordinate_exists(*it) == true)
                      {
                        Vector<double> zeta(2, 0.0);
                        clone_el_pt->interpolated_zeta_on_face(
                          *it, my_bound, s, zeta);

                        created_node_pt->set_coordinates_on_boundary(*it, zeta);
                      }
                    }

                    // Make sure that we add the node to the mesh
                    mesh_pt->add_node_pt(created_node_pt);
                  }
                } // End of periodic case
                // Otherwise the node is not periodic, so just set the
                // pointer to the neighbours node
                else
                {
                  this->node_pt(jnod) = created_node_pt;
                }
                // Node has been created
                node_done = true;
              } // Node does not exist in son of neighbouring element
            } // Node does not exist in neighbouring element
          } // Node does not exist in this element

          // Node has not been built anywhere ---> build it here
          if (!node_done)
          {
            // Firstly, we need to determine whether or not a node lies
            // on the boundary before building it, because
            // we actually assign a different type of node on boundaries.

            // Initially none
            int my_bound = Tree::OMEGA;
            // If we are on the left face
            if (i0 == 0)
            {
              my_bound = OcTreeNames::L;
            }
            // If we are on the right face
            else if (i0 == this->P_order - 1)
            {
              my_bound = OcTreeNames::R;
            }

            // If we are on the bottom face
            if (i1 == 0)
            {
              // If we already have already set the boundary, we're on an edge
              switch (my_bound)
              {
                case OcTreeNames::L:
                  my_bound = OcTreeNames::LD;
                  break;
                case OcTreeNames::R:
                  my_bound = OcTreeNames::RD;
                  break;
                // Boundary not set
                default:
                  my_bound = OcTreeNames::D;
                  break;
              }
            }
            // If we are the top face
            else if (i1 == this->P_order - 1)
            {
              // If we already have a boundary
              switch (my_bound)
              {
                case OcTreeNames::L:
                  my_bound = OcTreeNames::LU;
                  break;
                case OcTreeNames::R:
                  my_bound = OcTreeNames::RU;
                  break;
                default:
                  my_bound = OcTreeNames::U;
                  break;
              }
            }

            // If we are on the back face
            if (i2 == 0)
            {
              // If we already have already set the boundary, we're on an edge
              switch (my_bound)
              {
                case OcTreeNames::L:
                  my_bound = OcTreeNames::LB;
                  break;
                case OcTreeNames::LD:
                  my_bound = OcTreeNames::LDB;
                  break;
                case OcTreeNames::LU:
                  my_bound = OcTreeNames::LUB;
                  break;
                case OcTreeNames::R:
                  my_bound = OcTreeNames::RB;
                  break;
                case OcTreeNames::RD:
                  my_bound = OcTreeNames::RDB;
                  break;
                case OcTreeNames::RU:
                  my_bound = OcTreeNames::RUB;
                  break;
                case OcTreeNames::D:
                  my_bound = OcTreeNames::DB;
                  break;
                case OcTreeNames::U:
                  my_bound = OcTreeNames::UB;
                  break;
                // Boundary not set
                default:
                  my_bound = OcTreeNames::B;
                  break;
              }
            }
            // If we are the front face
            else if (i2 == this->P_order - 1)
            {
              // If we already have a boundary
              switch (my_bound)
              {
                case OcTreeNames::L:
                  my_bound = OcTreeNames::LF;
                  break;
                case OcTreeNames::LD:
                  my_bound = OcTreeNames::LDF;
                  break;
                case OcTreeNames::LU:
                  my_bound = OcTreeNames::LUF;
                  break;
                case OcTreeNames::R:
                  my_bound = OcTreeNames::RF;
                  break;
                case OcTreeNames::RD:
                  my_bound = OcTreeNames::RDF;
                  break;
                case OcTreeNames::RU:
                  my_bound = OcTreeNames::RUF;
                  break;
                case OcTreeNames::D:
                  my_bound = OcTreeNames::DF;
                  break;
                case OcTreeNames::U:
                  my_bound = OcTreeNames::UF;
                  break;
                default:
                  my_bound = OcTreeNames::F;
                  break;
              }
            }

            // Storage for the set of Mesh boundaries on which the
            // appropriate edge lives.
            // [New nodes should always be mid-edge nodes and therefore
            // only live on one boundary but just to play it safe...]
            std::set<unsigned> boundaries;
            // Only get the boundaries if we are at the edge of
            // an element. Nodes in the centre of an element cannot be
            // on Mesh boundaries
            if (my_bound != Tree::OMEGA)
              clone_el_pt->get_boundaries(my_bound, boundaries);

#ifdef PARANOID
            // Case where a new node lives on more than two boundaries
            // seems fishy enough to flag
            if (boundaries.size() > 2)
            {
              throw OomphLibError(
                "boundaries.size()>2 seems a bit strange..\n",
                "PRefineableQElement<3,INITIAL_NNODE_1D>::p_refine()",
                OOMPH_EXCEPTION_LOCATION);
            }
#endif

            // If the node lives on a mesh boundary,
            // then we need to create a boundary node
            if (boundaries.size() > 0)
            {
              // Create node and set the pointer to it from the element
              created_node_pt =
                this->construct_boundary_node(jnod, time_stepper_pt);

              // Now we need to work out whether to pin the values at
              // the new node based on the boundary conditions applied at
              // its Mesh boundary

              // Get the boundary conditions from the father
              Vector<int> bound_cons(this->ncont_interpolated_values());
              clone_el_pt->get_bcs(my_bound, bound_cons);

              // Loop over the values and pin, if necessary
              unsigned n_value = created_node_pt->nvalue();
              for (unsigned k = 0; k < n_value; k++)
              {
                if (bound_cons[k])
                {
                  created_node_pt->pin(k);
                }
              }

              // Solid node? If so, deal with the positional boundary
              // conditions:
              SolidNode* solid_node_pt =
                dynamic_cast<SolidNode*>(created_node_pt);
              if (solid_node_pt != 0)
              {
                // Get the positional boundary conditions from the father:
                unsigned n_dim = created_node_pt->ndim();
                Vector<int> solid_bound_cons(n_dim);
                RefineableSolidQElement<3>* clone_solid_el_pt =
                  dynamic_cast<RefineableSolidQElement<3>*>(clone_el_pt);
#ifdef PARANOID
                if (clone_solid_el_pt == 0)
                {
                  std::string error_message =
                    "We have a SolidNode outside a refineable SolidElement\n";
                  error_message +=
                    "during mesh refinement -- this doesn't make sense";

                  throw OomphLibError(
                    error_message,
                    "PRefineableQElement<3,INITIAL_NNODE_1D>::p_refine()",
                    OOMPH_EXCEPTION_LOCATION);
                }
#endif
                clone_solid_el_pt->get_solid_bcs(my_bound, solid_bound_cons);

                // Loop over the positions and pin, if necessary
                for (unsigned k = 0; k < n_dim; k++)
                {
                  if (solid_bound_cons[k])
                  {
                    solid_node_pt->pin_position(k);
                  }
                }
              } // End of if solid_node_pt

              // Next, we Update the boundary lookup schemes

              // Check if we need to add nodes to the mesh
              if (mesh_pt != 0)
              {
                // Loop over the boundaries stored in the set
                for (std::set<unsigned>::iterator it = boundaries.begin();
                     it != boundaries.end();
                     ++it)
                {
                  // Add the node to the boundary
                  mesh_pt->add_boundary_node(*it, created_node_pt);

                  // If we have set an intrinsic coordinate on this
                  // mesh boundary then it must also be interpolated on
                  // the new node
                  // Now interpolate the intrinsic boundary coordinate
                  if (mesh_pt->boundary_coordinate_exists(*it) == true)
                  {
                    Vector<double> zeta(2);
                    clone_el_pt->interpolated_zeta_on_face(
                      *it, my_bound, s, zeta);

                    created_node_pt->set_coordinates_on_boundary(*it, zeta);
                  }
                }
              }
            }
            // Otherwise the node is not on a Mesh boundary and
            // we create a normal "bulk" node
            else
            {
              // Create node and set the pointer to it from the element
              created_node_pt = this->construct_node(jnod, time_stepper_pt);
            }

            // Now we set the position and values at the newly created node

            // In the first instance use macro element or FE representation
            // to create past and present nodal positions.
            // (THIS STEP SHOULD NOT BE SKIPPED FOR ALGEBRAIC
            // ELEMENTS AS NOT ALL OF THEM NECESSARILY IMPLEMENT
            // NONTRIVIAL NODE UPDATE FUNCTIONS. CALLING
            // THE NODE UPDATE FOR SUCH ELEMENTS/NODES WILL LEAVE
            // THEIR NODAL POSITIONS WHERE THEY WERE (THIS IS APPROPRIATE
            // ONCE THEY HAVE BEEN GIVEN POSITIONS) BUT WILL
            // NOT ASSIGN SENSIBLE INITIAL POSITONS!

            // Loop over # of history values
            for (unsigned t = 0; t < ntstorage; t++)
            {
              // Get position from father element -- this uses the macro
              // element representation if appropriate. If the node
              // turns out to be a hanging node later on, then
              // its position gets adjusted in line with its
              // hanging node interpolation.
              Vector<double> x_prev(3);
              clone_el_pt->get_x(t, s, x_prev);

              // Set previous positions of the new node
              for (unsigned i = 0; i < 3; i++)
              {
                created_node_pt->x(t, i) = x_prev[i];
              }
            }

            // Loop over all history values
            for (unsigned t = 0; t < ntstorage; t++)
            {
              // Get values from father element
              // Note: get_interpolated_values() sets Vector size itself.
              Vector<double> prev_values;
              clone_el_pt->get_interpolated_values(t, s, prev_values);

              // Initialise the values at the new node
              unsigned n_value = created_node_pt->nvalue();
              for (unsigned k = 0; k < n_value; k++)
              {
                created_node_pt->set_value(t, k, prev_values[k]);
              }
            }

            // Add new node to mesh (if requested)
            if (mesh_pt != 0)
            {
              mesh_pt->add_node_pt(created_node_pt);
            }

            AlgebraicElementBase* alg_el_pt =
              dynamic_cast<AlgebraicElementBase*>(this);

            // If we do have an algebraic element
            if (alg_el_pt != 0)
            {
              std::string error_message =
                "Have not implemented p-refinement for";
              error_message += "Algebraic p-refineable elements yet\n";

              throw OomphLibError(
                error_message,
                "PRefineableQElement<3,INITIAL_NNODE_1D>::p_refine()",
                OOMPH_EXCEPTION_LOCATION);
            }

          } // End of case when we build the node ourselves

          // Check if the element is an algebraic element
          AlgebraicElementBase* alg_el_pt =
            dynamic_cast<AlgebraicElementBase*>(this);

          // If the element is an algebraic element, setup
          // node position (past and present) from algebraic node update
          // function. This over-writes previous assingments that
          // were made based on the macro-element/FE representation.
          // NOTE: YES, THIS NEEDS TO BE CALLED REPEATEDLY IF THE
          // NODE IS MEMBER OF MULTIPLE ELEMENTS: THEY ALL ASSIGN
          // THE SAME NODAL POSITIONS BUT WE NEED TO ADD THE REMESH
          // INFO FOR *ALL* ROOT ELEMENTS!
          if (alg_el_pt != 0)
          {
            // Build algebraic node update info for new node
            // This sets up the node update data for all node update
            // functions that are shared by all nodes in the father
            // element
            alg_el_pt->setup_algebraic_node_update(
              this->node_pt(jnod), s, clone_el_pt);
          }

        } // End of vertical loop over nodes in element (i2)

      } // End of horizontal loop over nodes in element (i1)

    } // End of horizontal loop over nodes in element (i0)


    // Loop over all nodes in element again, to re-set the positions
    // This must be done using the new element's macro-element
    // representation, rather than the old version which may be
    // of a different p-order!
    for (unsigned i0 = 0; i0 < this->P_order; i0++)
    {
      // Get the fractional position of the node in the direction of s[0]
      s_fraction[0] = this->local_one_d_fraction_of_node(i0, 0);
      // Local coordinate
      s[0] = s_lo[0] + (s_hi[0] - s_lo[0]) * s_fraction[0];

      for (unsigned i1 = 0; i1 < this->P_order; i1++)
      {
        // Get the fractional position of the node in the direction of s[1]
        s_fraction[1] = this->local_one_d_fraction_of_node(i1, 1);
        // Local coordinate
        s[1] = s_lo[1] + (s_hi[1] - s_lo[1]) * s_fraction[1];

        for (unsigned i2 = 0; i2 < this->P_order; i2++)
        {
          // Get the fractional position of the node in the direction of s[2]
          s_fraction[2] = this->local_one_d_fraction_of_node(i2, 2);
          // Local coordinate
          s[2] = s_lo[2] + (s_hi[2] - s_lo[2]) * s_fraction[2];

          // Local node number
          jnod = i0 + this->P_order * i1 + this->P_order * this->P_order * i2;

          // Loop over # of history values
          for (unsigned t = 0; t < ntstorage; t++)
          {
            // Get position from father element -- this uses the macro
            // element representation if appropriate. If the node
            // turns out to be a hanging node later on, then
            // its position gets adjusted in line with its
            // hanging node interpolation.
            Vector<double> x_prev(3);
            this->get_x(t, s, x_prev);

            // Set previous positions of the new node
            for (unsigned i = 0; i < 3; i++)
            {
              this->node_pt(jnod)->x(t, i) = x_prev[i];
            }
          }
        }
      }
    }


    // If the element is a MacroElementNodeUpdateElement, set
    // the update parameters for the current element's nodes --
    // all this needs is the vector of (pointers to the)
    // geometric objects that affect the MacroElement-based
    // node update -- this needs to be done to set the node
    // update info for newly created nodes
    MacroElementNodeUpdateElementBase* clone_m_el_pt =
      dynamic_cast<MacroElementNodeUpdateElementBase*>(clone_el_pt);
    if (clone_m_el_pt != 0)
    {
      // Get vector of geometric objects from father (construct vector
      // via copy operation)
      Vector<GeomObject*> geom_object_pt(clone_m_el_pt->geom_object_pt());

      // Cast current element to MacroElementNodeUpdateElement:
      MacroElementNodeUpdateElementBase* m_el_pt =
        dynamic_cast<MacroElementNodeUpdateElementBase*>(this);

#ifdef PARANOID
      if (m_el_pt == 0)
      {
        std::string error_message =
          "Failed to cast to MacroElementNodeUpdateElementBase*\n";
        error_message +=
          "Strange -- if my clone is a MacroElementNodeUpdateElement\n";
        error_message += "then I should be too....\n";

        throw OomphLibError(
          error_message,
          "PRefineableQElement<3,INITIAL_NNODE_1D>::p_refine()",
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // Build update info by passing vector of geometric objects:
      // This sets the current element to be the update element
      // for all of the element's nodes -- this is reversed
      // if the element is ever un-refined in the father element's
      // rebuild_from_sons() function which overwrites this
      // assignment to avoid nasty segmentation faults that occur
      // when a node tries to update itself via an element that no
      // longer exists...
      m_el_pt->set_node_update_info(geom_object_pt);
    }

    // Not necessary to delete the old nodes since all original nodes are in the
    // current mesh and so will be pruned as part of the mesh adaption process.

    // Do any further-build required
    this->further_build();
  }

  //=======================================================================
  /// Shape functions for PRefineableQElement<DIM>
  //=======================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<3, INITIAL_NNODE_1D>::shape(const Vector<double>& s,
                                                       Shape& psi) const
  {
    switch (P_order)
    {
      case 2:
      {
        // Call the OneDimensional Shape functions
        OneDimensionalLegendreShape<2>::calculate_nodal_positions();
        OneDimensionalLegendreShape<2> psi1(s[0]), psi2(s[1]), psi3(s[2]);

        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < P_order; i++)
        {
          for (unsigned j = 0; j < P_order; j++)
          {
            for (unsigned k = 0; k < P_order; k++)
            {
              psi(P_order * P_order * i + P_order * j + k) =
                psi3[i] * psi2[j] * psi1[k];
            }
          }
        }
        break;
      }
      case 3:
      {
        // Call the OneDimensional Shape functions
        OneDimensionalLegendreShape<3>::calculate_nodal_positions();
        OneDimensionalLegendreShape<3> psi1(s[0]), psi2(s[1]), psi3(s[2]);

        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < P_order; i++)
        {
          for (unsigned j = 0; j < P_order; j++)
          {
            for (unsigned k = 0; k < P_order; k++)
            {
              psi(P_order * P_order * i + P_order * j + k) =
                psi3[i] * psi2[j] * psi1[k];
            }
          }
        }
        break;
      }
      case 4:
      {
        // Call the OneDimensional Shape functions
        OneDimensionalLegendreShape<4>::calculate_nodal_positions();
        OneDimensionalLegendreShape<4> psi1(s[0]), psi2(s[1]), psi3(s[2]);

        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < P_order; i++)
        {
          for (unsigned j = 0; j < P_order; j++)
          {
            for (unsigned k = 0; k < P_order; k++)
            {
              psi(P_order * P_order * i + P_order * j + k) =
                psi3[i] * psi2[j] * psi1[k];
            }
          }
        }
        break;
      }
      case 5:
      {
        // Call the OneDimensional Shape functions
        OneDimensionalLegendreShape<5>::calculate_nodal_positions();
        OneDimensionalLegendreShape<5> psi1(s[0]), psi2(s[1]), psi3(s[2]);

        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < P_order; i++)
        {
          for (unsigned j = 0; j < P_order; j++)
          {
            for (unsigned k = 0; k < P_order; k++)
            {
              psi(P_order * P_order * i + P_order * j + k) =
                psi3[i] * psi2[j] * psi1[k];
            }
          }
        }
        break;
      }
      case 6:
      {
        // Call the OneDimensional Shape functions
        OneDimensionalLegendreShape<6>::calculate_nodal_positions();
        OneDimensionalLegendreShape<6> psi1(s[0]), psi2(s[1]), psi3(s[2]);

        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < P_order; i++)
        {
          for (unsigned j = 0; j < P_order; j++)
          {
            for (unsigned k = 0; k < P_order; k++)
            {
              psi(P_order * P_order * i + P_order * j + k) =
                psi3[i] * psi2[j] * psi1[k];
            }
          }
        }
        break;
      }
      case 7:
      {
        // Call the OneDimensional Shape functions
        OneDimensionalLegendreShape<7>::calculate_nodal_positions();
        OneDimensionalLegendreShape<7> psi1(s[0]), psi2(s[1]), psi3(s[2]);

        // Now let's loop over the nodal points in the element
        // and copy the values back in
        for (unsigned i = 0; i < P_order; i++)
        {
          for (unsigned j = 0; j < P_order; j++)
          {
            for (unsigned k = 0; k < P_order; k++)
            {
              psi(P_order * P_order * i + P_order * j + k) =
                psi3[i] * psi2[j] * psi1[k];
            }
          }
        }
        break;
      }
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       polynomial order for shape functions.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
  }

  //=======================================================================
  /// Derivatives of shape functions for PRefineableQElement<DIM>
  //=======================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<3, INITIAL_NNODE_1D>::dshape_local(
    const Vector<double>& s, Shape& psi, DShape& dpsids) const
  {
    switch (P_order)
    {
      case 2:
      {
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<2>::calculate_nodal_positions();
        OneDimensionalLegendreShape<2> psi1(s[0]), psi2(s[1]), psi3(s[2]);
        OneDimensionalLegendreDShape<2> dpsi1ds(s[0]), dpsi2ds(s[1]),
          dpsi3ds(s[2]);

        // Index for the shape functions
        unsigned index = 0;
        // Loop over shape functions in element
        for (unsigned i = 0; i < P_order; i++)
        {
          for (unsigned j = 0; j < P_order; j++)
          {
            for (unsigned k = 0; k < P_order; k++)
            {
              // Assign the values
              dpsids(index, 0) = psi3[i] * psi2[j] * dpsi1ds[k];
              dpsids(index, 1) = psi3[i] * dpsi2ds[j] * psi1[k];
              dpsids(index, 2) = dpsi3ds[i] * psi2[j] * psi1[k];
              psi[index] = psi3[i] * psi2[j] * psi1[k];
              // Increment the index
              ++index;
            }
          }
        }
        break;
      }
      case 3:
      {
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<3>::calculate_nodal_positions();
        OneDimensionalLegendreShape<3> psi1(s[0]), psi2(s[1]), psi3(s[2]);
        OneDimensionalLegendreDShape<3> dpsi1ds(s[0]), dpsi2ds(s[1]),
          dpsi3ds(s[2]);

        // Index for the shape functions
        unsigned index = 0;
        // Loop over shape functions in element
        for (unsigned i = 0; i < P_order; i++)
        {
          for (unsigned j = 0; j < P_order; j++)
          {
            for (unsigned k = 0; k < P_order; k++)
            {
              // Assign the values
              dpsids(index, 0) = psi3[i] * psi2[j] * dpsi1ds[k];
              dpsids(index, 1) = psi3[i] * dpsi2ds[j] * psi1[k];
              dpsids(index, 2) = dpsi3ds[i] * psi2[j] * psi1[k];
              psi[index] = psi3[i] * psi2[j] * psi1[k];
              // Increment the index
              ++index;
            }
          }
        }
        break;
      }
      case 4:
      {
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<4>::calculate_nodal_positions();
        OneDimensionalLegendreShape<4> psi1(s[0]), psi2(s[1]), psi3(s[2]);
        OneDimensionalLegendreDShape<4> dpsi1ds(s[0]), dpsi2ds(s[1]),
          dpsi3ds(s[2]);

        // Index for the shape functions
        unsigned index = 0;
        // Loop over shape functions in element
        for (unsigned i = 0; i < P_order; i++)
        {
          for (unsigned j = 0; j < P_order; j++)
          {
            for (unsigned k = 0; k < P_order; k++)
            {
              // Assign the values
              dpsids(index, 0) = psi3[i] * psi2[j] * dpsi1ds[k];
              dpsids(index, 1) = psi3[i] * dpsi2ds[j] * psi1[k];
              dpsids(index, 2) = dpsi3ds[i] * psi2[j] * psi1[k];
              psi[index] = psi3[i] * psi2[j] * psi1[k];
              // Increment the index
              ++index;
            }
          }
        }
        break;
      }
      case 5:
      {
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<5>::calculate_nodal_positions();
        OneDimensionalLegendreShape<5> psi1(s[0]), psi2(s[1]), psi3(s[2]);
        OneDimensionalLegendreDShape<5> dpsi1ds(s[0]), dpsi2ds(s[1]),
          dpsi3ds(s[2]);

        // Index for the shape functions
        unsigned index = 0;
        // Loop over shape functions in element
        for (unsigned i = 0; i < P_order; i++)
        {
          for (unsigned j = 0; j < P_order; j++)
          {
            for (unsigned k = 0; k < P_order; k++)
            {
              // Assign the values
              dpsids(index, 0) = psi3[i] * psi2[j] * dpsi1ds[k];
              dpsids(index, 1) = psi3[i] * dpsi2ds[j] * psi1[k];
              dpsids(index, 2) = dpsi3ds[i] * psi2[j] * psi1[k];
              psi[index] = psi3[i] * psi2[j] * psi1[k];
              // Increment the index
              ++index;
            }
          }
        }
        break;
      }
      case 6:
      {
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<6>::calculate_nodal_positions();
        OneDimensionalLegendreShape<6> psi1(s[0]), psi2(s[1]), psi3(s[2]);
        OneDimensionalLegendreDShape<6> dpsi1ds(s[0]), dpsi2ds(s[1]),
          dpsi3ds(s[2]);

        // Index for the shape functions
        unsigned index = 0;
        // Loop over shape functions in element
        for (unsigned i = 0; i < P_order; i++)
        {
          for (unsigned j = 0; j < P_order; j++)
          {
            for (unsigned k = 0; k < P_order; k++)
            {
              // Assign the values
              dpsids(index, 0) = psi3[i] * psi2[j] * dpsi1ds[k];
              dpsids(index, 1) = psi3[i] * dpsi2ds[j] * psi1[k];
              dpsids(index, 2) = dpsi3ds[i] * psi2[j] * psi1[k];
              psi[index] = psi3[i] * psi2[j] * psi1[k];
              // Increment the index
              ++index;
            }
          }
        }
        break;
      }
      case 7:
      {
        // Call the shape functions and derivatives
        OneDimensionalLegendreShape<7>::calculate_nodal_positions();
        OneDimensionalLegendreShape<7> psi1(s[0]), psi2(s[1]), psi3(s[2]);
        OneDimensionalLegendreDShape<7> dpsi1ds(s[0]), dpsi2ds(s[1]),
          dpsi3ds(s[2]);

        // Index for the shape functions
        unsigned index = 0;
        // Loop over shape functions in element
        for (unsigned i = 0; i < P_order; i++)
        {
          for (unsigned j = 0; j < P_order; j++)
          {
            for (unsigned k = 0; k < P_order; k++)
            {
              // Assign the values
              dpsids(index, 0) = psi3[i] * psi2[j] * dpsi1ds[k];
              dpsids(index, 1) = psi3[i] * dpsi2ds[j] * psi1[k];
              dpsids(index, 2) = dpsi3ds[i] * psi2[j] * psi1[k];
              psi[index] = psi3[i] * psi2[j] * psi1[k];
              // Increment the index
              ++index;
            }
          }
        }
        break;
      }
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       polynomial order for shape functions.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
  }

  //=======================================================================
  /// Second derivatives of shape functions for PRefineableQElement<DIM>
  ///  d2psids(i,0) = \f$ d^2 \psi_j / d s^2 \f$
  //=======================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<3, INITIAL_NNODE_1D>::d2shape_local(
    const Vector<double>& s, Shape& psi, DShape& dpsids, DShape& d2psids) const
  {
    std::ostringstream error_message;
    error_message
      << "\nd2shape_local currently not implemented for this element\n";
    throw OomphLibError(
      error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  //=======================================================================
  /// Rebuild the element from nodes found in its sons
  /// Adjusts its p-order to be the maximum of its sons' p-orders
  //=======================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<3, INITIAL_NNODE_1D>::rebuild_from_sons(
    Mesh*& mesh_pt)
  {
    // Get p-orders of sons
    unsigned n_sons = this->tree_pt()->nsons();
    Vector<unsigned> son_p_order(n_sons);
    unsigned max_son_p_order = 0;
    for (unsigned ison = 0; ison < n_sons; ison++)
    {
      PRefineableElement* el_pt = dynamic_cast<PRefineableElement*>(
        this->tree_pt()->son_pt(ison)->object_pt());
      son_p_order[ison] = el_pt->p_order();
      if (son_p_order[ison] > max_son_p_order)
        max_son_p_order = son_p_order[ison];
    }

    unsigned old_Nnode = this->nnode();
    unsigned old_P_order = this->p_order();
    // Set p-order of the element
    this->p_order() = max_son_p_order;

    // Change integration scheme
    delete this->integral_pt();
    switch (this->p_order())
    {
      case 2:
        this->set_integration_scheme(new GaussLobattoLegendre<3, 2>);
        break;
      case 3:
        this->set_integration_scheme(new GaussLobattoLegendre<3, 3>);
        break;
      case 4:
        this->set_integration_scheme(new GaussLobattoLegendre<3, 4>);
        break;
      case 5:
        this->set_integration_scheme(new GaussLobattoLegendre<3, 5>);
        break;
      case 6:
        this->set_integration_scheme(new GaussLobattoLegendre<3, 6>);
        break;
      case 7:
        this->set_integration_scheme(new GaussLobattoLegendre<3, 7>);
        break;
      default:
        std::ostringstream error_message;
        error_message << "\nERROR: Exceeded maximum polynomial order for";
        error_message << "\n       integration scheme.\n";
        throw OomphLibError(
          error_message.str(),
          "PRefineableQElement<3,INITIAL_NNODE_1D>::rebuild_from_sons()",
          OOMPH_EXCEPTION_LOCATION);
    }

    // Back up pointers to old nodes before they are lost
    Vector<Node*> old_node_pt(old_Nnode);
    for (unsigned n = 0; n < old_Nnode; n++)
    {
      old_node_pt[n] = this->node_pt(n);
    }

    // Allocate new space for Nodes (at the element level)
    this->set_n_node(this->p_order() * this->p_order() * this->p_order());

    // Copy vertex nodes which were populated in the pre-build
    this->node_pt(0) = old_node_pt[0];
    this->node_pt(this->p_order() - 1) = old_node_pt[old_P_order - 1];
    this->node_pt(this->p_order() * (this->p_order() - 1)) =
      old_node_pt[(old_P_order) * (old_P_order - 1)];
    this->node_pt(this->p_order() * this->p_order() - 1) =
      old_node_pt[(old_P_order) * (old_P_order)-1];
    this->node_pt(this->p_order() * this->p_order() * (this->p_order() - 1)) =
      old_node_pt[old_P_order * old_P_order * (old_P_order - 1)];
    this->node_pt((this->p_order() * this->p_order() + 1) *
                  (this->p_order() - 1)) =
      old_node_pt[(old_P_order * old_P_order + 1) * (old_P_order - 1)];
    this->node_pt(this->p_order() * (this->p_order() + 1) *
                  (this->p_order() - 1)) =
      old_node_pt[old_P_order * (old_P_order + 1) * (old_P_order - 1)];
    this->node_pt(this->p_order() * this->p_order() * this->p_order() - 1) =
      old_node_pt[old_P_order * old_P_order * old_P_order - 1];

    // Copy midpoint nodes from sons if new p-order is odd
    if (this->p_order() % 2 == 1)
    {
      // Work out which is midpoint node
      unsigned n_p = this->p_order();
      unsigned m = (n_p - 1) / 2;

      unsigned s0space = m;
      unsigned s1space = m * n_p;
      unsigned s2space = m * n_p * n_p;

      // Back face
      // DB edge
      this->node_pt(1 * s0space + 0 * s1space + 0 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::RDB)->object_pt())
          ->vertex_node_pt(0);
      // LB edge
      this->node_pt(0 * s0space + 1 * s1space + 0 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::LUB)->object_pt())
          ->vertex_node_pt(0);
      // B centre
      this->node_pt(1 * s0space + 1 * s1space + 0 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::RUB)->object_pt())
          ->vertex_node_pt(0);
      // RB edge
      this->node_pt(2 * s0space + 1 * s1space + 0 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::RUB)->object_pt())
          ->vertex_node_pt(1);
      // UB edge
      this->node_pt(1 * s0space + 2 * s1space + 0 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::RUB)->object_pt())
          ->vertex_node_pt(2);

      // Mid-way between between back and front faces
      // LD edge
      this->node_pt(0 * s0space + 0 * s1space + 1 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::LDF)->object_pt())
          ->vertex_node_pt(0);
      // D centre
      this->node_pt(1 * s0space + 0 * s1space + 1 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::LDF)->object_pt())
          ->vertex_node_pt(1);
      // RD edge
      this->node_pt(2 * s0space + 0 * s1space + 1 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::RDF)->object_pt())
          ->vertex_node_pt(1);
      // L centre
      this->node_pt(0 * s0space + 1 * s1space + 1 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::LUF)->object_pt())
          ->vertex_node_pt(0);
      // Centre
      this->node_pt(1 * s0space + 1 * s1space + 1 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::RUF)->object_pt())
          ->vertex_node_pt(0);
      // R centre
      this->node_pt(2 * s0space + 1 * s1space + 1 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::RUF)->object_pt())
          ->vertex_node_pt(1);
      // LU edge
      this->node_pt(0 * s0space + 2 * s1space + 1 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::LUF)->object_pt())
          ->vertex_node_pt(2);
      // U center
      this->node_pt(1 * s0space + 2 * s1space + 1 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::RUF)->object_pt())
          ->vertex_node_pt(2);
      // RU edge
      this->node_pt(2 * s0space + 2 * s1space + 1 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::RUF)->object_pt())
          ->vertex_node_pt(3);

      // Front face
      // DF edge
      this->node_pt(1 * s0space + 0 * s1space + 2 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::LDF)->object_pt())
          ->vertex_node_pt(5);
      // LF edge
      this->node_pt(0 * s0space + 1 * s1space + 2 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::LUF)->object_pt())
          ->vertex_node_pt(4);
      // F centre
      this->node_pt(1 * s0space + 1 * s1space + 2 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::RUF)->object_pt())
          ->vertex_node_pt(4);
      // RF edge
      this->node_pt(2 * s0space + 1 * s1space + 2 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::RUF)->object_pt())
          ->vertex_node_pt(5);
      // UF edge
      this->node_pt(1 * s0space + 2 * s1space + 2 * s2space) =
        dynamic_cast<RefineableQElement<3>*>(
          octree_pt()->son_pt(OcTreeNames::RUF)->object_pt())
          ->vertex_node_pt(6);
    }

    // The timestepper should be the same for all nodes and node 0 should
    // never be deleted.
    if (this->node_pt(0) == 0)
    {
      throw OomphLibError(
        "The Corner node (0) does not exist",
        "PRefineableQElement<3,INITIAL_NNODE_1D>::rebuild_from_sons()",
        OOMPH_EXCEPTION_LOCATION);
    }

    TimeStepper* time_stepper_pt = this->node_pt(0)->time_stepper_pt();
    unsigned ntstorage = time_stepper_pt->ntstorage();

    unsigned jnod = 0;
    Vector<double> s_fraction(3), s(3);
    // Loop over the nodes in the element
    unsigned n_p = this->nnode_1d();
    for (unsigned i0 = 0; i0 < n_p; i0++)
    {
      // Get the fractional position of the node
      s_fraction[0] = this->local_one_d_fraction_of_node(i0, 0);
      // Local coordinate
      s[0] = -1.0 + 2.0 * s_fraction[0];

      for (unsigned i1 = 0; i1 < n_p; i1++)
      {
        // Get the fractional position of the node in the direction of s[1]
        s_fraction[1] = this->local_one_d_fraction_of_node(i1, 1);
        // Local coordinate in father element
        s[1] = -1.0 + 2.0 * s_fraction[1];

        for (unsigned i2 = 0; i2 < n_p; i2++)
        {
          // Get the fractional position of the node in the direction of s[2]
          s_fraction[2] = this->local_one_d_fraction_of_node(i2, 2);
          // Local coordinate in father element
          s[2] = -1.0 + 2.0 * s_fraction[2];

          // Set the local node number
          jnod = i0 + n_p * i1 + n_p * n_p * i2;

          // Initialise flag: So far, this node hasn't been built
          // or copied yet
          bool node_done = false;

          // Get the pointer to the node in this element, returns NULL
          // if there is not node
          Node* created_node_pt = this->get_node_at_local_coordinate(s);

          // Does this node already exist in this element?
          //----------------------------------------------
          if (created_node_pt != 0)
          {
            // Copy node across
            this->node_pt(jnod) = created_node_pt;

            // Node has been created by copy
            node_done = true;
          }
          // Node does not exist in this element but might already
          //------------------------------------------------------
          // have been created by neighbouring elements
          //-------------------------------------------
          else
          {
            // Was the node created by one of its neighbours
            // Whether or not the node lies on an edge can be calculated
            // by from the fractional position
            bool is_periodic = false;
            created_node_pt =
              node_created_by_neighbour(s_fraction, is_periodic);

            // If the node was so created, assign the pointers
            if (created_node_pt != 0)
            {
              // If the node is periodic
              if (is_periodic)
              {
                throw OomphLibError("Cannot handle periodic nodes yet",
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
              // Non-periodic case, just set the pointer
              else
              {
                this->node_pt(jnod) = created_node_pt;
              }
              // Node has been created
              node_done = true;
            }
          } // Node does not exist in this element

          // Node has not been built anywhere ---> build it here
          if (!node_done)
          {
            // First, find the son element in which the node should live

            // Find coordinates in the sons
            Vector<double> s_in_son(3);
            // using namespace OcTreeNames;
            int son = -10;
            // On the left
            if (s_fraction[0] < 0.5)
            {
              // On the bottom
              if (s_fraction[1] < 0.5)
              {
                // On the back
                if (s_fraction[2] < 0.5)
                {
                  // It's the LDB son
                  son = OcTreeNames::LDB;
                  s_in_son[0] = -1.0 + 4.0 * s_fraction[0];
                  s_in_son[1] = -1.0 + 4.0 * s_fraction[1];
                  s_in_son[2] = -1.0 + 4.0 * s_fraction[2];
                }
                // On the front
                else
                {
                  // It's the LDF son
                  son = OcTreeNames::LDF;
                  s_in_son[0] = -1.0 + 4.0 * s_fraction[0];
                  s_in_son[1] = -1.0 + 4.0 * s_fraction[1];
                  s_in_son[2] = -1.0 + 4.0 * (s_fraction[2] - 0.5);
                }
              }
              // On the top
              else
              {
                // On the back
                if (s[2] < 0.0)
                {
                  // It's the LUB son
                  son = OcTreeNames::LUB;
                  s_in_son[0] = -1.0 + 4.0 * s_fraction[0];
                  s_in_son[1] = -1.0 + 4.0 * (s_fraction[1] - 0.5);
                  s_in_son[2] = -1.0 + 4.0 * s_fraction[2];
                }
                // On the front
                else
                {
                  // It's the LUF son
                  son = OcTreeNames::LUF;
                  s_in_son[0] = -1.0 + 4.0 * s_fraction[0];
                  s_in_son[1] = -1.0 + 4.0 * (s_fraction[1] - 0.5);
                  s_in_son[2] = -1.0 + 4.0 * (s_fraction[2] - 0.5);
                }
              }
            }
            // On the right
            else
            {
              // On the bottom
              if (s[1] < 0.0)
              {
                // On the back
                if (s[2] < 0.0)
                {
                  // It's the RDB son
                  son = OcTreeNames::RDB;
                  s_in_son[0] = -1.0 + 4.0 * (s_fraction[0] - 0.5);
                  s_in_son[1] = -1.0 + 4.0 * s_fraction[1];
                  s_in_son[2] = -1.0 + 4.0 * s_fraction[2];
                }
                // On the front
                else
                {
                  // It's the RDF son
                  son = OcTreeNames::RDF;
                  s_in_son[0] = -1.0 + 4.0 * (s_fraction[0] - 0.5);
                  s_in_son[1] = -1.0 + 4.0 * s_fraction[1];
                  s_in_son[2] = -1.0 + 4.0 * (s_fraction[2] - 0.5);
                }
              }
              // On the top
              else
              {
                // On the back
                if (s[2] < 0.0)
                {
                  // It's the RUB son
                  son = OcTreeNames::RUB;
                  s_in_son[0] = -1.0 + 4.0 * (s_fraction[0] - 0.5);
                  s_in_son[1] = -1.0 + 4.0 * (s_fraction[1] - 0.5);
                  s_in_son[2] = -1.0 + 4.0 * s_fraction[2];
                }
                // On the front
                else
                {
                  // It's the RUF son
                  son = OcTreeNames::RUF;
                  s_in_son[0] = -1.0 + 4.0 * (s_fraction[0] - 0.5);
                  s_in_son[1] = -1.0 + 4.0 * (s_fraction[1] - 0.5);
                  s_in_son[2] = -1.0 + 4.0 * (s_fraction[2] - 0.5);
                }
              }
            }

            // Get the pointer to the son element in which the new node
            // would live
            PRefineableQElement<3, INITIAL_NNODE_1D>* son_el_pt =
              dynamic_cast<PRefineableQElement<3, INITIAL_NNODE_1D>*>(
                this->tree_pt()->son_pt(son)->object_pt());

            // If we are rebuilding, then worry about the boundary conditions
            // Find the boundary of the node
            // Initially none
            int boundary = Tree::OMEGA;
            // If we are on the left face
            if (i0 == 0)
            {
              boundary = OcTreeNames::L;
            }
            // If we are on the right face
            else if (i0 == n_p - 1)
            {
              boundary = OcTreeNames::R;
            }

            // If we are on the bottom face
            if (i1 == 0)
            {
              // If we already have already set the boundary, we're on an edge
              switch (boundary)
              {
                case OcTreeNames::L:
                  boundary = OcTreeNames::LD;
                  break;
                case OcTreeNames::R:
                  boundary = OcTreeNames::RD;
                  break;
                // Boundary not set
                default:
                  boundary = OcTreeNames::D;
                  break;
              }
            }
            // If we are the top face
            else if (i1 == n_p - 1)
            {
              // If we already have a boundary
              switch (boundary)
              {
                case OcTreeNames::L:
                  boundary = OcTreeNames::LU;
                  break;
                case OcTreeNames::R:
                  boundary = OcTreeNames::RU;
                  break;
                default:
                  boundary = OcTreeNames::U;
                  break;
              }
            }

            // If we are on the back face
            if (i2 == 0)
            {
              // If we already have already set the boundary, we're on an edge
              switch (boundary)
              {
                case OcTreeNames::L:
                  boundary = OcTreeNames::LB;
                  break;
                case OcTreeNames::LD:
                  boundary = OcTreeNames::LDB;
                  break;
                case OcTreeNames::LU:
                  boundary = OcTreeNames::LUB;
                  break;
                case OcTreeNames::R:
                  boundary = OcTreeNames::RB;
                  break;
                case OcTreeNames::RD:
                  boundary = OcTreeNames::RDB;
                  break;
                case OcTreeNames::RU:
                  boundary = OcTreeNames::RUB;
                  break;
                case OcTreeNames::D:
                  boundary = OcTreeNames::DB;
                  break;
                case OcTreeNames::U:
                  boundary = OcTreeNames::UB;
                  break;
                // Boundary not set
                default:
                  boundary = OcTreeNames::B;
                  break;
              }
            }
            // If we are the front face
            else if (i2 == n_p - 1)
            {
              // If we already have a boundary
              switch (boundary)
              {
                case OcTreeNames::L:
                  boundary = OcTreeNames::LF;
                  break;
                case OcTreeNames::LD:
                  boundary = OcTreeNames::LDF;
                  break;
                case OcTreeNames::LU:
                  boundary = OcTreeNames::LUF;
                  break;
                case OcTreeNames::R:
                  boundary = OcTreeNames::RF;
                  break;
                case OcTreeNames::RD:
                  boundary = OcTreeNames::RDF;
                  break;
                case OcTreeNames::RU:
                  boundary = OcTreeNames::RUF;
                  break;
                case OcTreeNames::D:
                  boundary = OcTreeNames::DF;
                  break;
                case OcTreeNames::U:
                  boundary = OcTreeNames::UF;
                  break;
                default:
                  boundary = OcTreeNames::F;
                  break;
              }
            }

            // set of boundaries that this edge in the son lives on
            std::set<unsigned> boundaries;

            // Now get the boundary conditions from the son
            // The boundaries will be common to the son because there can be
            // no rotations here
            if (boundary != Tree::OMEGA)
            {
              son_el_pt->get_boundaries(boundary, boundaries);
            }

            // If the node lives on a boundary:
            // Construct a boundary node,
            // Get boundary conditions and
            // update all lookup schemes
            if (boundaries.size() > 0)
            {
              // Construct the new node
              created_node_pt =
                this->construct_boundary_node(jnod, time_stepper_pt);

              // Get the boundary conditions from the son
              Vector<int> bound_cons(this->ncont_interpolated_values());
              son_el_pt->get_bcs(boundary, bound_cons);

              // Loop over the values and pin if necessary
              unsigned nval = created_node_pt->nvalue();
              for (unsigned k = 0; k < nval; k++)
              {
                if (bound_cons[k])
                {
                  created_node_pt->pin(k);
                }
              }

              // Solid node? If so, deal with the positional boundary
              // conditions:
              SolidNode* solid_node_pt =
                dynamic_cast<SolidNode*>(created_node_pt);
              if (solid_node_pt != 0)
              {
                // Get the positional boundary conditions from the father:
                unsigned n_dim = created_node_pt->ndim();
                Vector<int> solid_bound_cons(n_dim);
                RefineableSolidQElement<3>* son_solid_el_pt =
                  dynamic_cast<RefineableSolidQElement<3>*>(son_el_pt);
#ifdef PARANOID
                if (son_solid_el_pt == 0)
                {
                  std::string error_message =
                    "We have a SolidNode outside a refineable SolidElement\n";
                  error_message +=
                    "during mesh refinement -- this doesn't make sense\n";

                  throw OomphLibError(error_message,
                                      "PRefineableQElement<3,INITIAL_NNODE_1D>:"
                                      ":rebuild_from_sons()",
                                      OOMPH_EXCEPTION_LOCATION);
                }
#endif
                son_solid_el_pt->get_solid_bcs(boundary, solid_bound_cons);

                // Loop over the positions and pin, if necessary
                for (unsigned k = 0; k < n_dim; k++)
                {
                  if (solid_bound_cons[k])
                  {
                    solid_node_pt->pin_position(k);
                  }
                }
              } // End of if solid_node_pt


              // Next we update the boundary look-up schemes
              // Loop over the boundaries stored in the set
              for (std::set<unsigned>::iterator it = boundaries.begin();
                   it != boundaries.end();
                   ++it)
              {
                // Add the node to the boundary
                mesh_pt->add_boundary_node(*it, created_node_pt);

                // If we have set an intrinsic coordinate on this
                // mesh boundary then it must also be interpolated on
                // the new node
                // Now interpolate the intrinsic boundary coordinate
                if (mesh_pt->boundary_coordinate_exists(*it) == true)
                {
                  Vector<double> zeta(2);
                  son_el_pt->interpolated_zeta_on_face(
                    *it, boundary, s_in_son, zeta);

                  created_node_pt->set_coordinates_on_boundary(*it, zeta);
                }
              }
            }
            // Otherwise the node is not on a Mesh boundary
            // and we create a normal "bulk" node
            else
            {
              // Construct the new node
              created_node_pt = this->construct_node(jnod, time_stepper_pt);
            }

            // Now we set the position and values at the newly created node

            // In the first instance use macro element or FE representation
            // to create past and present nodal positions.
            // (THIS STEP SHOULD NOT BE SKIPPED FOR ALGEBRAIC
            // ELEMENTS AS NOT ALL OF THEM NECESSARILY IMPLEMENT
            // NONTRIVIAL NODE UPDATE FUNCTIONS. CALLING
            // THE NODE UPDATE FOR SUCH ELEMENTS/NODES WILL LEAVE
            // THEIR NODAL POSITIONS WHERE THEY WERE (THIS IS APPROPRIATE
            // ONCE THEY HAVE BEEN GIVEN POSITIONS) BUT WILL
            // NOT ASSIGN SENSIBLE INITIAL POSITONS!

            // Loop over # of history values
            // Loop over # of history values
            for (unsigned t = 0; t < ntstorage; t++)
            {
              using namespace QuadTreeNames;
              // Get the position from the son
              Vector<double> x_prev(3);

              // Now let's fill in the value
              son_el_pt->get_x(t, s_in_son, x_prev);
              for (unsigned i = 0; i < 3; i++)
              {
                created_node_pt->x(t, i) = x_prev[i];
              }
            }

            // Now set up the values
            // Loop over all history values
            for (unsigned t = 0; t < ntstorage; t++)
            {
              // Get values from father element
              // Note: get_interpolated_values() sets Vector size itself.
              Vector<double> prev_values;
              son_el_pt->get_interpolated_values(t, s_in_son, prev_values);

              // Initialise the values at the new node
              for (unsigned k = 0; k < created_node_pt->nvalue(); k++)
              {
                created_node_pt->set_value(t, k, prev_values[k]);
              }
            }

            // Add the node to the mesh
            mesh_pt->add_node_pt(created_node_pt);

            // Check if the element is an algebraic element
            AlgebraicElementBase* alg_el_pt =
              dynamic_cast<AlgebraicElementBase*>(this);

            // If we do have an algebraic element
            if (alg_el_pt != 0)
            {
              std::string error_message =
                "Have not implemented rebuilding from sons for";
              error_message += "Algebraic p-refineable elements yet\n";

              throw OomphLibError(
                error_message,
                "PRefineableQElement<3,INITIAL_NNODE_1D>::rebuild_from_sons()",
                OOMPH_EXCEPTION_LOCATION);
            }

          } // End of the case when we build the node ourselves
        }
      }
    }
  }

  //=================================================================
  /// Check inter-element continuity of
  /// - nodal positions
  /// - (nodally) interpolated function values
  /// Overloaded to not check differences in the value. Mortaring
  /// doesn't enforce strong continuity between elements.
  //====================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<3, INITIAL_NNODE_1D>::check_integrity(
    double& max_error)
  {
    // Not yet implemented -- doing nothing for now.

    // Dummy set max_error to 0
    max_error = 0.0;
    return;
  }

  //=================================================================
  /// Internal function to set up the hanging nodes on a particular
  /// edge of the element. Implements the mortar method.
  //=================================================================
  template<unsigned INITIAL_NNODE_1D>
  void PRefineableQElement<3, INITIAL_NNODE_1D>::oc_hang_helper(
    const int& value_id, const int& my_face, std::ofstream& output_hangfile)
  {
    using namespace OcTreeNames;

    Vector<unsigned> translate_s(3);
    Vector<double> s_lo_neigh(3);
    Vector<double> s_hi_neigh(3);
    int neigh_face, diff_level;
    bool in_neighbouring_tree;

    // Find pointer to neighbour in this direction
    OcTree* neigh_pt;
    neigh_pt = this->octree_pt()->gteq_face_neighbour(my_face,
                                                      translate_s,
                                                      s_lo_neigh,
                                                      s_hi_neigh,
                                                      neigh_face,
                                                      diff_level,
                                                      in_neighbouring_tree);

    // Work out master/dependent faces
    //----------------------------

    // Set up booleans
    // bool h_type_master = false;
    bool h_type_dependent = false;
    // bool p_type_master = false;
    bool p_type_dependent = false;

    // Neighbour exists
    if (neigh_pt != 0)
    {
      // Check if neighbour is bigger than me
      if (diff_level != 0)
      {
        // Dependent at h-type non-conformity
        h_type_dependent = true;
      }
      // Check if neighbour is the same size as me
      else if (neigh_pt->nsons() == 0)
      {
        // Neighbour is same size as me
        // Find p-orders of each element
        unsigned my_p_order =
          dynamic_cast<PRefineableQElement<3, INITIAL_NNODE_1D>*>(this)
            ->p_order();
        unsigned neigh_p_order =
          dynamic_cast<PRefineableQElement<3, INITIAL_NNODE_1D>*>(
            neigh_pt->object_pt())
            ->p_order();

        // Check for p-type non-conformity
        if (neigh_p_order == my_p_order)
        {
          // At a conforming interface
        }
        else if (neigh_p_order < my_p_order)
        {
          // Dependent at p-type non-conformity
          p_type_dependent = true;
        }
        else
        {
          // Master at p-type non-conformity
          // p_type_master = true;
        }
      }
      // Neighbour must be smaller than me
      else
      {
        // Master at h-type non-conformity
        // h_type_master = true;
      }
    }
    else
    {
      // Face is on a boundary
    }

    // Work out master/dependent edges
    //----------------------------
    if (h_type_dependent || p_type_dependent || true)
    {
      // Get edges of face and the face that shares that edge
      Vector<unsigned> face_edge(4), face_edge_other_face(4);
      switch (my_face)
      {
        case L:
          face_edge[0] = LB;
          face_edge_other_face[0] = B;
          face_edge[1] = LF;
          face_edge_other_face[1] = F;
          face_edge[2] = LD;
          face_edge_other_face[2] = D;
          face_edge[3] = LU;
          face_edge_other_face[3] = U;
          break;
        case R:
          face_edge[0] = RB;
          face_edge_other_face[0] = B;
          face_edge[1] = RF;
          face_edge_other_face[1] = F;
          face_edge[2] = RD;
          face_edge_other_face[2] = D;
          face_edge[3] = RU;
          face_edge_other_face[3] = U;
          break;
        case U:
          face_edge[0] = UB;
          face_edge_other_face[0] = B;
          face_edge[1] = UF;
          face_edge_other_face[1] = F;
          face_edge[2] = LU;
          face_edge_other_face[2] = L;
          face_edge[3] = RU;
          face_edge_other_face[3] = R;
          break;
        case D:
          face_edge[0] = DB;
          face_edge_other_face[0] = B;
          face_edge[1] = DF;
          face_edge_other_face[1] = F;
          face_edge[2] = LD;
          face_edge_other_face[2] = L;
          face_edge[3] = RD;
          face_edge_other_face[3] = R;
          break;
        case B:
          face_edge[0] = DB;
          face_edge_other_face[0] = D;
          face_edge[1] = UB;
          face_edge_other_face[1] = U;
          face_edge[2] = LB;
          face_edge_other_face[2] = L;
          face_edge[3] = RB;
          face_edge_other_face[3] = R;
          break;
        case F:
          face_edge[0] = DF;
          face_edge_other_face[0] = D;
          face_edge[1] = UF;
          face_edge_other_face[1] = U;
          face_edge[2] = LF;
          face_edge_other_face[2] = L;
          face_edge[3] = RF;
          face_edge_other_face[3] = R;
          break;
        default:
          throw OomphLibError(
            "my_face not L, R, D, U, B, F\n",
            "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
            OOMPH_EXCEPTION_LOCATION);
      }

      // Loop over edges of face
      for (unsigned i = 0; i < 4; i++)
      {
        // Get edge
        unsigned my_edge = face_edge[i];

        // Separate storage for edge mortaring
        OcTree* edge_neigh_pt = 0;
        Vector<unsigned> edge_translate_s(3);
        Vector<double> edge_s_lo_neigh(3);
        Vector<double> edge_s_hi_neigh(3);
        int neigh_edge = 0, edge_diff_level = 0;
        unsigned edge_p_order =
          dynamic_cast<PRefineableQElement<3, INITIAL_NNODE_1D>*>(this)
            ->p_order();

        // Temporary storage to keep track of master edge
        Vector<unsigned> tmp_edge_translate_s(3);
        Vector<double> tmp_edge_s_lo_neigh(3);
        Vector<double> tmp_edge_s_hi_neigh(3);
        int tmp_neigh_edge, tmp_edge_diff_level;

        // Initially return the zero-th true edge neighbour
        unsigned i_root_edge_neighbour = 0;

        // Initialise the total number of true edge neighbours
        unsigned nroot_edge_neighbour = 0;

        // Flag to keep track of master status of my_edge
        bool my_edge_is_master = true;
        // unsigned master_edge_index=0;

        // Keep searching until we've found the node or until we've checked
        // all available edge neighbours
        bool keep_searching = true;
        bool search_faces = false;
        bool first_face_searched = false;
        // unsigned index=0;
        while (keep_searching)
        {
          // Pointer to edge neighbouring OcTree
          OcTree* tmp_edge_neigh_pt;

          // Looking in edge neighbours that are also face neighbours
          // if(index==0 || index==1)
          if (search_faces)
          {
            if (!first_face_searched)
            {
              // Find pointer to neighbour in this direction
              tmp_edge_neigh_pt =
                this->octree_pt()->gteq_face_neighbour(my_face,
                                                       tmp_edge_translate_s,
                                                       tmp_edge_s_lo_neigh,
                                                       tmp_edge_s_hi_neigh,
                                                       tmp_neigh_edge,
                                                       tmp_edge_diff_level,
                                                       in_neighbouring_tree);

              // Mark first face as searched
              first_face_searched = true;
            }
            else
            {
              // Find pointer to neighbour in this direction
              tmp_edge_neigh_pt =
                this->octree_pt()->gteq_face_neighbour(face_edge_other_face[i],
                                                       tmp_edge_translate_s,
                                                       tmp_edge_s_lo_neigh,
                                                       tmp_edge_s_hi_neigh,
                                                       tmp_neigh_edge,
                                                       tmp_edge_diff_level,
                                                       in_neighbouring_tree);

              // Search is finally exhausted
              keep_searching = false;
            }
          }
          // Looking in a true edge neighbour
          else
          {
            // Find pointer to neighbour in this direction
            tmp_edge_neigh_pt =
              this->octree_pt()->gteq_true_edge_neighbour(my_edge,
                                                          i_root_edge_neighbour,
                                                          nroot_edge_neighbour,
                                                          tmp_edge_translate_s,
                                                          tmp_edge_s_lo_neigh,
                                                          tmp_edge_s_hi_neigh,
                                                          tmp_neigh_edge,
                                                          tmp_edge_diff_level);
          }

          // Set up booleans
          // bool h_type_edge_master = false;
          // bool h_type_edge_dependent  = false;
          // bool p_type_edge_master = false;
          // bool p_type_edge_dependent  = false;

          // Flag to check if we have a new edge master
          bool new_edge_master = false;

          // Edge neighbour exists
          if (tmp_edge_neigh_pt != 0)
          {
            // Check if neighbour is bigger than my biggest neighbour
            if (tmp_edge_diff_level < edge_diff_level)
            {
              // Dependent at h-type non-conformity
              // h_type_edge_dependent = true;
              new_edge_master = true;
              // Update edge_diff_level and p-order
              edge_diff_level = tmp_edge_diff_level;
              edge_p_order =
                dynamic_cast<PRefineableQElement<3, INITIAL_NNODE_1D>*>(
                  tmp_edge_neigh_pt->object_pt())
                  ->p_order();
            }
            // Check if neighbour is the same size as my biggest neighbour
            else if (tmp_edge_diff_level == edge_diff_level &&
                     tmp_edge_neigh_pt->nsons() == 0)
            {
              // Neighbour is same size as me
              // Find p-orders of each element
              // unsigned my_p_order =
              // dynamic_cast<PRefineableQElement<3,INITIAL_NNODE_1D>*>
              //  (this)->p_order();
              unsigned tmp_edge_neigh_p_order =
                dynamic_cast<PRefineableQElement<3, INITIAL_NNODE_1D>*>(
                  tmp_edge_neigh_pt->object_pt())
                  ->p_order();

              // Check for p-type non-conformity
              if (tmp_edge_neigh_p_order == edge_p_order)
              {
                // At a conforming interface
              }
              else if (tmp_edge_neigh_p_order < edge_p_order)
              {
                // Dependent at p-type non-conformity
                // p_type_edge_dependent = true;
                new_edge_master = true;
                // Update edge_diff_level and p-order
                edge_diff_level = tmp_edge_diff_level;
                edge_p_order =
                  dynamic_cast<PRefineableQElement<3, INITIAL_NNODE_1D>*>(
                    tmp_edge_neigh_pt->object_pt())
                    ->p_order();
              }
              else
              {
                // Master at p-type non-conformity
                // p_type_edge_master = true;
              }
            }
            // Neighbour must be smaller than me
            else
            {
              // Master at h-type non-conformity
              // h_type_edge_master = true;
            }
          }
          else
          {
            // Edge is on a boundary
          }

          // Update master neighbour information
          if (new_edge_master)
          {
            // Store new master edge information
            edge_neigh_pt = tmp_edge_neigh_pt;
            edge_translate_s = tmp_edge_translate_s;
            edge_s_lo_neigh = tmp_edge_s_lo_neigh;
            edge_s_hi_neigh = tmp_edge_s_hi_neigh;
            neigh_edge = tmp_neigh_edge;
            edge_diff_level = tmp_edge_diff_level;

            // Set this edge as dependent
            my_edge_is_master = false;
          }

          // Keep searching, but only if there are further edge neighbours
          // Try next root edge neighbour
          i_root_edge_neighbour++;

          // Increment counter
          // index++;

          // Have we exhausted the search over true edge neighbours
          // if (index>2 && i_root_edge_neighbour>=nroot_edge_neighbour)
          if (i_root_edge_neighbour >= nroot_edge_neighbour)
          {
            // keep_searching = false;
            // Extend search to face neighours (these are edge neighbours too)
            search_faces = true;
          }

        } // End of while keep searching over all face and true edge neighbours

        // Now do edge mortaring to enforce the mortar vertex matching condition
        if (!my_edge_is_master)
        {
          // Compute the active coordinate index along the this side of mortar
          unsigned active_coord_index;
          switch (my_edge)
          {
            case DB:
            case DF:
            case UB:
            case UF:
              active_coord_index = 0;
              break;
            case LB:
            case RB:
            case LF:
            case RF:
              active_coord_index = 1;
              break;
            case LD:
            case RD:
            case LU:
            case RU:
              active_coord_index = 2;
              break;
            default:
              throw OomphLibError(
                "Cannot transform coordinates",
                "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
                OOMPH_EXCEPTION_LOCATION);
          }

          // Get pointer to neighbouring master element (in p-refineable form)
          PRefineableQElement<3, INITIAL_NNODE_1D>* edge_neigh_obj_pt;
          edge_neigh_obj_pt =
            dynamic_cast<PRefineableQElement<3, INITIAL_NNODE_1D>*>(
              edge_neigh_pt->object_pt());

          // Create vector of master and dependent nodes
          //----------------------------------------
          Vector<Node*> master_node_pt, dependent_node_pt;
          Vector<unsigned> master_node_number, dependent_node_number;
          Vector<Vector<double>> dependent_node_s_fraction;

          // Number of nodes in one dimension
          const unsigned my_n_p = this->ninterpolating_node_1d(value_id);
          const unsigned neigh_n_p =
            edge_neigh_obj_pt->ninterpolating_node_1d(value_id);

          // Storage for pointers to the nodes and their numbers along the
          // master edge
          unsigned neighbour_node_number = 0;
          Node* neighbour_node_pt = 0;

          // Loop over nodes along the edge
          bool master_is_not_edge = false;
          for (unsigned i0 = 0; i0 < neigh_n_p; i0++)
          {
            const unsigned s0space = 1;
            const unsigned s1space = neigh_n_p;
            const unsigned s2space = neigh_n_p * neigh_n_p;

            // Find the neighbour's node
            switch (neigh_edge)
            {
              case DB:
                neighbour_node_number = i0 * s0space;
                neighbour_node_pt = edge_neigh_obj_pt->interpolating_node_pt(
                  neighbour_node_number, value_id);
                break;
              case UB:
                neighbour_node_number =
                  (neigh_n_p - 1) * s1space + i0 * s0space;
                neighbour_node_pt = edge_neigh_obj_pt->interpolating_node_pt(
                  neighbour_node_number, value_id);
                break;
              case DF:
                neighbour_node_number =
                  (neigh_n_p - 1) * s2space + i0 * s0space;
                neighbour_node_pt = edge_neigh_obj_pt->interpolating_node_pt(
                  neighbour_node_number, value_id);
                break;
              case UF:
                neighbour_node_number = (neigh_n_p - 1) * s1space +
                                        (neigh_n_p - 1) * s2space +
                                        i0 * s0space;
                neighbour_node_pt = edge_neigh_obj_pt->interpolating_node_pt(
                  neighbour_node_number, value_id);
                break;

              case LB:
                neighbour_node_number = i0 * s1space;
                neighbour_node_pt = edge_neigh_obj_pt->interpolating_node_pt(
                  neighbour_node_number, value_id);
                break;
              case RB:
                neighbour_node_number =
                  (neigh_n_p - 1) * s0space + i0 * s1space;
                neighbour_node_pt = edge_neigh_obj_pt->interpolating_node_pt(
                  neighbour_node_number, value_id);
                break;
              case LF:
                neighbour_node_number =
                  (neigh_n_p - 1) * s2space + i0 * s1space;
                neighbour_node_pt = edge_neigh_obj_pt->interpolating_node_pt(
                  neighbour_node_number, value_id);
                break;
              case RF:
                neighbour_node_number = (neigh_n_p - 1) * s0space +
                                        (neigh_n_p - 1) * s2space +
                                        i0 * s1space;
                neighbour_node_pt = edge_neigh_obj_pt->interpolating_node_pt(
                  neighbour_node_number, value_id);
                break;

              case LD:
                neighbour_node_number = i0 * s2space;
                neighbour_node_pt = edge_neigh_obj_pt->interpolating_node_pt(
                  neighbour_node_number, value_id);
                break;
              case RD:
                neighbour_node_number =
                  (neigh_n_p - 1) * s0space + i0 * s2space;
                neighbour_node_pt = edge_neigh_obj_pt->interpolating_node_pt(
                  neighbour_node_number, value_id);
                break;
              case LU:
                neighbour_node_number =
                  (neigh_n_p - 1) * s1space + i0 * s2space;
                neighbour_node_pt = edge_neigh_obj_pt->interpolating_node_pt(
                  neighbour_node_number, value_id);
                break;
              case RU:
                neighbour_node_number = (neigh_n_p - 1) * s0space +
                                        (neigh_n_p - 1) * s1space +
                                        i0 * s2space;
                neighbour_node_pt = edge_neigh_obj_pt->interpolating_node_pt(
                  neighbour_node_number, value_id);
                break;

              default:
                // Master 'edge' may be a face instead, so no need to throw an
                // error
                master_is_not_edge = true;
            }

            if (master_is_not_edge) break;

            // Set node as master node
            master_node_number.push_back(neighbour_node_number);
            master_node_pt.push_back(neighbour_node_pt);
          }

          // Now if edge is really a face (from an edge neighbour that isn't
          // a true edge neighbour) each node on the face is (potentially) a
          // master
          if (master_is_not_edge)
          {
            // Loop over nodes along the face
            for (unsigned i0 = 0; i0 < neigh_n_p; i0++)
            {
              // Loop over nodes along the face
              for (unsigned i1 = 0; i1 < neigh_n_p; i1++)
              {
                const unsigned s0space = 1;
                const unsigned s1space = neigh_n_p;
                const unsigned s2space = neigh_n_p * neigh_n_p;

                // Find the neighbour's node
                switch (neigh_edge)
                {
                  case B:
                    neighbour_node_number = i0 * s0space + i1 * s1space;
                    neighbour_node_pt =
                      edge_neigh_obj_pt->interpolating_node_pt(
                        neighbour_node_number, value_id);
                    break;
                  case F:
                    neighbour_node_number =
                      (neigh_n_p - 1) * s2space + i0 * s0space + i1 * s1space;
                    neighbour_node_pt =
                      edge_neigh_obj_pt->interpolating_node_pt(
                        neighbour_node_number, value_id);
                    break;
                  case D:
                    neighbour_node_number = i0 * s0space + i1 * s2space;
                    neighbour_node_pt =
                      edge_neigh_obj_pt->interpolating_node_pt(
                        neighbour_node_number, value_id);
                    break;
                  case U:
                    neighbour_node_number =
                      (neigh_n_p - 1) * s1space + i0 * s0space + i1 * s2space;
                    neighbour_node_pt =
                      edge_neigh_obj_pt->interpolating_node_pt(
                        neighbour_node_number, value_id);
                    break;
                  case L:
                    neighbour_node_number = i0 * s1space + i1 * s2space;
                    neighbour_node_pt =
                      edge_neigh_obj_pt->interpolating_node_pt(
                        neighbour_node_number, value_id);
                    break;
                  case R:
                    neighbour_node_number =
                      (neigh_n_p - 1) * s0space + i0 * s1space + i1 * s2space;
                    neighbour_node_pt =
                      edge_neigh_obj_pt->interpolating_node_pt(
                        neighbour_node_number, value_id);
                    break;

                  default:
                    throw OomphLibError("neigh_edge not recognised\n",
                                        "PRefineableQElement<3,INITIAL_NNODE_"
                                        "1D>::oc_hang_helper()",
                                        OOMPH_EXCEPTION_LOCATION);
                }

                // Set node as master node
                master_node_number.push_back(neighbour_node_number);
                master_node_pt.push_back(neighbour_node_pt);
              }
            }
          }

          // Storage for pointers to the local nodes and their numbers along my
          // edge
          unsigned local_node_number = 0;
          Node* local_node_pt = 0;

          // Loop over the nodes along my edge
          for (unsigned i0 = 0; i0 < my_n_p; i0++)
          {
            // Storage for the fractional position of the node in the element
            Vector<double> s_fraction(3);

            const unsigned s0space = 1;
            const unsigned s1space = my_n_p;
            const unsigned s2space = my_n_p * my_n_p;

            // Find the local node and the fractional position of the node
            // which depends on the edge, of course
            switch (my_edge)
            {
              case DB:
                s_fraction[0] =
                  local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
                s_fraction[1] = 0.0;
                s_fraction[2] = 0.0;
                local_node_number = i0 * s0space;
                local_node_pt =
                  this->interpolating_node_pt(local_node_number, value_id);
                break;
              case UB:
                s_fraction[0] =
                  local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
                s_fraction[1] = 1.0;
                s_fraction[2] = 0.0;
                local_node_number = (my_n_p - 1) * s1space + i0 * s0space;
                local_node_pt =
                  this->interpolating_node_pt(local_node_number, value_id);
                break;
              case DF:
                s_fraction[0] =
                  local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
                s_fraction[1] = 0.0;
                s_fraction[2] = 1.0;
                local_node_number = (my_n_p - 1) * s2space + i0 * s0space;
                local_node_pt =
                  this->interpolating_node_pt(local_node_number, value_id);
                break;
              case UF:
                s_fraction[0] =
                  local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
                s_fraction[1] = 1.0;
                s_fraction[2] = 1.0;
                local_node_number = (my_n_p - 1) * s1space +
                                    (my_n_p - 1) * s2space + i0 * s0space;
                local_node_pt =
                  this->interpolating_node_pt(local_node_number, value_id);
                break;

              case LB:
                s_fraction[0] = 0.0;
                s_fraction[1] =
                  local_one_d_fraction_of_interpolating_node(i0, 1, value_id);
                s_fraction[2] = 0.0;
                local_node_number = i0 * s1space;
                local_node_pt =
                  this->interpolating_node_pt(local_node_number, value_id);
                break;
              case RB:
                s_fraction[0] = 1.0;
                s_fraction[1] =
                  local_one_d_fraction_of_interpolating_node(i0, 1, value_id);
                s_fraction[2] = 0.0;
                local_node_number = (my_n_p - 1) * s0space + i0 * s1space;
                local_node_pt =
                  this->interpolating_node_pt(local_node_number, value_id);
                break;
              case LF:
                s_fraction[0] = 0.0;
                s_fraction[1] =
                  local_one_d_fraction_of_interpolating_node(i0, 1, value_id);
                s_fraction[2] = 1.0;
                local_node_number = (my_n_p - 1) * s2space + i0 * s1space;
                local_node_pt =
                  this->interpolating_node_pt(local_node_number, value_id);
                break;
              case RF:
                s_fraction[0] = 1.0;
                s_fraction[1] =
                  local_one_d_fraction_of_interpolating_node(i0, 1, value_id);
                s_fraction[2] = 1.0;
                local_node_number = (my_n_p - 1) * s0space +
                                    (my_n_p - 1) * s2space + i0 * s1space;
                local_node_pt =
                  this->interpolating_node_pt(local_node_number, value_id);
                break;

              case LD:
                s_fraction[0] = 0.0;
                s_fraction[1] = 0.0;
                s_fraction[2] =
                  local_one_d_fraction_of_interpolating_node(i0, 2, value_id);
                local_node_number = i0 * s2space;
                local_node_pt =
                  this->interpolating_node_pt(local_node_number, value_id);
                break;
              case RD:
                s_fraction[0] = 1.0;
                s_fraction[1] = 0.0;
                s_fraction[2] =
                  local_one_d_fraction_of_interpolating_node(i0, 2, value_id);
                local_node_number = (my_n_p - 1) * s0space + i0 * s2space;
                local_node_pt =
                  this->interpolating_node_pt(local_node_number, value_id);
                break;
              case LU:
                s_fraction[0] = 0.0;
                s_fraction[1] = 1.0;
                s_fraction[2] =
                  local_one_d_fraction_of_interpolating_node(i0, 2, value_id);
                local_node_number = (my_n_p - 1) * s1space + i0 * s2space;
                local_node_pt =
                  this->interpolating_node_pt(local_node_number, value_id);
                break;
              case RU:
                s_fraction[0] = 1.0;
                s_fraction[1] = 1.0;
                s_fraction[2] =
                  local_one_d_fraction_of_interpolating_node(i0, 2, value_id);
                local_node_number = (my_n_p - 1) * s0space +
                                    (my_n_p - 1) * s1space + i0 * s2space;
                local_node_pt =
                  this->interpolating_node_pt(local_node_number, value_id);
                break;

              default:
                throw OomphLibError(
                  "my_edge not recognised\n",
                  "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
                  OOMPH_EXCEPTION_LOCATION);
            }

            // Add node to vector of dependent element nodes
            dependent_node_number.push_back(local_node_number);
            dependent_node_pt.push_back(local_node_pt);

            // Store node's local fraction
            dependent_node_s_fraction.push_back(s_fraction);
          }

          // Store the number of dependent and master nodes for use later
          const unsigned n_dependent_nodes = dependent_node_pt.size();
          const unsigned n_master_nodes = master_node_pt.size();
          const unsigned dependent_element_nnode_1d = my_n_p;
          const unsigned master_element_nnode_1d = neigh_n_p;

          // Storage for master shapes
          Shape master_shapes(edge_neigh_obj_pt->ninterpolating_node(value_id));

          // Get master and dependent nodal positions
          //-------------------------------------
          Vector<double> dependent_nodal_position;
          Vector<double> dependent_weight(dependent_element_nnode_1d);
          Orthpoly::gll_nodes(dependent_element_nnode_1d,
                              dependent_nodal_position,
                              dependent_weight);
          Vector<double> master_nodal_position;
          Vector<double> master_weight(master_element_nnode_1d);
          Orthpoly::gll_nodes(
            master_element_nnode_1d, master_nodal_position, master_weight);

          // Apply the (1D) vertex matching condition
          //-----------------------------------------
          // Vertiex matching is ensured automatically in cases where there is a
          // node at each end of the mortar that is shared between the master
          // and dependent elements. Where this is not the case, the vertex
          // matching condition must be enforced by constraining the value of
          // the unknown at the node on the dependent side to be the same as the
          // value at that location in the master.

          // Store positions of the mortar vertex/non-vertex nodes in the
          // dependent element
          const unsigned n_mortar_vertices = 2;
          Vector<unsigned> vertex_pos(n_mortar_vertices);
          vertex_pos[0] = 0;
          vertex_pos[1] = this->ninterpolating_node_1d(value_id) - 1;
          Vector<unsigned> non_vertex_pos(my_n_p - n_mortar_vertices);
          for (unsigned i = 0; i < my_n_p - n_mortar_vertices; i++)
          {
            non_vertex_pos[i] = i + 1;
          }

          // If the node is on a master edge, we may be setting the
          // hanging info incorrectly. Hanging schemes (if they are
          // required) for such nodes are instead computed by the
          // dependent edge. We dont't want to overwrite them here!
          std::vector<bool> mortar_vertex_on_master_edge(n_mortar_vertices);
          for (unsigned v = 0; v < n_mortar_vertices; v++)
          {
            // Check if each node is on the master edge
            mortar_vertex_on_master_edge[v] = true;
            unsigned non_extreme_coordinate = 0;
            Vector<double> s_in_neigh(3);
            for (unsigned i = 0; i < 3; i++)
            {
              // Work out this node's location in the master
              s_in_neigh[i] =
                edge_s_lo_neigh[i] +
                dependent_node_s_fraction[vertex_pos[v]][edge_translate_s[i]] *
                  (edge_s_hi_neigh[i] - edge_s_lo_neigh[i]);

              // Check if local coordinate in master element takes non-extreme
              // value
              if (std::fabs(std::fabs(s_in_neigh[i]) - 1.0) > 1.0e-14)
              {
                non_extreme_coordinate++;
                if (non_extreme_coordinate > 1)
                {
                  mortar_vertex_on_master_edge[v] = false;
                  break;
                }
              }
            }
          }

          // Now work out if my edge coincides with the master edge
          bool my_edge_coincides_with_master = true;
          for (unsigned v = 0; v < n_mortar_vertices; v++)
          {
            my_edge_coincides_with_master =
              my_edge_coincides_with_master && mortar_vertex_on_master_edge[v];
          }

          // Check if we need to apply the (1D) vertex matching condition at the
          // mortar vertices. This is trivially satisfied if the node is shared
          // with the master, but if not then we need to constrain it.
          for (unsigned v = 0; v < n_mortar_vertices; v++)
          {
            // Don't make hanging node if my edge doesn't coincide with
            // the master edge *and* this node is on the master edge!
            // (We are not in a position to determine its hanging status.)
            if ((!my_edge_coincides_with_master) &&
                mortar_vertex_on_master_edge[v])
              continue;

            // Search master node storage for the node
            bool node_is_shared = false;
            for (unsigned i = 0; i < master_node_pt.size(); i++)
            {
              if (dependent_node_pt[vertex_pos[v]] == master_node_pt[i])
              {
                node_is_shared = true;
                break;
              }
            }

            // If the node is not shared then we must constrain its value by
            // setting up a hanging scheme
            if (!node_is_shared)
            {
              // Calculate weights. These are just the master shapes evaluated
              // at this dependent node's position

              // Work out this node's location in the master
              Vector<double> s_in_neigh(3);
              for (unsigned i = 0; i < 3; i++)
              {
                s_in_neigh[i] = edge_s_lo_neigh[i] +
                                dependent_node_s_fraction[vertex_pos[v]]
                                                         [edge_translate_s[i]] *
                                  (edge_s_hi_neigh[i] - edge_s_lo_neigh[i]);
              }

              // Get master shapes at dependent nodal positions
              edge_neigh_obj_pt->interpolating_basis(
                s_in_neigh, master_shapes, value_id);

              // Remove any existing hanging node info
              // (This may be a bit wasteful, but guarantees correctness)
              dependent_node_pt[vertex_pos[v]]->set_nonhanging();

              // Don't include master nodes with zero weights
              Vector<unsigned> master_node_zero_weight;
              for (unsigned m = 0; m < n_master_nodes; m++)
              {
                // Compare weights to some (small) tolerance
                if (std::fabs(master_shapes[master_node_number[m]]) < 1.0e-14)
                {
                  // Store
                  master_node_zero_weight.push_back(m);
                }
              }

              // Set up hanging scheme for this node
              HangInfo* explicit_hang_pt =
                new HangInfo(n_master_nodes - master_node_zero_weight.size());
              unsigned mindex = 0;
              for (unsigned m = 0; m < n_master_nodes; m++)
              {
                // Check that master doesn't have zero weight
                bool skip = false;
                for (unsigned i = 0; i < master_node_zero_weight.size(); i++)
                {
                  if (m == master_node_zero_weight[i]) skip = true;
                }

                // Add pointer and weight to hang info
                if (!skip)
                  explicit_hang_pt->set_master_node_pt(
                    mindex++,
                    master_node_pt[m],
                    master_shapes[master_node_number[m]]);
              }

              /// / Set up hanging scheme for this node
              // HangInfo* explicit_hang_pt = new HangInfo(n_master_nodes);
              // for(unsigned m=0; m<n_master_nodes; m++)
              // {
              //  explicit_hang_pt->set_master_node_pt(m,master_node_pt[m],master_shapes[master_node_number[m]]);
              // }

              // Make node hang
              dependent_node_pt[vertex_pos[v]]->set_hanging_pt(explicit_hang_pt,
                                                               -1);

              /// / Print out hanging scheme
              // std::cout << "Hanging node: "
              //          << dependent_node_pt[vertex_pos[v]]->x(0) << "  "
              //          << dependent_node_pt[vertex_pos[v]]->x(1) << "  "
              //          << dependent_node_pt[vertex_pos[v]]->x(2) << "  "
              //          << std::endl;
              // for(unsigned m=0; m<explicit_hang_pt->nmaster(); m++)
              // {
              //  std::cout << "   m = " << m << ": "
              //            << explicit_hang_pt->master_node_pt(m)->x(0) << "  "
              //            << explicit_hang_pt->master_node_pt(m)->x(1) << "  "
              //            << explicit_hang_pt->master_node_pt(m)->x(2) << "  "
              //            << "w = " << explicit_hang_pt->master_weight(m) << "
              //            "
              //            << std::endl;
              // }

              // Check there are no zero weights at this stage
              for (unsigned m = 0; m < explicit_hang_pt->nmaster(); m++)
              {
                if (std::fabs(explicit_hang_pt->master_weight(m)) < 1.0e-14)
                {
                  throw OomphLibError(
                    "Master has zero weight!",
                    "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
                    OOMPH_EXCEPTION_LOCATION);
                }
              }
            }
          }

          // Check that there are actually dependent nodes for which we still
          // need to construct a hanging scheme. If not then there is nothing
          // more to do.
          if (n_dependent_nodes - n_mortar_vertices > 0)
          {
            // Assemble mass matrix for mortar
            //--------------------------------
            Vector<double> psi(n_dependent_nodes - n_mortar_vertices);
            Vector<double> diag_M(n_dependent_nodes - n_mortar_vertices);
            Vector<Vector<double>> shared_node_M(n_mortar_vertices);
            for (unsigned i = 0; i < shared_node_M.size(); i++)
            {
              shared_node_M[i].resize(n_dependent_nodes - n_mortar_vertices);
            }

            // Fill in part corresponding to dependent nodal positions (unknown)
            for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
            {
              // Use L'Hosptal's rule:
              psi[i] =
                pow(-1.0, int((dependent_element_nnode_1d - 1) - i - 1)) *
                -Orthpoly::ddlegendre(
                  dependent_element_nnode_1d - 1,
                  dependent_nodal_position[non_vertex_pos[i]]);
              // Put in contribution
              diag_M[i] = psi[i] * dependent_weight[non_vertex_pos[i]];
            }

            // Fill in part corresponding to dependent element vertices (known)
            for (unsigned v = 0; v < shared_node_M.size(); v++)
            {
              for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices;
                   i++)
              {
                // Check if denominator is zero
                if (std::fabs(dependent_nodal_position[non_vertex_pos[i]] -
                              dependent_nodal_position[vertex_pos[v]]) >=
                    1.0e-8)
                {
                  // We're ok
                  psi[i] =
                    pow(-1.0, int((dependent_element_nnode_1d - 1) - i - 1)) *
                    Orthpoly::dlegendre(
                      dependent_element_nnode_1d - 1,
                      dependent_nodal_position[vertex_pos[v]]) /
                    (dependent_nodal_position[non_vertex_pos[i]] -
                     dependent_nodal_position[vertex_pos[v]]);
                }
                // Check if numerator is zero
                else if (std::fabs(Orthpoly::dlegendre(
                           dependent_element_nnode_1d - 1,
                           dependent_nodal_position[vertex_pos[v]])) < 1.0e-8)
                {
                  // We can use l'hopital's rule
                  psi[i] =
                    pow(-1.0, int((dependent_element_nnode_1d - 1) - i - 1)) *
                    -Orthpoly::ddlegendre(
                      dependent_element_nnode_1d - 1,
                      dependent_nodal_position[non_vertex_pos[i]]);
                }
                else
                {
                  // We can't use l'hopital's rule
                  throw OomphLibError(
                    "Cannot use l'Hopital's rule. Dividing by zero is not "
                    "allowed!",
                    "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
                    OOMPH_EXCEPTION_LOCATION);
                }
                // Put in contribution
                shared_node_M[v][i] = psi[i] * dependent_weight[vertex_pos[v]];
              }
            }

            // Assemble local projection matrix for mortar
            //--------------------------------------------

            // Have only one local projection matrix because there is just one
            // master
            Vector<Vector<double>> P(n_dependent_nodes - n_mortar_vertices);
            for (unsigned i = 0; i < P.size(); i++)
            {
              P[i].resize(n_master_nodes, 0.0);
            }

            // Storage for local coordinate
            Vector<double> s(3);

            // Sum contributions from master element shapes (quadrature).
            // The order of quadrature must be high enough to evaluate a
            // polynomial of order N_s + N_m - 3 exactly, where N_s =
            // n_dependent_nodes, N_m = n_master_nodes. (Use pointers for the
            // quadrature knots and weights so that data is not unnecessarily
            // copied)
            // unsigned quadrature_order =
            // std::max(dependent_element_nnode_1d,master_element_nnode_1d);
            Vector<double>*quadrature_knot, *quadrature_weight;
            if (dependent_element_nnode_1d >= master_element_nnode_1d)
            {
              // Use the same quadrature order as the dependent element (me)
              quadrature_knot = &dependent_nodal_position;
              quadrature_weight = &dependent_weight;
            }
            else
            {
              // Use the same quadrature order as the master element (neighbour)
              quadrature_knot = &master_nodal_position;
              quadrature_weight = &master_weight;
            }

            // Quadrature loop
            for (unsigned q = 0; q < (*quadrature_weight).size(); q++)
            {
              // Evaluate mortar test functions at the quadrature knots in the
              // dependent
              // s[active_coord_index] = (*quadrature_knot)[q];
              Vector<double> s_on_mortar(1);
              s_on_mortar[0] = (*quadrature_knot)[q];

              // Get psi
              for (unsigned k = 0; k < n_dependent_nodes - n_mortar_vertices;
                   k++)
              {
                // Check if denominator is zero
                if (std::fabs(dependent_nodal_position[non_vertex_pos[k]] -
                              s_on_mortar[0]) >= 1.0e-08)
                {
                  // We're ok
                  psi[k] =
                    pow(-1.0, int((dependent_element_nnode_1d - 1) - k - 1)) *
                    Orthpoly::dlegendre(dependent_element_nnode_1d - 1,
                                        s_on_mortar[0]) /
                    (dependent_nodal_position[non_vertex_pos[k]] -
                     s_on_mortar[0]);
                }
                // Check if numerator is zero
                else if (std::fabs(Orthpoly::dlegendre(
                           dependent_element_nnode_1d - 1, s_on_mortar[0])) <
                         1.0e-8)
                {
                  // We can use l'Hopital's rule
                  psi[k] =
                    pow(-1.0, int((dependent_element_nnode_1d - 1) - k - 1)) *
                    -Orthpoly::ddlegendre(dependent_element_nnode_1d - 1,
                                          s_on_mortar[0]);
                }
                else
                {
                  // We can't use l'hopital's rule
                  throw OomphLibError(
                    "Cannot use l'Hopital's rule. Dividing by zero is not "
                    "allowed!",
                    "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
                    OOMPH_EXCEPTION_LOCATION);
                }
              }

              // Convert coordinate on mortar to local fraction in dependent
              // element
              Vector<double> s_fraction(3);
              for (unsigned i = 0; i < 3; i++)
              {
                s_fraction[i] = (i == active_coord_index) ?
                                  0.5 * (s_on_mortar[0] + 1.0) :
                                  dependent_node_s_fraction[vertex_pos[0]][i];
              }

              // Project active coordinate into master element
              Vector<double> s_in_neigh(3);
              for (unsigned i = 0; i < 3; i++)
              {
                s_in_neigh[i] = edge_s_lo_neigh[i] +
                                s_fraction[edge_translate_s[i]] *
                                  (edge_s_hi_neigh[i] - edge_s_lo_neigh[i]);
              }

              // Evaluate master shapes at projections of local quadrature knots
              edge_neigh_obj_pt->interpolating_basis(
                s_in_neigh, master_shapes, value_id);

              // Populate local projection matrix
              for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices;
                   i++)
              {
                for (unsigned j = 0; j < n_master_nodes; j++)
                {
                  P[i][j] += master_shapes[master_node_number[j]] * psi[i] *
                             (*quadrature_weight)[q];
                }
              }
            }

            /// / Print out local projection matrix
            // std::cout << "P_e:" << std::endl;
            // for(unsigned i=0; i<P.size(); i++)
            // {
            //  for(unsigned j=0; j<P[i].size(); j++)
            //   {
            //    std::cout << "  " << P[i][j];
            //   }
            // }
            // std::cout << std::endl;

            // Assemble global projection matrices for mortar
            //-----------------------------------------------
            // Need to subtract contributions from the "known unknowns"
            // corresponding to the nodes at the vertices of the mortar

            // Assemble contributions from mortar vertex nodes
            for (unsigned v = 0; v < n_mortar_vertices; v++)
            {
              // Convert coordinate on mortar to local fraction in dependent
              // element
              Vector<double> s_fraction(3);
              for (unsigned i = 0; i < 3; i++)
              {
                s_fraction[i] =
                  (i == active_coord_index) ?
                    0.5 * (dependent_nodal_position[vertex_pos[v]] + 1.0) :
                    dependent_node_s_fraction[vertex_pos[0]][i];
              }

              // Project active coordinate into master element
              Vector<double> s_in_neigh(3);
              for (unsigned i = 0; i < 3; i++)
              {
                s_in_neigh[i] = edge_s_lo_neigh[i] +
                                s_fraction[edge_translate_s[i]] *
                                  (edge_s_hi_neigh[i] - edge_s_lo_neigh[i]);
              }

              // Get master shapes at dependent nodal positions
              edge_neigh_obj_pt->interpolating_basis(
                s_in_neigh, master_shapes, value_id);

              for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices;
                   i++)
              {
                for (unsigned k = 0; k < n_master_nodes; k++)
                {
                  P[i][k] -=
                    master_shapes[master_node_number[k]] * shared_node_M[v][i];
                }
              }
            }

            /// / Print out global projection matrix
            // std::cout << "P:" << std::endl;
            // for(unsigned i=0; i<P.size(); i++)
            // {
            //  for(unsigned j=0; j<P[i].size(); j++)
            //   {
            //    std::cout << "  " << P[i][j];
            //   }
            //  std::cout << std::endl;
            // }
            // std::cout << std::endl;

            // Solve mortar system
            //--------------------
            for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
            {
              for (unsigned j = 0; j < n_master_nodes; j++)
              {
                P[i][j] /= diag_M[i];
              }
            }

            /// / Print out solved global projection matrix
            // std::cout << "solved P:" << std::endl;
            // for(unsigned i=0; i<P.size(); i++)
            // {
            //  for(unsigned j=0; j<P[i].size(); j++)
            //   {
            //    std::cout << "  " << P[i][j];
            //   }
            //  std::cout << std::endl;
            // }
            // std::cout << std::endl;

            // Create and populate structures to hold the hanging info
            //--------------------------------------------------------
            Vector<HangInfo*> hang_info_pt(n_dependent_nodes);
            for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
            {
              // Don't include master nodes with zero weights
              Vector<unsigned> master_node_zero_weight;
              for (unsigned m = 0; m < n_master_nodes; m++)
              {
                // Compare weights to some (small) tolerance
                if (std::fabs(P[i][m]) < 1.0e-14)
                {
                  // Store
                  master_node_zero_weight.push_back(m);
                }
              }

              // Set up hanging scheme for this node
              hang_info_pt[i] =
                new HangInfo(n_master_nodes - master_node_zero_weight.size());
              unsigned mindex = 0;
              for (unsigned m = 0; m < n_master_nodes; m++)
              {
                // Check that master doesn't have zero weight
                bool skip = false;
                for (unsigned k = 0; k < master_node_zero_weight.size(); k++)
                {
                  if (m == master_node_zero_weight[k]) skip = true;
                }

                // Add pointer and weight to hang info
                if (!skip)
                  hang_info_pt[i]->set_master_node_pt(
                    mindex++, master_node_pt[m], P[i][m]);
              }
            }

            /// / Create structures to hold the hanging info
            /// /-------------------------------------------
            // Vector<HangInfo*> hang_info_pt(n_dependent_nodes);
            // for (unsigned i=0; i<n_dependent_nodes-n_mortar_vertices; i++)
            // {
            //  hang_info_pt[i] = new HangInfo(n_master_nodes);
            // }
            //
            /// / Copy information to hanging nodes
            /// /----------------------------------
            // for(unsigned i=0; i<n_dependent_nodes-n_mortar_vertices; i++)
            // {
            //  for(unsigned j=0; j<n_master_nodes; j++)
            //   {
            //    hang_info_pt[i]->set_master_node_pt(j,master_node_pt[j],P[i][j]);
            //   }
            // }

            // Set pointers to hanging info
            //-----------------------------
            for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
            {
              // Check that the node shoule actually hang.
              // That is, if the polynomial orders of the elements at a p-type
              // non-conormity are both odd then the middle node on the edge is
              // shared but a hanging scheme has been constructed for it.
              bool node_is_really_shared = false;
              for (unsigned m = 0; m < hang_info_pt[i]->nmaster(); m++)
              {
                // We can simply check if the hanging scheme lists itself as a
                // master
                if (hang_info_pt[i]->master_node_pt(m) ==
                    dependent_node_pt[non_vertex_pos[i]])
                {
                  node_is_really_shared = true;

#ifdef PARANOID
                  // Also check the corresponding weight: it should be 1
                  if (std::fabs(hang_info_pt[i]->master_weight(m) - 1.0) >
                      1.0e-06)
                  {
                    throw OomphLibError("Something fishy here -- with shared "
                                        "node at a mortar vertex",
                                        "PRefineableQElemen<2,INITIAL_NNODE_1D>"
                                        "t::quad_hang_helper()",
                                        OOMPH_EXCEPTION_LOCATION);
                  }
#endif
                }
              }

              // Now we can make the node hang, if it isn't a shared node
              if (!node_is_really_shared)
              {
                dependent_node_pt[non_vertex_pos[i]]->set_hanging_pt(
                  hang_info_pt[i], -1);

                /// / Print out hanging scheme
                // std::cout << "Hanging node: "
                //          << dependent_node_pt[non_vertex_pos[i]]->x(0) << " "
                //          << dependent_node_pt[non_vertex_pos[i]]->x(1) << " "
                //          << dependent_node_pt[non_vertex_pos[i]]->x(2) << " "
                //          << std::endl;
                // for(unsigned m=0; m<hang_info_pt[i]->nmaster(); m++)
                // {
                //  std::cout << "   m = " << m << ": "
                //            << hang_info_pt[i]->master_node_pt(m)->x(0) << " "
                //            << hang_info_pt[i]->master_node_pt(m)->x(1) << " "
                //            << hang_info_pt[i]->master_node_pt(m)->x(2) << " "
                //            << "w = " << hang_info_pt[i]->master_weight(m) <<
                //            "  "
                //            << std::endl;
                // }
              }
            }

          } // End of case where there are still dependent nodes

        } // End of edge is dependent

      } // End of loop over face edges

    } // End of if face is dependent

    // Now do the mortaring
    //---------------------
    if (h_type_dependent || p_type_dependent)
    {
      // Compute the active coordinate indices along the this side of mortar
      Vector<unsigned> active_coord_index(2);
      if (my_face == B || my_face == F)
      {
        active_coord_index[0] = 0;
        active_coord_index[1] = 1;
      }
      else if (my_face == D || my_face == U)
      {
        active_coord_index[0] = 0;
        active_coord_index[1] = 2;
      }
      else if (my_face == L || my_face == R)
      {
        active_coord_index[0] = 1;
        active_coord_index[1] = 2;
      }
      else
      {
        throw OomphLibError(
          "Cannot transform coordinates",
          "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
          OOMPH_EXCEPTION_LOCATION);
      }

      // Get pointer to neighbouring master element (in p-refineable form)
      PRefineableQElement<3, INITIAL_NNODE_1D>* neigh_obj_pt;
      neigh_obj_pt = dynamic_cast<PRefineableQElement<3, INITIAL_NNODE_1D>*>(
        neigh_pt->object_pt());

      // Create vector of master and dependent nodes
      //----------------------------------------
      Vector<Node*> master_node_pt, dependent_node_pt;
      Vector<unsigned> master_node_number, dependent_node_number;
      Vector<Vector<double>> dependent_node_s_fraction;

      // Number of nodes in one dimension
      const unsigned my_n_p = this->ninterpolating_node_1d(value_id);
      const unsigned neigh_n_p = neigh_obj_pt->ninterpolating_node_1d(value_id);

      // Storage for pointers to the nodes and their numbers along the master
      // edge
      unsigned neighbour_node_number = 0;
      Node* neighbour_node_pt = 0;

      // Loop over nodes on the face
      for (unsigned i0 = 0; i0 < neigh_n_p; i0++)
      {
        for (unsigned i1 = 0; i1 < neigh_n_p; i1++)
        {
          const unsigned s0space = 1;
          const unsigned s1space = neigh_n_p;
          const unsigned s2space = neigh_n_p * neigh_n_p;

          // Find the neighbour's node
          switch (neigh_face)
          {
            case B:
              neighbour_node_number = i0 * s0space + i1 * s1space;
              neighbour_node_pt = neigh_obj_pt->interpolating_node_pt(
                neighbour_node_number, value_id);
              break;

            case F:
              neighbour_node_number =
                (neigh_n_p - 1) * s2space + i0 * s0space + i1 * s1space;
              neighbour_node_pt = neigh_obj_pt->interpolating_node_pt(
                neighbour_node_number, value_id);
              break;

            case D:
              neighbour_node_number = i0 * s0space + i1 * s2space;
              neighbour_node_pt = neigh_obj_pt->interpolating_node_pt(
                neighbour_node_number, value_id);
              break;

            case U:
              neighbour_node_number =
                (neigh_n_p - 1) * s1space + i0 * s0space + i1 * s2space;
              neighbour_node_pt = neigh_obj_pt->interpolating_node_pt(
                neighbour_node_number, value_id);
              break;

            case L:
              neighbour_node_number = i0 * s1space + i1 * s2space;
              neighbour_node_pt = neigh_obj_pt->interpolating_node_pt(
                neighbour_node_number, value_id);
              break;

            case R:
              neighbour_node_number =
                (neigh_n_p - 1) * s0space + i0 * s1space + i1 * s2space;
              neighbour_node_pt = neigh_obj_pt->interpolating_node_pt(
                neighbour_node_number, value_id);
              break;

            default:
              throw OomphLibError(
                "my_face not L, R, D, U, B, F\n",
                "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
                OOMPH_EXCEPTION_LOCATION);
          }

          // Set node as master node
          master_node_number.push_back(neighbour_node_number);
          master_node_pt.push_back(neighbour_node_pt);
        }
      }

      // Storage for pointers to the local nodes and their numbers along my edge
      unsigned local_node_number = 0;
      Node* local_node_pt = 0;

      // Loop over the nodes along my edge
      for (unsigned i0 = 0; i0 < my_n_p; i0++)
      {
        for (unsigned i1 = 0; i1 < my_n_p; i1++)
        {
          // Storage for the fractional position of the node in the element
          Vector<double> s_fraction(3);

          const unsigned s0space = 1;
          const unsigned s1space = my_n_p;
          const unsigned s2space = my_n_p * my_n_p;

          // Find the local node and the fractional position of the node
          // which depends on the edge, of course
          switch (my_face)
          {
            case B:
              s_fraction[0] =
                local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
              s_fraction[1] =
                local_one_d_fraction_of_interpolating_node(i1, 1, value_id);
              s_fraction[2] = 0.0;
              local_node_number = i0 * s0space + i1 * s1space;
              local_node_pt =
                this->interpolating_node_pt(local_node_number, value_id);
              break;

            case F:
              s_fraction[0] =
                local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
              s_fraction[1] =
                local_one_d_fraction_of_interpolating_node(i1, 1, value_id);
              s_fraction[2] = 1.0;
              local_node_number =
                (my_n_p - 1) * s2space + i0 * s0space + i1 * s1space;
              local_node_pt =
                this->interpolating_node_pt(local_node_number, value_id);
              break;

            case D:
              s_fraction[0] =
                local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
              s_fraction[1] = 0.0;
              s_fraction[2] =
                local_one_d_fraction_of_interpolating_node(i1, 2, value_id);
              local_node_number = i0 * s0space + i1 * s2space;
              local_node_pt =
                this->interpolating_node_pt(local_node_number, value_id);
              break;

            case U:
              s_fraction[0] =
                local_one_d_fraction_of_interpolating_node(i0, 0, value_id);
              s_fraction[1] = 1.0;
              s_fraction[2] =
                local_one_d_fraction_of_interpolating_node(i1, 2, value_id);
              local_node_number =
                (my_n_p - 1) * s1space + i0 * s0space + i1 * s2space;
              local_node_pt =
                this->interpolating_node_pt(local_node_number, value_id);
              break;

            case L:
              s_fraction[0] = 0.0;
              s_fraction[1] =
                local_one_d_fraction_of_interpolating_node(i0, 1, value_id);
              s_fraction[2] =
                local_one_d_fraction_of_interpolating_node(i1, 2, value_id);
              local_node_number = i0 * s1space + i1 * s2space;
              local_node_pt =
                this->interpolating_node_pt(local_node_number, value_id);
              break;

            case R:
              s_fraction[0] = 1.0;
              s_fraction[1] =
                local_one_d_fraction_of_interpolating_node(i0, 1, value_id);
              s_fraction[2] =
                local_one_d_fraction_of_interpolating_node(i1, 2, value_id);
              local_node_number =
                (my_n_p - 1) * s0space + i0 * s1space + i1 * s2space;
              local_node_pt =
                this->interpolating_node_pt(local_node_number, value_id);
              break;

            default:
              throw OomphLibError(
                "my_face not L, R, D, U, B, F\n",
                "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
                OOMPH_EXCEPTION_LOCATION);
          }

          // Add node to vector of dependent element nodes
          dependent_node_number.push_back(local_node_number);
          dependent_node_pt.push_back(local_node_pt);

          // Store node's local fraction
          dependent_node_s_fraction.push_back(s_fraction);
        }
      }

      // Store the number of dependent and master nodes for use later
      const unsigned n_dependent_nodes = dependent_node_pt.size();
      const unsigned n_master_nodes = master_node_pt.size();
      const unsigned dependent_element_nnode_1d = my_n_p;
      const unsigned master_element_nnode_1d = neigh_n_p;

      /// / Print out dependent and master node coords
      // std::cout << "Dependent nodes on face: " <<
      // OcTree::Direct_string[my_face] << std::endl; for(unsigned i=0;
      // i<dependent_node_pt.size(); i++)
      // {
      //  std::cout << i << ": "
      //            << dependent_node_pt[i]->x(0) << "  "
      //            << dependent_node_pt[i]->x(1) << "  "
      //            << dependent_node_pt[i]->x(2) << "  "
      //            << std::endl;
      // }
      // std::cout << "Master nodes on face: " << OcTree::Direct_string[my_face]
      // << std::endl; for(unsigned i=0; i<master_node_pt.size(); i++)
      // {
      //  std::cout << i << ": "
      //            << master_node_pt[i]->x(0) << "  "
      //            << master_node_pt[i]->x(1) << "  "
      //            << master_node_pt[i]->x(2) << "  "
      //            << std::endl;
      // }

      // Storage for master shapes
      Shape master_shapes(neigh_obj_pt->ninterpolating_node(value_id));

      // Get master and dependent nodal positions
      //-------------------------------------
      // Get in 1D
      Vector<double> dependent_nodal_position_1d;
      Vector<double> dependent_weight_1d(dependent_element_nnode_1d);
      Orthpoly::gll_nodes(dependent_element_nnode_1d,
                          dependent_nodal_position_1d,
                          dependent_weight_1d);
      Vector<double> master_nodal_position_1d;
      Vector<double> master_weight_1d(master_element_nnode_1d);
      Orthpoly::gll_nodes(
        master_element_nnode_1d, master_nodal_position_1d, master_weight_1d);

      // Storage for 2D
      Vector<Vector<double>> dependent_nodal_position(
        dependent_element_nnode_1d * dependent_element_nnode_1d);
      for (unsigned i = 0; i < dependent_nodal_position.size(); i++)
      {
        dependent_nodal_position[i].resize(2);
      }
      Vector<double> dependent_weight(dependent_element_nnode_1d *
                                      dependent_element_nnode_1d);
      Vector<Vector<double>> master_nodal_position(master_element_nnode_1d *
                                                   master_element_nnode_1d);
      for (unsigned i = 0; i < master_nodal_position.size(); i++)
      {
        master_nodal_position[i].resize(2);
      }
      Vector<double> master_weight(master_element_nnode_1d *
                                   master_element_nnode_1d);

      // Fill in coordinates and weights in 2D
      unsigned dependent_index = 0;
      for (unsigned i = 0; i < dependent_element_nnode_1d; i++)
      {
        for (unsigned j = 0; j < dependent_element_nnode_1d; j++)
        {
          dependent_nodal_position[dependent_index][0] =
            dependent_nodal_position_1d[i];
          dependent_nodal_position[dependent_index][1] =
            dependent_nodal_position_1d[j];
          dependent_weight[dependent_index] =
            dependent_weight_1d[i] * dependent_weight_1d[j];
          dependent_index++;
        }
      }
      unsigned master_index = 0;
      for (unsigned i = 0; i < master_element_nnode_1d; i++)
      {
        for (unsigned j = 0; j < master_element_nnode_1d; j++)
        {
          master_nodal_position[master_index][0] = master_nodal_position_1d[i];
          master_nodal_position[master_index][1] = master_nodal_position_1d[j];
          master_weight[master_index] =
            master_weight_1d[i] * master_weight_1d[j];
          master_index++;
        }
      }

      // Apply the vertex matching condition
      //------------------------------------
      // Vertiex matching is ensured automatically in cases where there is a
      // node at each end of the mortar that is shared between the master and
      // dependent elements. Where this is not the case, the vertex matching
      // condition must be enforced by constraining the value of the unknown at
      // the node on the dependent side to be the same as the value at that
      // location in the master.

      // Store positions of the mortar vertex/non-vertex nodes in the dependent
      // element
      // const unsigned n_mortar_vertices = 4;
      // Vector<unsigned> vertex_pos(n_mortar_vertices);
      // vertex_pos[0] = 0;
      // vertex_pos[1] = my_n_p-1;
      // vertex_pos[2] = my_n_p*(my_n_p-1);
      // vertex_pos[3] = my_n_p*my_n_p-1;
      Vector<unsigned> non_vertex_pos((dependent_element_nnode_1d - 2) *
                                      (dependent_element_nnode_1d - 2));
      unsigned nvindex = 0;
      for (unsigned i = 1; i < dependent_element_nnode_1d - 1; i++)
      {
        for (unsigned j = 1; j < dependent_element_nnode_1d - 1; j++)
        {
          non_vertex_pos[nvindex++] = i * dependent_element_nnode_1d + j;
        }
      }
      Vector<unsigned> vertex_pos;
      for (unsigned i = 0; i < n_dependent_nodes; i++)
      {
        // Check if node number is in the non-vertex storage
        bool node_is_vertex = true;
        for (unsigned j = 0; j < non_vertex_pos.size(); j++)
        {
          if (i == non_vertex_pos[j])
          {
            // Node is not a vertex
            node_is_vertex = false;
            break;
          }
        }
        // If we get here and the node is a vertex then add it's index
        if (node_is_vertex)
        {
          vertex_pos.push_back(i);
        }
      }
      // Store number of mortar vertices
      const unsigned n_mortar_vertices = vertex_pos.size();

      // Check that there are actually dependent nodes for which we still need
      // to construct a hanging scheme. If not then there is nothing more to do.
      if (n_dependent_nodes - n_mortar_vertices > 0)
      {
        // Assemble mass matrix for mortar
        //--------------------------------
        Vector<double> psi(n_dependent_nodes - n_mortar_vertices);
        Vector<double> diag_M(n_dependent_nodes - n_mortar_vertices);
        Vector<Vector<double>> shared_node_M(n_mortar_vertices);
        for (unsigned i = 0; i < shared_node_M.size(); i++)
        {
          shared_node_M[i].resize(n_dependent_nodes - n_mortar_vertices);
        }

        // Fill in part corresponding to dependent nodal positions (unknown)
        for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
        {
          // Mortar test functions in 2D are just the cross product of the 1D
          // test functions Initialise to 1
          psi[i] = 1.0;
          // Take product in each direction
          for (unsigned dir = 0; dir < 2; dir++)
          {
            unsigned index1d =
              (dir == 0) ? i : i % (dependent_element_nnode_1d - 2);
            // Use L'Hosptal's rule:
            psi[i] *=
              pow(-1.0, int((dependent_element_nnode_1d - 1) - index1d - 1)) *
              -Orthpoly::ddlegendre(
                dependent_element_nnode_1d - 1,
                dependent_nodal_position[non_vertex_pos[i]][dir]);
          }
          // Put in contribution
          diag_M[i] = psi[i] * dependent_weight[non_vertex_pos[i]];
        }

        /// / Print out diag(M)
        // std::cout << "diag(M):" << std::endl;
        // for(unsigned i=0; i<diag_M.size(); i++)
        // {
        //  std::cout << "  " << diag_M[i];
        // }
        // std::cout << std::endl;

        // Fill in part corresponding to dependent element vertices (known)
        for (unsigned v = 0; v < shared_node_M.size(); v++)
        {
          for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
          {
            // Mortar test functions in 2D are just the cross product of the 1D
            // test functions Initialise to 1
            psi[i] = 1.0;
            // Take product in each direction
            for (unsigned dir = 0; dir < 2; dir++)
            {
              unsigned index1d =
                (dir == 0) ? i : i % (dependent_element_nnode_1d - 2);
              // Check if denominator is zero
              if (std::fabs(dependent_nodal_position[non_vertex_pos[i]][dir] -
                            dependent_nodal_position[vertex_pos[v]][dir]) >=
                  1.0e-8)
              {
                // We're ok
                psi[i] *=
                  pow(-1.0,
                      int((dependent_element_nnode_1d - 1) - index1d - 1)) *
                  Orthpoly::dlegendre(
                    dependent_element_nnode_1d - 1,
                    dependent_nodal_position[vertex_pos[v]][dir]) /
                  (dependent_nodal_position[non_vertex_pos[i]][dir] -
                   dependent_nodal_position[vertex_pos[v]][dir]);
              }
              // Check if numerator is zero
              else if (std::fabs(Orthpoly::dlegendre(
                         dependent_element_nnode_1d - 1,
                         dependent_nodal_position[vertex_pos[v]][dir])) <
                       1.0e-8)
              {
                // We can use l'hopital's rule
                psi[i] *=
                  pow(-1.0,
                      int((dependent_element_nnode_1d - 1) - index1d - 1)) *
                  -Orthpoly::ddlegendre(
                    dependent_element_nnode_1d - 1,
                    dependent_nodal_position[non_vertex_pos[i]][dir]);
              }
              else
              {
                // We can't use l'hopital's rule
                throw OomphLibError(
                  "Cannot use l'Hopital's rule. Dividing by zero is not "
                  "allowed!",
                  "PRefineableQElement<3,INITIAL_NNODE_1D>::quad_hang_helper()",
                  OOMPH_EXCEPTION_LOCATION);
              }
            }
            // Put in contribution
            shared_node_M[v][i] = psi[i] * dependent_weight[vertex_pos[v]];
          }
        }

        /// / Print out diag(M)
        // std::cout << "shared node M:" << std::endl;
        // for(unsigned i=0; i<shared_node_M.size(); i++)
        // {
        //  for(unsigned j=0; j<shared_node_M[i].size(); j++)
        //   {
        //    std::cout << "  " << shared_node_M[i][j];
        //   }
        // }
        // std::cout << std::endl;

        // Assemble local projection matrix for mortar
        //--------------------------------------------

        // Have only one local projection matrix because there is just one
        // master
        Vector<Vector<double>> P(n_dependent_nodes - n_mortar_vertices);
        for (unsigned i = 0; i < P.size(); i++)
        {
          P[i].resize(n_master_nodes, 0.0);
        }

        // Storage for local coordinate
        Vector<double> s(3);

        // Sum contributions from master element shapes (quadrature).
        // The order of quadrature must be high enough to evaluate a polynomial
        // of order N_s + N_m - 3 exactly, where N_s = n_dependent_nodes, N_m =
        // n_master_nodes.
        // (Use pointers for the quadrature knots and weights so that
        // data is not unnecessarily copied)
        // unsigned quadrature_order =
        // std::max(dependent_element_nnode_1d,master_element_nnode_1d);
        Vector<Vector<double>>* quadrature_knot;
        Vector<double>* quadrature_weight;
        if (dependent_element_nnode_1d >= master_element_nnode_1d)
        {
          // Use the same quadrature order as the dependent element (me)
          quadrature_knot = &dependent_nodal_position;
          quadrature_weight = &dependent_weight;
        }
        else
        {
          // Use the same quadrature order as the master element (neighbour)
          quadrature_knot = &master_nodal_position;
          quadrature_weight = &master_weight;
        }

        // Quadrature loop
        for (unsigned q = 0; q < (*quadrature_weight).size(); q++)
        {
          // Evaluate mortar test functions at the quadrature knots in the
          // dependent
          Vector<double> s_on_mortar(2);
          for (unsigned i = 0; i < 2; i++)
          {
            s_on_mortar[i] = (*quadrature_knot)[q][i];
          }

          // Get psi
          for (unsigned k = 0; k < n_dependent_nodes - n_mortar_vertices; k++)
          {
            // Mortar test functions in 2D are just the cross product of the 1D
            // test functions Initialise to 1
            psi[k] = 1.0;
            // Take product in each direction
            for (unsigned dir = 0; dir < 2; dir++)
            {
              unsigned index1d =
                (dir == 0) ? k : k % (dependent_element_nnode_1d - 2);
              // Check if denominator is zero
              if (std::fabs(dependent_nodal_position[non_vertex_pos[k]][dir] -
                            s_on_mortar[dir]) >= 1.0e-08)
              {
                // We're ok
                psi[k] *=
                  pow(-1.0,
                      int((dependent_element_nnode_1d - 1) - index1d - 1)) *
                  Orthpoly::dlegendre(dependent_element_nnode_1d - 1,
                                      s_on_mortar[dir]) /
                  (dependent_nodal_position[non_vertex_pos[k]][dir] -
                   s_on_mortar[dir]);
              }
              // Check if numerator is zero
              else if (std::fabs(Orthpoly::dlegendre(
                         dependent_element_nnode_1d - 1, s_on_mortar[dir])) <
                       1.0e-8)
              {
                // We can use l'Hopital's rule
                psi[k] *=
                  pow(-1.0,
                      int((dependent_element_nnode_1d - 1) - index1d - 1)) *
                  -Orthpoly::ddlegendre(dependent_element_nnode_1d - 1,
                                        s_on_mortar[dir]);
              }
              else
              {
                // We can't use l'hopital's rule
                throw OomphLibError(
                  "Cannot use l'Hopital's rule. Dividing by zero is not "
                  "allowed!",
                  "PRefineableQElement<3,INITIAL_NNODE_1D>::quad_hang_helper()",
                  OOMPH_EXCEPTION_LOCATION);
              }
            }
          }

          // Convert coordinate on mortar to local fraction in dependent element
          Vector<double> s_fraction(3);
          for (unsigned i = 0; i < 3; i++)
          {
            if (i == active_coord_index[0])
            {
              s_fraction[i] = 0.5 * (s_on_mortar[0] + 1.0);
            }
            else if (i == active_coord_index[1])
            {
              s_fraction[i] = 0.5 * (s_on_mortar[1] + 1.0);
            }
            else
            {
              s_fraction[i] = dependent_node_s_fraction[vertex_pos[0]][i];
            }
          }

          // Project active coordinate into master element
          Vector<double> s_in_neigh(3);
          for (unsigned i = 0; i < 3; i++)
          {
            s_in_neigh[i] = s_lo_neigh[i] + s_fraction[translate_s[i]] *
                                              (s_hi_neigh[i] - s_lo_neigh[i]);
          }

          // Evaluate master shapes at projections of local quadrature knots
          neigh_obj_pt->interpolating_basis(
            s_in_neigh, master_shapes, value_id);

          // Populate local projection matrix
          for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
          {
            for (unsigned j = 0; j < n_master_nodes; j++)
            {
              P[i][j] += master_shapes[master_node_number[j]] * psi[i] *
                         (*quadrature_weight)[q];
            }
          }
        }

        /// / Print out local projection matrix
        // std::cout << "P_e:" << std::endl;
        // for(unsigned i=0; i<P.size(); i++)
        // {
        //  for(unsigned j=0; j<P[i].size(); j++)
        //   {
        //    std::cout << "  " << P[i][j];
        //   }
        // }
        // std::cout << std::endl;

        // Assemble global projection matrices for mortar
        //-----------------------------------------------
        // Need to subtract contributions from the "known unknowns"
        // corresponding to the nodes at the vertices of the mortar

        // Assemble contributions from mortar vertex nodes
        for (unsigned v = 0; v < n_mortar_vertices; v++)
        {
          // Convert coordinate on mortar to local fraction in dependent element
          Vector<double> s_fraction(3);
          for (unsigned i = 0; i < 3; i++)
          {
            if (i == active_coord_index[0])
            {
              s_fraction[i] =
                0.5 * (dependent_nodal_position[vertex_pos[v]][0] + 1.0);
            }
            else if (i == active_coord_index[1])
            {
              s_fraction[i] =
                0.5 * (dependent_nodal_position[vertex_pos[v]][1] + 1.0);
            }
            else
            {
              s_fraction[i] = dependent_node_s_fraction[vertex_pos[0]][i];
            }
          }

          // Project active coordinate into master element
          Vector<double> s_in_neigh(3);
          for (unsigned i = 0; i < 3; i++)
          {
            s_in_neigh[i] = s_lo_neigh[i] + s_fraction[translate_s[i]] *
                                              (s_hi_neigh[i] - s_lo_neigh[i]);
          }

          // Get master shapes at dependent nodal positions
          neigh_obj_pt->interpolating_basis(
            s_in_neigh, master_shapes, value_id);

          for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
          {
            for (unsigned k = 0; k < n_master_nodes; k++)
            {
              P[i][k] -=
                master_shapes[master_node_number[k]] * shared_node_M[v][i];
            }
          }
        }

        /// / Print out global projection matrix
        // std::cout << "P:" << std::endl;
        // for(unsigned i=0; i<P.size(); i++)
        // {
        //  for(unsigned j=0; j<P[i].size(); j++)
        //   {
        //    std::cout << "  " << P[i][j];
        //   }
        // }
        // std::cout << std::endl;

        // Solve mortar system
        //--------------------
        for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
        {
          for (unsigned j = 0; j < n_master_nodes; j++)
          {
            P[i][j] /= diag_M[i];
          }
        }

        /// / Print out solved matrix
        // std::cout << "solved P:" << std::endl;
        // for(unsigned i=0; i<P.size(); i++)
        // {
        //  for(unsigned j=0; j<P[i].size(); j++)
        //   {
        //    std::cout << "  " << P[i][j];
        //   }
        // }
        // std::cout << std::endl;

        // Create structures to hold the hanging info
        //-------------------------------------------
        Vector<HangInfo*> hang_info_pt(n_dependent_nodes);
        for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
        {
          hang_info_pt[i] = new HangInfo(n_master_nodes);
        }

        // Copy information to hanging nodes
        //----------------------------------
        for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
        {
          for (unsigned j = 0; j < n_master_nodes; j++)
          {
            hang_info_pt[i]->set_master_node_pt(j, master_node_pt[j], P[i][j]);
          }
        }

        // Set pointers to hanging info
        //-----------------------------
        for (unsigned i = 0; i < n_dependent_nodes - n_mortar_vertices; i++)
        {
          // Check that the node shoule actually hang.
          // That is, if the polynomial orders of the elements at a p-type
          // non-conormity are both odd then the middle node on the edge is
          // shared but a hanging scheme has been constructed for it.
          bool node_is_really_shared = false;
          for (unsigned m = 0; m < hang_info_pt[i]->nmaster(); m++)
          {
            // We can simply check if the hanging scheme lists itself as a
            // master
            if (hang_info_pt[i]->master_node_pt(m) ==
                dependent_node_pt[non_vertex_pos[i]])
            {
              node_is_really_shared = true;

#ifdef PARANOID
              // Also check the corresponding weight: it should be 1
              if (std::fabs(hang_info_pt[i]->master_weight(m) - 1.0) > 1.0e-06)
              {
                throw OomphLibError(
                  "Something fishy here -- with shared node at a mortar vertex",
                  "PRefineableQElement<3,INITIAL_NNODE_1D>::quad_hang_helper()",
                  OOMPH_EXCEPTION_LOCATION);
              }
#endif
            }
          }

          // Now we can make the node hang, if it isn't a shared node
          if (!node_is_really_shared)
          {
            dependent_node_pt[non_vertex_pos[i]]->set_hanging_pt(
              hang_info_pt[i], -1);
          }
        }

      } // End of case where there are still dependent nodes

#ifdef PARANOID
      // Check all dependent nodes, hanging or otherwise
      for (unsigned i = 0; i < n_dependent_nodes; i++)
      {
        // Check that weights sum to 1 for those that hang
        if (dependent_node_pt[i]->is_hanging())
        {
          // Check that weights sum to 1 and no zero weights
          double weight_sum = 0.0;
          bool zero_weight = false;
          for (unsigned m = 0;
               m < dependent_node_pt[i]->hanging_pt()->nmaster();
               m++)
          {
            weight_sum += dependent_node_pt[i]->hanging_pt()->master_weight(m);
            if (std::fabs(dependent_node_pt[i]->hanging_pt()->master_weight(
                  m)) < 1.0e-14)
            {
              zero_weight = true;
              oomph_info << "In the hanging scheme for dependent node " << i
                         << ", master node " << m << " has weight "
                         << dependent_node_pt[i]->hanging_pt()->master_weight(m)
                         << " < 1.0e-14" << std::endl;
            }
          }

          // Warn if not
          if (std::fabs(weight_sum - 1.0) > 1.0e-08)
          {
            oomph_info << "Sum of master weights: " << weight_sum << std::endl;
            OomphLibWarning(
              "Weights in hanging scheme do not sum to 1",
              "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
              OOMPH_EXCEPTION_LOCATION);
          }
          if (zero_weight)
          {
            OomphLibWarning(
              "Zero weights present in hanging schemes",
              "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
              OOMPH_EXCEPTION_LOCATION);
          }

          // Also check that dependent nodes do not have themselves as masters
          // bool dependent_should_not_be_hanging = false;
          for (unsigned m = 0;
               m < dependent_node_pt[i]->hanging_pt()->nmaster();
               m++)
          {
            if (dependent_node_pt[i] ==
                dependent_node_pt[i]->hanging_pt()->master_node_pt(m))
            {
              // This shouldn't happen!
              throw OomphLibError(
                "Dependent node has itself as a master!",
                "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
                OOMPH_EXCEPTION_LOCATION);
            }
            if (dependent_node_pt[i]
                  ->hanging_pt()
                  ->master_node_pt(m)
                  ->is_hanging())
            {
              // Check if node is master of master
              Node* new_nod_pt =
                dependent_node_pt[i]->hanging_pt()->master_node_pt(m);
              for (unsigned mm = 0; mm < new_nod_pt->hanging_pt()->nmaster();
                   mm++)
              {
                if (dependent_node_pt[i] ==
                    new_nod_pt->hanging_pt()->master_node_pt(mm))
                {
                  std::cout << "++++++++++++++++++++++++++++++++++++++++"
                            << std::endl;
                  std::cout
                    << "          Dependent node: " << dependent_node_pt[i]
                    << " = " << dependent_node_pt[i]->x(0) << " "
                    << dependent_node_pt[i]->x(1) << " "
                    << dependent_node_pt[i]->x(2) << " " << std::endl;
                  std::cout
                    << "         Master node: "
                    << dependent_node_pt[i]->hanging_pt()->master_node_pt(m)
                    << " = "
                    << dependent_node_pt[i]->hanging_pt()->master_node_pt(m)->x(
                         0)
                    << " "
                    << dependent_node_pt[i]->hanging_pt()->master_node_pt(m)->x(
                         1)
                    << " "
                    << dependent_node_pt[i]->hanging_pt()->master_node_pt(m)->x(
                         2)
                    << " " << std::endl;
                  std::cout << "Master's master node: "
                            << dependent_node_pt[i]
                                 ->hanging_pt()
                                 ->master_node_pt(m)
                                 ->hanging_pt()
                                 ->master_node_pt(mm)
                            << " = "
                            << dependent_node_pt[i]
                                 ->hanging_pt()
                                 ->master_node_pt(m)
                                 ->hanging_pt()
                                 ->master_node_pt(mm)
                                 ->x(0)
                            << " "
                            << dependent_node_pt[i]
                                 ->hanging_pt()
                                 ->master_node_pt(m)
                                 ->hanging_pt()
                                 ->master_node_pt(mm)
                                 ->x(1)
                            << " "
                            << dependent_node_pt[i]
                                 ->hanging_pt()
                                 ->master_node_pt(m)
                                 ->hanging_pt()
                                 ->master_node_pt(mm)
                                 ->x(2)
                            << " " << std::endl;


                  // Print out hanging scheme
                  std::cout << "Hanging node: " << dependent_node_pt[i]->x(0)
                            << "  " << dependent_node_pt[i]->x(1) << "  "
                            << dependent_node_pt[i]->x(2) << "  " << std::endl;
                  for (unsigned m_tmp = 0;
                       m_tmp < dependent_node_pt[i]->hanging_pt()->nmaster();
                       m_tmp++)
                  {
                    std::cout
                      << "   m = " << m_tmp << ": "
                      << dependent_node_pt[i]
                           ->hanging_pt()
                           ->master_node_pt(m_tmp)
                           ->x(0)
                      << "  "
                      << dependent_node_pt[i]
                           ->hanging_pt()
                           ->master_node_pt(m_tmp)
                           ->x(1)
                      << "  "
                      << dependent_node_pt[i]
                           ->hanging_pt()
                           ->master_node_pt(m_tmp)
                           ->x(2)
                      << "  "
                      << "w = "
                      << dependent_node_pt[i]->hanging_pt()->master_weight(
                           m_tmp)
                      << "  " << std::endl;
                  }

                  // Print out hanging scheme
                  std::cout << "Master node " << m
                            << " of Hanging node: " << new_nod_pt->x(0) << "  "
                            << new_nod_pt->x(1) << "  " << new_nod_pt->x(2)
                            << "  " << std::endl;
                  for (unsigned mm_tmp = 0;
                       mm_tmp < new_nod_pt->hanging_pt()->nmaster();
                       mm_tmp++)
                  {
                    std::cout
                      << "   mm = " << mm_tmp << ": "
                      << new_nod_pt->hanging_pt()->master_node_pt(mm_tmp)->x(0)
                      << "  "
                      << new_nod_pt->hanging_pt()->master_node_pt(mm_tmp)->x(1)
                      << "  "
                      << new_nod_pt->hanging_pt()->master_node_pt(mm_tmp)->x(2)
                      << "  "
                      << "w = "
                      << new_nod_pt->hanging_pt()->master_weight(mm_tmp) << "  "
                      << std::endl;
                  }

                  // This shouldn't happen!
                  throw OomphLibError(
                    "Dependent node has itself as a master of a master!",
                    "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
                    OOMPH_EXCEPTION_LOCATION);
                }
              }
            }
          }
        }
        else
        {
          // Check that this node is shared with the master element if it
          // isn't hanging
          bool is_master = false;
          for (unsigned n = 0; n < n_master_nodes; n++)
          {
            if (dependent_node_pt[i] == master_node_pt[n])
            {
              // Node is a master
              is_master = true;
              break;
            }
          }

          if (!is_master)
          {
            /// / Throw error
            // std::ostringstream error_string;
            // error_string
            // << "This node in the dependent element is neither" << std::endl
            // << "hanging or shared with a master element." << std::endl;
            //
            // throw OomphLibError(
            //       error_string.str(),
            //       "PRefineableQElement<3,INITIAL_NNODE_1D>::quad_hang_helper()",
            //       OOMPH_EXCEPTION_LOCATION);
          }
        }
      }
#endif

      // Finally, Loop over all dependent nodes and fine-tune their positions
      //-----------------------------------------------------------------
      // Here we simply set the node's positions to be consistent
      // with the hanging scheme. This is not strictly necessary
      // because it is done in the mesh adaptation before the node
      // becomes non-hanging later on. We make no attempt to ensure
      // (strong) continuity in the position across the mortar.
      for (unsigned i = 0; i < n_dependent_nodes; i++)
      {
        // Only fine-tune hanging nodes
        if (dependent_node_pt[i]->is_hanging())
        {
          // If we are doing the position, then
          if (value_id == -1)
          {
            // Get the local coordinate of this dependent node
            Vector<double> s_local(3);
            this->local_coordinate_of_node(dependent_node_number[i], s_local);

            // Get the position from interpolation in this element via
            // the hanging scheme
            Vector<double> x_in_neighb(3);
            this->interpolated_x(s_local, x_in_neighb);

            // Fine adjust the coordinates (macro map will pick up boundary
            // accurately but will lead to different element edges)
            dependent_node_pt[i]->x(0) = x_in_neighb[0];
            dependent_node_pt[i]->x(1) = x_in_neighb[1];
            dependent_node_pt[i]->x(2) = x_in_neighb[2];
          }
        }
      }
    } // End of case where this interface is to be mortared
  }

  //===================================================================
  // Build required templates
  //===================================================================
  template class PRefineableQElement<1, 2>;
  template class PRefineableQElement<1, 3>;
  template class PRefineableQElement<1, 4>;

  template class PRefineableQElement<2, 2>;
  template class PRefineableQElement<2, 3>;
  template class PRefineableQElement<2, 4>;

  template class PRefineableQElement<3, 2>;
  template class PRefineableQElement<3, 3>;
  template class PRefineableQElement<3, 4>;

} // namespace oomph
