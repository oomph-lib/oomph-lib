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
// Functions for the ElementWithMovingNode class
#include "element_with_moving_nodes.h"
#include "geom_objects.h"
#include "algebraic_elements.h"

namespace oomph
{
  /////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////
  // Functions for the ElementWithMovingNodes class
  ///////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////

  //====================================================================
  /// Return a set of all geometric data associated with the element's node
  /// update function
  //======================================================================
  void ElementWithMovingNodes::assemble_set_of_all_geometric_data(
    std::set<Data*>& unique_geom_data_pt)
  {
    // First clear the set (just in case)
    unique_geom_data_pt.clear();

    // Get number of nodes
    const unsigned n_node = this->nnode();

    // Loop over all nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // Cache pointer to the Node
      Node* const nod_pt = this->node_pt(n);

      // Is the node hanging
      const bool node_is_hanging = nod_pt->is_hanging();

      // Default number of master nodes
      unsigned nmaster = 1;

      // Default: Node isn't hanging so it's its own master node
      Node* master_node_pt = nod_pt;

      // Cache the hanging point
      HangInfo* hang_info_pt = 0;

      // Find the number of master nodes if the node is hanging
      if (node_is_hanging)
      {
        hang_info_pt = nod_pt->hanging_pt();
        nmaster = hang_info_pt->nmaster();
      }

      // Loop over all master nodes
      for (unsigned imaster = 0; imaster < nmaster; imaster++)
      {
        // Get the master node
        if (node_is_hanging)
        {
          master_node_pt = hang_info_pt->master_node_pt(imaster);
        }


        // Find the number of data
        const unsigned n_geom_data = master_node_pt->ngeom_data();
        // If there are geometric data add them to the set
        if (n_geom_data > 0)
        {
          // Get vector of geometric data involved in the geometric
          // change of this node
          Data** node_geom_data_pt = master_node_pt->all_geom_data_pt();

          for (unsigned i = 0; i < n_geom_data; i++)
          {
            unique_geom_data_pt.insert(node_geom_data_pt[i]);
          }
        }

        // Find the number of geometric objects
        unsigned n_geom_obj = master_node_pt->ngeom_object();

        // If there are geometric objects, add them to the set
        if (n_geom_obj > 0)
        {
          // Get vector of geometric objects involved in the default
          // update function for this (master) node.
          // Vector is constructed by copy operation.
          GeomObject** geom_object_pt = master_node_pt->all_geom_object_pt();

          // Loop over the geometric objects
          for (unsigned i = 0; i < n_geom_obj; i++)
          {
            // Get the next geometric object
            GeomObject* geom_obj_pt = geom_object_pt[i];

            // Number of items of geometric data that the geometric
            // object depends on
            unsigned n_geom_data = geom_obj_pt->ngeom_data();

            // Loop over geometric data and add to set (use set to ensure
            // that each one is only counted once)
            for (unsigned idata = 0; idata < n_geom_data; idata++)
            {
              unique_geom_data_pt.insert(geom_obj_pt->geom_data_pt(idata));
            }
          }
        } // End of add geom object loop
      }
    }
  }

  //=================================================================
  /// Construct the vector of (unique) geometric data
  //=================================================================
  void ElementWithMovingNodes::complete_setup_of_dependencies()
  {
    // This set will hold the pointers to all the unique (geometric) Data that
    // affects the shape of this element
    std::set<Data*> unique_geom_data_pt;
    // Assemble that data
    this->assemble_set_of_all_geometric_data(unique_geom_data_pt);

    // Resize storage for the pointers to the Data items that are
    // involved in the element's node update operation.
    Geom_data_pt.clear();

    // Loop over all the unique remaining Data items involved in the
    // node update operations
    typedef std::set<Data*>::iterator IT;
    for (IT it = unique_geom_data_pt.begin(); it != unique_geom_data_pt.end();
         it++)
    {
      Geom_data_pt.push_back(*it);
    }
  }

  //==================================================================
  /// Function to describe the local dofs of the element[s]. The ostream
  /// specifies the output stream to which the description
  /// is written; the string stores the currently
  /// assembled output that is ultimately written to the
  /// output stream by Data::describe_dofs(...); it is typically
  /// built up incrementally as we descend through the
  /// call hierarchy of this function when called from
  /// Problem::describe_dofs(...)
  //==================================================================
  void ElementWithMovingNodes::describe_local_dofs(
    std::ostream& out, const std::string& current_string) const
  {
    // Call the standard finite element classification.
    FiniteElement::describe_local_dofs(out, current_string);

    // Set the number of data
    const unsigned n_geom_data = ngeom_data();

    // Loop over the node update data
    for (unsigned i = 0; i < n_geom_data; i++)
    {
      // Pointer to geometric Data
      Data* data_pt = Geom_data_pt[i];

      std::stringstream conversion;
      conversion << " of Geometric Data " << i << current_string;
      std::string in(conversion.str());
      data_pt->describe_dofs(out, in);
    }
  }


  //==================================================================
  /// Assign local equation numbers for the geometric data associated
  /// with the element.
  //==================================================================
  void ElementWithMovingNodes::assign_all_generic_local_eqn_numbers(
    const bool& store_local_dof_pt)
  {
    // Get local number of dofs so far
    unsigned local_eqn_number = this->ndof();

    // Set the number of data
    const unsigned n_geom_data = ngeom_data();

    // Reset number of geometric dofs
    Ngeom_dof = 0;

    // If we have any geometric data
    if (n_geom_data > 0)
    {
      // Work out total number of values involved
      // Initialise from the first object
      unsigned n_total_values = Geom_data_pt[0]->nvalue();

      // Add the values from the other data
      for (unsigned i = 1; i < n_geom_data; i++)
      {
        n_total_values += Geom_data_pt[i]->nvalue();
      }

      // If allocated delete the old storage
      if (Geometric_data_local_eqn)
      {
        delete[] Geometric_data_local_eqn[0];
        delete[] Geometric_data_local_eqn;
      }

      // If there are no values, we are done, null out the storage and
      // return
      if (n_total_values == 0)
      {
        Geometric_data_local_eqn = 0;
        return;
      }

      // Resize the storage for the geometric data local equation numbers
      // Firstly allocate the row pointers
      Geometric_data_local_eqn = new int*[n_geom_data];

      // Now allocate storage for all the equation numbers
      Geometric_data_local_eqn[0] = new int[n_total_values];

      // Initially all local equations are unclassified
      for (unsigned i = 0; i < n_total_values; i++)
      {
        Geometric_data_local_eqn[0][i] = Data::Is_unclassified;
      }

      // Loop over the remaining rows and set their pointers
      for (unsigned i = 1; i < n_geom_data; ++i)
      {
        // Initially set the pointer to the i-th row to the pointer
        // to the i-1th row
        Geometric_data_local_eqn[i] = Geometric_data_local_eqn[i - 1];

        // Now increase the row pointer by the number of values
        // stored at the i-1th geometric data
        Geometric_data_local_eqn[i] += Geom_data_pt[i - 1]->nvalue();
      }

      // A local queue to store the global equation numbers
      std::deque<unsigned long> global_eqn_number_queue;

      // Loop over the node update data
      for (unsigned i = 0; i < n_geom_data; i++)
      {
        // Pointer to geometric Data
        Data* data_pt = Geom_data_pt[i];

        // Loop over values at this Data item
        unsigned n_value = data_pt->nvalue();
        for (unsigned j = 0; j < n_value; j++)
        {
          // Get global equation number
          long eqn_number = data_pt->eqn_number(j);

          // If equation number positive
          if (eqn_number >= 0)
          {
            // Add the global equation number to our queue
            global_eqn_number_queue.push_back(eqn_number);
            // Add pointer to the dof to the queue if required
            if (store_local_dof_pt)
            {
              GeneralisedElement::Dof_pt_deque.push_back(data_pt->value_pt(j));
            }

            // Add to local value
            Geometric_data_local_eqn[i][j] = local_eqn_number;
            local_eqn_number++;

            // Bump up number of geometric dofs
            Ngeom_dof++;
          }
          else
          {
            // Set the local scheme to be pinned
            Geometric_data_local_eqn[i][j] = Data::Is_pinned;
          }
        }
      }

      // Now add our global equations numbers to the internal element storage
      this->add_global_eqn_numbers(global_eqn_number_queue,
                                   GeneralisedElement::Dof_pt_deque);
      // Clear the memory used in the deque
      if (store_local_dof_pt)
      {
        std::deque<double*>().swap(GeneralisedElement::Dof_pt_deque);
      }
    }
  }


  //==================================================================
  /// Calculate the node-update--related entries in the
  /// Jacobian. The vector passed
  /// in residuals has to contain the nonlinear residuals,
  /// evaluated for the current values of the unknowns, in
  /// case FDing is used to computed the derivatives.
  //==================================================================
  void ElementWithMovingNodes::fill_in_jacobian_from_geometric_data(
    Vector<double>& residuals, DenseMatrix<double>& jacobian)
  {
    if (!Bypass_fill_in_jacobian_from_geometric_data)
    {
      // Get number of Data items involved in node update operations
      const unsigned n_geometric_data = ngeom_data();

      // If there is nothing to be done, then leave
      if (n_geometric_data == 0) return;

      // Number of dofs
      const unsigned n_dof = this->ndof();

      // Number of nodes
      unsigned n_nod = this->nnode();

      // If there are no dofs, return
      if (n_nod == 0) return;

      // Get nodal dimension from first node
      const unsigned dim_nod = node_pt(0)->ndim();

      // Number of shape controlling nodes for nonrefineable elements
      unsigned n_shape_controlling_node = nnode();

      // Are we dealing with a refineable element?
      RefineableElement* ref_el_pt = dynamic_cast<RefineableElement*>(this);
      if (ref_el_pt != 0)
      {
        // Adjust number of shape controlling nodes
        n_shape_controlling_node = ref_el_pt->nshape_controlling_nodes();
      }

      // How are we going to evaluate the shape derivs?
      unsigned method = 0;
      if (Method_for_shape_derivs == Shape_derivs_by_direct_fd)
      {
        method = 0;
      }
      else if (Method_for_shape_derivs == Shape_derivs_by_chain_rule)
      {
        method = 1;
      }
      else if (Method_for_shape_derivs == Shape_derivs_by_fastest_method)
      {
        // Direct FD-ing of residuals w.r.t. geometric dofs is likely to be
        // faster if there are fewer geometric dofs than total nodal coordinates
        // (nodes x dim) in element:
        if (Ngeom_dof < (n_shape_controlling_node * dim_nod))
        {
          method = 0;
        }
        else
        {
          method = 1;
        }
      }

      // Choose method
      //===============
      switch (method)
      {
          // Direct FD:
          //-----------
        case 0:

        {
          // Create newres vector
          Vector<double> newres(n_dof);

          // Use the default finite difference step
          const double fd_step = GeneralisedElement::Default_fd_jacobian_step;

          // Integer storage for the local unknown
          int local_unknown = 0;

          // Loop over the Data items that affect the node update operations
          for (unsigned i = 0; i < n_geometric_data; i++)
          {
            // Loop over values
            unsigned n_value = Geom_data_pt[i]->nvalue();
            for (unsigned j = 0; j < n_value; j++)
            {
              local_unknown = geometric_data_local_eqn(i, j);

              // If the value is free
              if (local_unknown >= 0)
              {
                // Get a pointer to the geometric data value
                double* value_pt = Geom_data_pt[i]->value_pt(j);

                // Save the old value
                double old_var = *value_pt;

                // Increment the variable
                *value_pt += fd_step;

                // Update the whole element (Bit inefficient)
                this->node_update();

                // Calculate the new residuals
                this->get_residuals(newres);

                // Now do finite differences
                for (unsigned m = 0; m < n_dof; m++)
                {
                  // Stick the entry into the Jacobian matrix
                  jacobian(m, local_unknown) =
                    (newres[m] - residuals[m]) / fd_step;
                }

                // Reset the variable
                *value_pt = old_var;

                // We're relying on the total node update in the next loop
              }
            }
          }

          // Node update the element one final time to get things back to
          // the original state
          this->node_update();
        }

        break;

        // Chain rule
        //-----------
        case 1:

        {
          // Get derivatives of residuals w.r.t. all nodal coordinates
          RankThreeTensor<double> dresidual_dnodal_coordinates(
            n_dof, dim_nod, n_shape_controlling_node, 0.0);

          // Use FD-version in base class?
          if (Evaluate_dresidual_dnodal_coordinates_by_fd)
          {
            if (ref_el_pt != 0)
            {
              ref_el_pt->RefineableElement::get_dresidual_dnodal_coordinates(
                dresidual_dnodal_coordinates);
            }
            else
            {
              FiniteElement::get_dresidual_dnodal_coordinates(
                dresidual_dnodal_coordinates);
            }
          }
          // Otherwise use the overloaded analytical version in derived
          // class (if it exists -- if it doesn't this just drops through
          // to the default implementation in FiniteElement).
          else
          {
            this->get_dresidual_dnodal_coordinates(
              dresidual_dnodal_coordinates);
          }

          // Get derivatives of nodal coordinates w.r.t. geometric dofs
          RankThreeTensor<double> dnodal_coordinates_dgeom_dofs(
            n_dof, dim_nod, n_shape_controlling_node, 0.0);

          get_dnodal_coordinates_dgeom_dofs(dnodal_coordinates_dgeom_dofs);

          // Assemble Jacobian via chain rule
          for (unsigned l = 0; l < n_dof; l++)
          {
            // Loop over the Data items that affect the node update operations
            for (unsigned i_data = 0; i_data < n_geometric_data; i_data++)
            {
              // Loop over values
              unsigned n_value = Geom_data_pt[i_data]->nvalue();
              for (unsigned j_val = 0; j_val < n_value; j_val++)
              {
                int k = geometric_data_local_eqn(i_data, j_val);

                // If the value is free
                if (k >= 0)
                {
                  jacobian(l, k) = 0.0;
                  for (unsigned i = 0; i < dim_nod; i++)
                  {
                    for (unsigned j = 0; j < n_shape_controlling_node; j++)
                    {
                      jacobian(l, k) += dresidual_dnodal_coordinates(l, i, j) *
                                        dnodal_coordinates_dgeom_dofs(k, i, j);
                    }
                  }
                }
              }
            }
          }
        }

        break;

        default:

          std::ostringstream error_message;
          error_message << "Never get here: method " << method;
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }
  }


  //======================================================================
  /// Compute derivatives of the nodal coordinates w.r.t.
  /// to the geometric dofs. Default implementation by FD can be overwritten
  /// for specific elements.
  /// dnodal_coordinates_dgeom_dofs(l,i,j) = dX_{ij} / d s_l
  //======================================================================
  void ElementWithMovingNodes::get_dnodal_coordinates_dgeom_dofs(
    RankThreeTensor<double>& dnodal_coordinates_dgeom_dofs)
  {
    // Get number of Data items involved in node update operations
    const unsigned n_geometric_data = ngeom_data();

    // If there is nothing to be done, then leave
    if (n_geometric_data == 0)
    {
      return;
    }

    // Number of nodes
    const unsigned n_nod = nnode();

    // If the element has no nodes (why??!!) return straightaway
    if (n_nod == 0) return;

    // Get dimension from first node
    unsigned dim_nod = node_pt(0)->ndim();

    // Number of shape controlling nodes for nonrefineable elements
    unsigned n_shape_controlling_node = n_nod;

    // Are we dealing with a refineable element?
    RefineableElement* ref_el_pt = dynamic_cast<RefineableElement*>(this);
    if (ref_el_pt != 0)
    {
      // Adjust number of shape controlling nodes
      n_shape_controlling_node = ref_el_pt->nshape_controlling_nodes();
    }

    // Current and advanced nodal positions
    DenseMatrix<double> pos(dim_nod, n_shape_controlling_node);

    // Shape controlling nodes
    std::map<Node*, unsigned> local_shape_controlling_node_lookup;

    // Refineable element:
    if (ref_el_pt != 0)
    {
      local_shape_controlling_node_lookup =
        ref_el_pt->shape_controlling_node_lookup();
    }
    // Non-refineable element: the nodes themselves
    else
    {
      unsigned count = 0;
      for (unsigned j = 0; j < n_nod; j++)
      {
        local_shape_controlling_node_lookup[node_pt(j)] = count;
        count++;
      }
    }

    // Loop over all shape-controlling nodes to backup their original position
    for (std::map<Node*, unsigned>::iterator it =
           local_shape_controlling_node_lookup.begin();
         it != local_shape_controlling_node_lookup.end();
         it++)
    {
      // Get node
      Node* nod_pt = it->first;

      // Get its number
      unsigned node_number = it->second;

      // Backup
      for (unsigned i = 0; i < dim_nod; i++)
      {
        pos(i, node_number) = nod_pt->position(i);
      }
    }


    // Integer storage for the local unknown
    int local_unknown = 0;

    // Use the default finite difference step
    const double fd_step = GeneralisedElement::Default_fd_jacobian_step;

    // Loop over the Data items that affect the node update operations
    for (unsigned i = 0; i < n_geometric_data; i++)
    {
      // Loop over values
      unsigned n_value = Geom_data_pt[i]->nvalue();
      for (unsigned j = 0; j < n_value; j++)
      {
        local_unknown = geometric_data_local_eqn(i, j);

        // If the value is free
        if (local_unknown >= 0)
        {
          // Get a pointer to the geometric data value
          double* value_pt = Geom_data_pt[i]->value_pt(j);

          // Save the old value
          double old_var = *value_pt;

          // Increment the variable
          *value_pt += fd_step;

          // Update the whole element
          this->node_update();

          // Loop over all shape-controlling nodes
          for (std::map<Node*, unsigned>::iterator it =
                 local_shape_controlling_node_lookup.begin();
               it != local_shape_controlling_node_lookup.end();
               it++)
          {
            // Get node
            Node* nod_pt = it->first;

            // Get its number
            unsigned node_number = it->second;

            // Get advanced position and FD
            for (unsigned ii = 0; ii < dim_nod; ii++)
            {
              dnodal_coordinates_dgeom_dofs(local_unknown, ii, node_number) =
                (nod_pt->position(ii) - pos(ii, node_number)) / fd_step;
            }
          }

          // Reset the variable
          *value_pt = old_var;

          // We're relying on the total node update in the next loop
        }
      }
    }
    // Node update the element one final time to get things back to
    // the original state
    this->node_update();
  }

} // namespace oomph
