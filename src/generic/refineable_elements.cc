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
// Non-inline member functions for refineable elements

// oomph-lib includes
#include "refineable_elements.h"
#include "shape.h"

namespace oomph
{
  //=====================================================================
  /// Destructor (needed here because of IBM xlC compiler under AIX)
  /// Delete the storage allocated for the local equations.
  //====================================================================
  RefineableElement::~RefineableElement()
  {
    delete[] Local_hang_eqn;
  }

  //========================================================================
  /// Max. allowed discrepancy in element integrity check
  //========================================================================
  double RefineableElement::Max_integrity_tolerance = 1.0e-8;

  //=========================================================================
  /// Helper function that is used to check that the value_id is in the range
  /// allowed by the element. The number of continuously interpolated values
  /// and the value_id must be passed as arguments.
  //========================================================================
  void RefineableElement::check_value_id(
    const int& n_continuously_interpolated_values, const int& value_id)
  {
    // If the value_id is more than the number of continuously interpolated
    // values or less than -1 (for the position), we have a problem
    if ((value_id < -1) || (value_id >= n_continuously_interpolated_values))
    {
      std::ostringstream error_stream;
      error_stream << "Value_id " << value_id << " is out of range."
                   << std::endl
                   << "It can only take the values -1 (position) "
                   << "or an integer in the range 0 to "
                   << n_continuously_interpolated_values << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  //=========================================================================
  /// Internal function that is used to assemble the jacobian of the mapping
  /// from local coordinates (s) to the eulerian coordinates (x), given the
  /// derivatives of the shape functions.
  //=========================================================================
  void RefineableElement::assemble_local_to_eulerian_jacobian(
    const DShape& dpsids, DenseMatrix<double>& jacobian) const
  {
    // Locally cache the elemental dimension
    const unsigned el_dim = dim();
    // The number of shape functions must be equal to the number
    // of nodes (by definition)
    const unsigned n_shape = nnode();
    // The number of shape function types must be equal to the number
    // of nodal position types (by definition)
    const unsigned n_shape_type = nnodal_position_type();

#ifdef PARANOID
    // Check for dimensional compatibility
    if (dim() != nodal_dimension())
    {
      std::ostringstream error_message;
      error_message << "Dimension mismatch" << std::endl;
      error_message << "The elemental dimension: " << dim()
                    << " must equal the nodal dimension: " << nodal_dimension()
                    << " for the jacobian of the mapping to be well-defined"
                    << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Loop over the rows of the jacobian
    for (unsigned i = 0; i < el_dim; i++)
    {
      // Loop over the columns of the jacobian
      for (unsigned j = 0; j < el_dim; j++)
      {
        // Zero the entry
        jacobian(i, j) = 0.0;
        // Loop over the shape functions
        for (unsigned l = 0; l < n_shape; l++)
        {
          for (unsigned k = 0; k < n_shape_type; k++)
          {
            // Jacobian is dx_j/ds_i, which is represented by the sum
            // over the dpsi/ds_i of the nodal points X j
            // Call the Hanging version of positions
            jacobian(i, j) += nodal_position_gen(l, k, j) * dpsids(l, k, i);
          }
        }
      }
    }
  }

  //=========================================================================
  /// Internal function that is used to assemble the jacobian of second
  /// derivatives of the the mapping from local coordinates (s) to the
  /// eulerian coordinates (x), given the second derivatives of the
  /// shape functions.
  //=========================================================================
  void RefineableElement::assemble_local_to_eulerian_jacobian2(
    const DShape& d2psids, DenseMatrix<double>& jacobian2) const
  {
    // Find the the dimension of the element
    const unsigned el_dim = dim();
    // Find the number of shape functions and shape functions types
    // Must be equal to the number of nodes and their position types
    // by the definition of the shape function.
    const unsigned n_shape = nnode();
    const unsigned n_shape_type = nnodal_position_type();
    // Find the number of second derivatives
    const unsigned n_row = N2deriv[el_dim];

    // Assemble the "jacobian" (d^2 x_j/ds_i^2) for second derivatives of
    // shape functions
    // Loop over the rows (number of second derivatives)
    for (unsigned i = 0; i < n_row; i++)
    {
      // Loop over the columns (element dimension
      for (unsigned j = 0; j < el_dim; j++)
      {
        // Zero the entry
        jacobian2(i, j) = 0.0;
        // Loop over the shape functions
        for (unsigned l = 0; l < n_shape; l++)
        {
          // Loop over the shape function types
          for (unsigned k = 0; k < n_shape_type; k++)
          {
            // Add the terms to the jacobian entry
            // Call the Hanging version of positions
            jacobian2(i, j) += nodal_position_gen(l, k, j) * d2psids(l, k, i);
          }
        }
      }
    }
  }

  //=====================================================================
  /// Assemble the covariant Eulerian base vectors and return them in the
  /// matrix interpolated_G. The derivatives of the shape functions with
  /// respect to the local coordinate should already have been calculated
  /// before calling this function
  //=====================================================================
  void RefineableElement::assemble_eulerian_base_vectors(
    const DShape& dpsids, DenseMatrix<double>& interpolated_G) const
  {
    // Find the number of nodes and position types
    const unsigned n_node = nnode();
    const unsigned n_position_type = nnodal_position_type();
    // Find the dimension of the node and element
    const unsigned n_dim_node = nodal_dimension();
    const unsigned n_dim_element = dim();

    // Loop over the dimensions and compute the entries of the
    // base vector matrix
    for (unsigned i = 0; i < n_dim_element; i++)
    {
      for (unsigned j = 0; j < n_dim_node; j++)
      {
        // Initialise the j-th component of the i-th base vector to zero
        interpolated_G(i, j) = 0.0;
        for (unsigned l = 0; l < n_node; l++)
        {
          for (unsigned k = 0; k < n_position_type; k++)
          {
            interpolated_G(i, j) +=
              nodal_position_gen(l, k, j) * dpsids(l, k, i);
          }
        }
      }
    }
  }


  //==========================================================================
  /// Calculate the mapping from local to eulerian coordinates
  /// assuming that the coordinates are aligned in the direction of the local
  /// coordinates, i.e. there are no cross terms and the jacobian is diagonal.
  /// The local derivatives are passed as dpsids and the jacobian and
  /// inverse jacobian are returned.
  //==========================================================================
  double RefineableElement::local_to_eulerian_mapping_diagonal(
    const DShape& dpsids,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& inverse_jacobian) const
  {
    // Find the dimension of the element
    const unsigned el_dim = dim();
    // Find the number of shape functions and shape functions types
    // Equal to the number of nodes and their position types by definition
    const unsigned n_shape = nnode();
    const unsigned n_shape_type = nnodal_position_type();

#ifdef PARANOID
    // Check for dimension compatibility
    if (dim() != nodal_dimension())
    {
      std::ostringstream error_message;
      error_message << "Dimension mismatch" << std::endl;
      error_message << "The elemental dimension: " << dim()
                    << " must equal the nodal dimension: " << nodal_dimension()
                    << " for the jacobian of the mapping to be well-defined"
                    << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // In this case we assume that there are no cross terms, that is
    // global coordinate 0 is always in the direction of local coordinate 0

    // Loop over the coordinates
    for (unsigned i = 0; i < el_dim; i++)
    {
      // Zero the jacobian and inverse jacobian entries
      for (unsigned j = 0; j < el_dim; j++)
      {
        jacobian(i, j) = 0.0;
        inverse_jacobian(i, j) = 0.0;
      }

      // Loop over the shape functions
      for (unsigned l = 0; l < n_shape; l++)
      {
        // Loop over the types of dof
        for (unsigned k = 0; k < n_shape_type; k++)
        {
          // Derivatives are always dx_{i}/ds_{i}
          jacobian(i, i) += nodal_position_gen(l, k, i) * dpsids(l, k, i);
        }
      }
    }

    // Now calculate the determinant of the matrix
    double det = 1.0;
    for (unsigned i = 0; i < el_dim; i++)
    {
      det *= jacobian(i, i);
    }

// Report if Matrix is singular, or negative
#ifdef PARANOID
    check_jacobian(det);
#endif

    // Calculate the inverse mapping (trivial in this case)
    for (unsigned i = 0; i < el_dim; i++)
    {
      inverse_jacobian(i, i) = 1.0 / jacobian(i, i);
    }

    // Return the value of the Jacobian
    return (det);
  }

  //========================================================================
  /// Deactivate the element by marking setting all local pointers to
  /// obsolete nodes to zero
  //=======================================================================
  void RefineableElement::deactivate_element()
  {
    // Find the number of nodes
    const unsigned n_node = nnode();
    // Loop over the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // If the node pointer has not already been nulled, but is obsolete
      // then null it
      if ((this->node_pt(n) != 0) && (this->node_pt(n)->is_obsolete()))
      {
        this->node_pt(n) = 0;
      }
    }
  }

  //=======================================================================
  /// Assign the local hang eqn.
  //=======================================================================
  void RefineableElement::assign_hanging_local_eqn_numbers(
    const bool& store_local_dof_pt)
  {
    // Find the number of nodes
    const unsigned n_node = nnode();

    // Check there are nodes!
    if (n_node > 0)
    {
      // Find the number of continuously interpolated values
      const unsigned n_cont_values = ncont_interpolated_values();

      // Delete the storage associated with the local hanging equation schemes
      delete[] Local_hang_eqn;
      // Allocate new storage
      Local_hang_eqn = new std::map<Node*, int>[n_cont_values];

      // Boolean that is set to true if there are hanging equation numbers
      bool hanging_eqn_numbers = false;

      // Maps that store whether the node's local equation number for a
      // particular value has already been assigned
      Vector<std::map<Node*, bool>> local_eqn_number_done(n_cont_values);

      // Get number of dofs assigned thus far
      unsigned local_eqn_number = ndof();

      // A local queue to store the global equation numbers
      std::deque<unsigned long> global_eqn_number_queue;

      // Now loop over all the nodes again to find the master nodes
      // external to the element
      for (unsigned n = 0; n < n_node; n++)
      {
        // Loop over the number of continuously-interpolated values at the node
        for (unsigned j = 0; j < n_cont_values; j++)
        {
          // If the node is hanging node in value j
          if (node_pt(n)->is_hanging(j))
          {
            // Get the pointer to the appropriate hanging info object
            HangInfo* hang_info_pt = node_pt(n)->hanging_pt(j);
            // Find the number of master nodes
            unsigned n_master = hang_info_pt->nmaster();
            // Loop over the master nodes
            for (unsigned m = 0; m < n_master; m++)
            {
              // Get the m-th master node
              Node* Master_node_pt = hang_info_pt->master_node_pt(m);

              // If the master node's value has not been considered already,
              // give it a local equation number
              if (local_eqn_number_done[j][Master_node_pt] == false)
              {
#ifdef PARANOID
                // Check that the value is stored at the master node
                // If not the master node must have be set up incorrectly
                if (Master_node_pt->nvalue() < j)
                {
                  std::ostringstream error_stream;
                  error_stream << "Master node for " << j
                               << "-th value only has "
                               << Master_node_pt->nvalue() << " stored values!"
                               << std::endl;

                  throw OomphLibError(error_stream.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
#endif

                // We now need to test whether the master node is actually
                // a local node, in which case its local equation number
                // will already have been assigned and will be stored
                // in nodal_local_eqn

                // Storage for the index of the local node
                // Initialised to n_node (beyond the possible indices of
                //"real" nodes)
                unsigned local_node_index = n_node;
                // Loop over the local nodes (again)
                for (unsigned n1 = 0; n1 < n_node; n1++)
                {
                  // If the master node is a local node
                  // get its index and break out of the loop
                  if (Master_node_pt == node_pt(n1))
                  {
                    local_node_index = n1;
                    break;
                  }
                }

                // Now we test whether the node was found
                if (local_node_index < n_node)
                {
                  // Copy the local equation number to the
                  // pointer-based look-up scheme
                  Local_hang_eqn[j][Master_node_pt] =
                    nodal_local_eqn(local_node_index, j);
                }
                // Otherwise it's a new master node
                else
                {
                  // Get the equation number
                  long eqn_number = Master_node_pt->eqn_number(j);
                  // If equation_number positive add to array
                  if (eqn_number >= 0)
                  {
                    // Add global equation number to the queue
                    global_eqn_number_queue.push_back(eqn_number);
                    // Add pointer to the dof to the queue if required
                    if (store_local_dof_pt)
                    {
                      GeneralisedElement::Dof_pt_deque.push_back(
                        Master_node_pt->value_pt(j));
                    }
                    // Add to pointer based scheme
                    Local_hang_eqn[j][Master_node_pt] = local_eqn_number;
                    // Increase number of local variables
                    local_eqn_number++;
                  }
                  // Otherwise the value is pinned
                  else
                  {
                    Local_hang_eqn[j][Master_node_pt] = Data::Is_pinned;
                  }
                  // There are now hanging equation numbers
                }
                // The value at this node has now been done
                local_eqn_number_done[j][Master_node_pt] = true;
                // There are hanging equation numbers
                hanging_eqn_numbers = true;
              }
            }
          }
        }
      } // End of second loop over nodes

      // Now add our global equations numbers to the internal element storage
      add_global_eqn_numbers(global_eqn_number_queue,
                             GeneralisedElement::Dof_pt_deque);
      // Clear the memory used in the deque
      if (store_local_dof_pt)
      {
        std::deque<double*>().swap(GeneralisedElement::Dof_pt_deque);
      }


      // If there are no hanging_eqn_numbers delete the (empty) stored maps
      if (!hanging_eqn_numbers)
      {
        delete[] Local_hang_eqn;
        Local_hang_eqn = 0;
      }


      // Setup map that associates a unique number with any of the nodes
      // that actively control the shape of the element (i.e. they are
      // either non-hanging nodes of this element or master nodes
      // of hanging nodes.
      unsigned count = 0;
      std::set<Node*> all_nodes;
      Shape_controlling_node_lookup.clear();
      for (unsigned j = 0; j < n_node; j++)
      {
        // Get node
        Node* nod_pt = node_pt(j);

        // If the node is hanging, consider master nodes
        if (nod_pt->is_hanging())
        {
          HangInfo* hang_info_pt = node_pt(j)->hanging_pt();
          unsigned n_master = hang_info_pt->nmaster();

          // Loop over the master nodes
          for (unsigned m = 0; m < n_master; m++)
          {
            Node* master_node_pt = hang_info_pt->master_node_pt(m);
            // Do we have this one already?
            unsigned old_size = all_nodes.size();
            all_nodes.insert(master_node_pt);
            if (all_nodes.size() > old_size)
            {
              Shape_controlling_node_lookup[master_node_pt] = count;
              count++;
            }
          }
        }
        // Not hanging: Consider the node itself
        else
        {
          // Do we have this one already?
          unsigned old_size = all_nodes.size();
          all_nodes.insert(nod_pt);
          if (all_nodes.size() > old_size)
          {
            Shape_controlling_node_lookup[nod_pt] = count;
            count++;
          }
        }
      }

    } // End of if nodes
  }

  //======================================================================
  /// The purpose of this function is to identify all possible
  /// Data that can affect the fields interpolated by the FiniteElement.
  /// This must be overloaded to include data from any hanging nodes
  /// correctly
  //=======================================================================
  void RefineableElement::identify_field_data_for_interactions(
    std::set<std::pair<Data*, unsigned>>& paired_field_data)
  {
    // Loop over all internal data
    const unsigned n_internal = this->ninternal_data();
    for (unsigned n = 0; n < n_internal; n++)
    {
      // Cache the data pointer
      Data* const dat_pt = this->internal_data_pt(n);
      // Find the number of data values stored in the data object
      const unsigned n_value = dat_pt->nvalue();
      // Add the index of each data value and the pointer to the set
      // of pairs
      for (unsigned i = 0; i < n_value; i++)
      {
        paired_field_data.insert(std::make_pair(dat_pt, i));
      }
    }

    // Loop over all the nodes
    const unsigned n_node = this->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      // Cache the node pointer
      Node* const nod_pt = this->node_pt(n);
      // Find the number of values stored at the node
      const unsigned n_value = nod_pt->nvalue();

      // Loop over the number of values
      for (unsigned i = 0; i < n_value; i++)
      {
        // If the node is hanging in the variable
        if (nod_pt->is_hanging(i))
        {
          // Cache the pointer to the HangInfo object for this variable
          HangInfo* const hang_pt = nod_pt->hanging_pt(i);
          // Get number of master nodes
          const unsigned nmaster = hang_pt->nmaster();

          // Loop over masters
          for (unsigned m = 0; m < nmaster; m++)
          {
            // Cache the pointer to the master node
            Node* const master_nod_pt = hang_pt->master_node_pt(m);

            // Under the assumption that the same data value is interpolated
            // by the hanging nodes, which it really must be
            // Add it to the paired field data
            paired_field_data.insert(std::make_pair(master_nod_pt, i));
          }
        }
        // Otherwise the node is not hanging in the variable, treat as normal
        else
        {
          // Add the index of each data value and the pointer to the set
          // of pairs
          paired_field_data.insert(std::make_pair(nod_pt, i));
        }
      }
    }
  }


  //==========================================================================
  /// Compute derivatives of elemental residual vector with respect
  /// to nodal coordinates. Default implementation by FD can be overwritten
  /// for specific elements.
  /// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
  /// This version is overloaded from the version in FiniteElement
  /// and takes hanging nodes into account -- j in the above loop
  /// loops over all the nodes that actively control the
  /// shape of the element (i.e. they are non-hanging or master nodes of
  /// hanging nodes in this element).
  //==========================================================================
  void RefineableElement::get_dresidual_dnodal_coordinates(
    RankThreeTensor<double>& dresidual_dnodal_coordinates)
  {
    // Number of nodes
    unsigned n_nod = nnode();

    // If the element has no nodes (why??!!) return straightaway
    if (n_nod == 0) return;

    // Get dimension from first node
    unsigned dim_nod = node_pt(0)->ndim();

    // Number of dofs
    unsigned n_dof = ndof();

    // Get reference residual
    Vector<double> res(n_dof);
    Vector<double> res_pls(n_dof);
    get_residuals(res);

    // FD step
    double eps_fd = GeneralisedElement::Default_fd_jacobian_step;

    // Do FD loop over all active nodes
    for (std::map<Node*, unsigned>::iterator it =
           Shape_controlling_node_lookup.begin();
         it != Shape_controlling_node_lookup.end();
         it++)
    {
      // Get node
      Node* nod_pt = it->first;

      // Get its number
      unsigned node_number = it->second;

      // Loop over coordinate directions
      for (unsigned i = 0; i < dim_nod; i++)
      {
        // Make backup
        double backup = nod_pt->x(i);

        // Do FD step. No node update required as we're
        // attacking the coordinate directly...
        nod_pt->x(i) += eps_fd;

        // Perform auxiliary node update function
        nod_pt->perform_auxiliary_node_update_fct();

        // Get advanced residual
        get_residuals(res_pls);

        // Fill in FD entries [Loop order is "wrong" here as l is the
        // slow index but this is in a function that's costly anyway
        // and gives us the fastest loop outside where these tensor
        // is actually used.]
        for (unsigned l = 0; l < n_dof; l++)
        {
          dresidual_dnodal_coordinates(l, i, node_number) =
            (res_pls[l] - res[l]) / eps_fd;
        }

        // Reset coordinate. No node update required as we're
        // attacking the coordinate directly...
        nod_pt->x(i) = backup;

        // Perform auxiliary node update function
        nod_pt->perform_auxiliary_node_update_fct();
      }
    }
  }


  //============================================================================
  /// This function calculates the entries of Jacobian matrix, used in
  /// the Newton method, associated with the nodal degrees of freedom.
  /// This is done by finite differences to handle the general case
  /// Overload the standard case to include hanging node case
  //==========================================================================
  void RefineableElement::fill_in_jacobian_from_nodal_by_fd(
    Vector<double>& residuals, DenseMatrix<double>& jacobian)
  {
    // Find the number of nodes
    const unsigned n_node = nnode();

    // If there are no nodes, return straight away
    if (n_node == 0)
    {
      return;
    }

    // Call the update function to ensure that the element is in
    // a consistent state before finite differencing starts
    update_before_nodal_fd();

    //  bool use_first_order_fd=false;

    // Find the number of dofs in the element
    const unsigned n_dof = ndof();

    // Create newres vector
    Vector<double> newres(n_dof); // , newres_minus(n_dof);

    // Used default value defined in GeneralisedElement
    const double fd_step = Default_fd_jacobian_step;

    // Integer storage for local unknowns
    int local_unknown = 0;

    // Loop over the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // Get the number of values stored at the node
      const unsigned n_value = node_pt(n)->nvalue();

      // Get a pointer to the local node
      Node* const local_node_pt = node_pt(n);

      // Loop over the number of values stored at the node
      for (unsigned i = 0; i < n_value; i++)
      {
        // If the nodal value is not hanging
        if (local_node_pt->is_hanging(i) == false)
        {
          // Get the equation number
          local_unknown = nodal_local_eqn(n, i);
          // If the variable is free
          if (local_unknown >= 0)
          {
            // Get a pointer to the nodal value
            double* const value_pt = local_node_pt->value_pt(i);

            // Save the old value of nodal value
            const double old_var = *value_pt;

            // Increment the nodal value
            *value_pt += fd_step;

            // Now update any dependent variables
            update_in_nodal_fd(i);

            // Calculate the new residuals
            get_residuals(newres);

            //          if (use_first_order_fd)
            {
              // Do forward finite differences
              for (unsigned m = 0; m < n_dof; m++)
              {
                // Stick the entry into the Jacobian matrix
                jacobian(m, local_unknown) =
                  (newres[m] - residuals[m]) / fd_step;
              }
            }
            //           else
            //            {
            //             //Take backwards step
            //             *value_pt = old_var - fd_step;

            //             //Calculate the new residuals at backward position
            //             get_residuals(newres_minus);

            //             //Do central finite differences
            //             for(unsigned m=0;m<n_dof;m++)
            //              {
            //               //Stick the entry into the Jacobian matrix
            //               jacobian(m,local_unknown) =
            //                (newres[m] - newres_minus[m])/(2.0*fd_step);
            //              }
            //            }

            // Reset the nodal value
            *value_pt = old_var;

            // Reset any dependent variables
            reset_in_nodal_fd(i);
          }
        }
        // Otherwise the value is hanging node
        else
        {
          // Get the local hanging infor
          HangInfo* hang_info_pt = local_node_pt->hanging_pt(i);
          // Loop over the master nodes
          const unsigned n_master = hang_info_pt->nmaster();
          for (unsigned m = 0; m < n_master; m++)
          {
            // Get the pointer to the master node
            Node* const master_node_pt = hang_info_pt->master_node_pt(m);

            // Get the number of the unknown
            local_unknown = local_hang_eqn(master_node_pt, i);
            // If the variable is free
            if (local_unknown >= 0)
            {
              // Get a pointer to the nodal value stored at the master node
              double* const value_pt = master_node_pt->value_pt(i);

              // Save the old value of the nodal value stored at the master node
              const double old_var = *value_pt;

              // Increment the nodal value stored at the master node
              *value_pt += fd_step;

              // Now update any dependent variables
              update_in_nodal_fd(i);

              // Calculate the new residuals
              get_residuals(newres);

              //            if (use_first_order_fd)
              {
                // Do forward finite differences
                for (unsigned m = 0; m < n_dof; m++)
                {
                  // Stick the entry into the Jacobian matrix
                  jacobian(m, local_unknown) =
                    (newres[m] - residuals[m]) / fd_step;
                }
              }
              //             else
              //              {
              //               //Take backwards step
              //               *value_pt = old_var - fd_step;

              //               //Calculate the new residuals at backward
              //               position get_residuals(newres_minus);

              //               //Do central finite differences
              //               for(unsigned m=0;m<n_dof;m++)
              //                {
              //                 //Stick the entry into the Jacobian matrix
              //                 jacobian(m,local_unknown) =
              //                  (newres[m] - newres_minus[m])/(2.0*fd_step);
              //                }
              //              }

              // Reset the value at the master node
              *value_pt = old_var;

              // Reset any dependent variables
              reset_in_nodal_fd(i);
            }
          } // End of loop over master nodes
        } // End of hanging node case
      } // End of loop over values
    } // End of loop over nodes

    // End of finite difference loop
    // Final reset of any dependent data
    reset_after_nodal_fd();
  }


  //=========================================================================
  /// Internal function that is used to assemble the jacobian of the mapping
  /// from local coordinates (s) to the lagrangian coordinates (xi), given the
  /// derivatives of the shape functions.
  //=========================================================================
  void RefineableSolidElement::assemble_local_to_lagrangian_jacobian(
    const DShape& dpsids, DenseMatrix<double>& jacobian) const
  {
    // Find the the dimension of the element
    const unsigned el_dim = dim();
    // Find the number of shape functions and shape functions types
    const unsigned n_shape = nnode();
    const unsigned n_shape_type = nnodal_lagrangian_type();

    // Loop over the rows of the jacobian
    for (unsigned i = 0; i < el_dim; i++)
    {
      // Loop over the columns of the jacobian
      for (unsigned j = 0; j < el_dim; j++)
      {
        // Zero the entry
        jacobian(i, j) = 0.0;
        // Loop over the shape functions
        for (unsigned l = 0; l < n_shape; l++)
        {
          for (unsigned k = 0; k < n_shape_type; k++)
          {
            // Jacobian is dx_j/ds_i, which is represented by the sum
            // over the dpsi/ds_i of the nodal points X j
            // Use the hanging version here
            jacobian(i, j) +=
              lagrangian_position_gen(l, k, j) * dpsids(l, k, i);
          }
        }
      }
    }
  }

  //=========================================================================
  /// Internal function that is used to assemble the jacobian of second
  /// derivatives of the the mapping from local coordinates (s) to the
  /// lagrangian coordinates (xi), given the second derivatives of the
  /// shape functions.
  //=========================================================================
  void RefineableSolidElement::assemble_local_to_lagrangian_jacobian2(
    const DShape& d2psids, DenseMatrix<double>& jacobian2) const
  {
    // Find the the dimension of the element
    const unsigned el_dim = dim();
    // Find the number of shape functions and shape functions types
    const unsigned n_shape = nnode();
    const unsigned n_shape_type = nnodal_lagrangian_type();
    // Find the number of second derivatives
    const unsigned n_row = N2deriv[el_dim];

    // Assemble the "jacobian" (d^2 x_j/ds_i^2) for second derivatives of
    // shape functions
    // Loop over the rows (number of second derivatives)
    for (unsigned i = 0; i < n_row; i++)
    {
      // Loop over the columns (element dimension
      for (unsigned j = 0; j < el_dim; j++)
      {
        // Zero the entry
        jacobian2(i, j) = 0.0;
        // Loop over the shape functions
        for (unsigned l = 0; l < n_shape; l++)
        {
          // Loop over the shape function types
          for (unsigned k = 0; k < n_shape_type; k++)
          {
            // Add the terms to the jacobian entry
            // Use the hanging version here
            jacobian2(i, j) +=
              lagrangian_position_gen(l, k, j) * d2psids(l, k, i);
          }
        }
      }
    }
  }

  //==========================================================================
  /// Calculate the mapping from local to lagrangian coordinates
  /// assuming that the coordinates are aligned in the direction of the local
  /// coordinates, i.e. there are no cross terms and the jacobian is diagonal.
  /// The local derivatives are passed as dpsids and the jacobian and
  /// inverse jacobian are returned.
  //==========================================================================
  double RefineableSolidElement::local_to_lagrangian_mapping_diagonal(
    const DShape& dpsids,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& inverse_jacobian) const
  {
    // Find the dimension of the element
    const unsigned el_dim = dim();
    // Find the number of shape functions and shape functions types
    const unsigned n_shape = nnode();
    const unsigned n_shape_type = nnodal_lagrangian_type();

    // In this case we assume that there are no cross terms, that is
    // global coordinate 0 is always in the direction of local coordinate 0

    // Loop over the coordinates
    for (unsigned i = 0; i < el_dim; i++)
    {
      // Zero the jacobian and inverse jacobian entries
      for (unsigned j = 0; j < el_dim; j++)
      {
        jacobian(i, j) = 0.0;
        inverse_jacobian(i, j) = 0.0;
      }

      // Loop over the shape functions
      for (unsigned l = 0; l < n_shape; l++)
      {
        // Loop over the types of dof
        for (unsigned k = 0; k < n_shape_type; k++)
        {
          // Derivatives are always dx_{i}/ds_{i}
          jacobian(i, i) += lagrangian_position_gen(l, k, i) * dpsids(l, k, i);
        }
      }
    }

    // Now calculate the determinant of the matrix
    double det = 1.0;
    for (unsigned i = 0; i < el_dim; i++)
    {
      det *= jacobian(i, i);
    }

// Report if Matrix is singular, or negative
#ifdef PARANOID
    check_jacobian(det);
#endif

    // Calculate the inverse mapping (trivial in this case)
    for (unsigned i = 0; i < el_dim; i++)
    {
      inverse_jacobian(i, i) = 1.0 / jacobian(i, i);
    }

    // Return the value of the Jacobian
    return (det);
  }


  //========================================================================
  /// The number of geometric data affecting a
  /// RefineableSolidFiniteElement is the positional Data of all
  /// non-hanging nodes plus the geometric Data of all distinct
  /// master nodes. Recomputed on the fly.
  //========================================================================
  unsigned RefineableSolidElement::ngeom_data() const
  {
    // Find the number of nodes
    const unsigned n_node = nnode();

    // Temporary storage for unique position data
    std::set<Data*> all_position_data_pt;

    // Now loop over all the nodes again to find the master nodes
    // of any hanging nodes that have not yet been assigned
    for (unsigned n = 0; n < n_node; n++)
    {
      // If the node is a hanging node
      if (node_pt(n)->is_hanging())
      {
        // Find the local hang info object
        HangInfo* hang_info_pt = node_pt(n)->hanging_pt();

        // Find the number of master nodes
        unsigned n_master = hang_info_pt->nmaster();

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // Get the m-th master node
          Node* Master_node_pt = hang_info_pt->master_node_pt(m);

          // Add to set
          all_position_data_pt.insert(
            dynamic_cast<SolidNode*>(Master_node_pt)->variable_position_pt());
        }
      }
      // Not hanging
      else
      {
        // Add node itself to set
        all_position_data_pt.insert(
          dynamic_cast<SolidNode*>(node_pt(n))->variable_position_pt());
      }

    } // End of loop over nodes

    // How many are there?
    return all_position_data_pt.size();
  }


  //========================================================================
  /// \short Return pointer to the j-th Data item that the object's
  /// shape depends on: Positional data of non-hanging nodes and
  /// positional data of master nodes. Recomputed on the fly.
  //========================================================================
  Data* RefineableSolidElement::geom_data_pt(const unsigned& j)
  {
    // Find the number of nodes
    const unsigned n_node = nnode();

    // Temporary storage for unique position data. Set and vector are
    // required to ensure uniqueness in enumeration on different processors.
    // Set checks uniqueness; Vector stores entries in predictable order.
    std::set<Data*> all_position_data_pt;
    Vector<Data*> all_position_data_vector_pt;

    // Number of entries in set before possibly adding new entry
    unsigned n_old = 0;

    // Now loop over all the nodes again to find the master nodes
    // of any hanging nodes that have not yet been assigned
    for (unsigned n = 0; n < n_node; n++)
    {
      // If the node is a hanging node
      if (node_pt(n)->is_hanging())
      {
        // Find the local hang info object
        HangInfo* hang_info_pt = node_pt(n)->hanging_pt();

        // Find the number of master nodes
        unsigned n_master = hang_info_pt->nmaster();

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // Get the m-th master node
          Node* Master_node_pt = hang_info_pt->master_node_pt(m);

          // Positional data
          Data* pos_data_pt =
            dynamic_cast<SolidNode*>(Master_node_pt)->variable_position_pt();

          // Add to set
          n_old = all_position_data_pt.size();
          all_position_data_pt.insert(pos_data_pt);

          // New entry?
          if (all_position_data_pt.size() > n_old)
          {
            all_position_data_vector_pt.push_back(pos_data_pt);
          }
        }
      }
      // Not hanging
      else
      {
        // Add node itself to set

        // Positional data
        Data* pos_data_pt =
          dynamic_cast<SolidNode*>(node_pt(n))->variable_position_pt();

        // Add to set
        n_old = all_position_data_pt.size();
        all_position_data_pt.insert(pos_data_pt);

        // New entry?
        if (all_position_data_pt.size() > n_old)
        {
          all_position_data_vector_pt.push_back(pos_data_pt);
        }
      }

    } // End of loop over nodes


    // Return j-th entry
    return all_position_data_vector_pt[j];
  }


  //========================================================================
  /// \short Specify Data that affects the geometry of the element
  /// by adding the position Data to the set that's passed in.
  /// (This functionality is required in FSI problems; set is used to
  /// avoid double counting). Refineable version includes hanging nodes
  //========================================================================
  void RefineableSolidElement::identify_geometric_data(
    std::set<Data*>& geometric_data_pt)
  {
    // Loop over the node update data and add to the set
    const unsigned n_node = this->nnode();
    for (unsigned j = 0; j < n_node; j++)
    {
      // If the node is a hanging node
      if (node_pt(j)->is_hanging())
      {
        // Find the local hang info object
        HangInfo* hang_info_pt = node_pt(j)->hanging_pt();

        // Find the number of master nodes
        unsigned n_master = hang_info_pt->nmaster();

        // Loop over the master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
          // Get the m-th master node
          Node* Master_node_pt = hang_info_pt->master_node_pt(m);

          // Add to set
          geometric_data_pt.insert(
            dynamic_cast<SolidNode*>(Master_node_pt)->variable_position_pt());
        }
      }
      // Not hanging
      else
      {
        // Add node itself to set
        geometric_data_pt.insert(
          dynamic_cast<SolidNode*>(node_pt(j))->variable_position_pt());
      }
    }
  }

  //========================================================================
  /// The standard equation numbering scheme for solid positions,
  /// so that hanging Node information is included.
  //========================================================================
  void RefineableSolidElement::assign_solid_hanging_local_eqn_numbers(
    const bool& store_local_dof_pt)
  {
    // Find the number of nodes
    const unsigned n_node = nnode();

    // Check there are nodes!
    if (n_node > 0)
    {
      // Wipe the local matrix maps, they will be assigned below
      Local_position_hang_eqn.clear();

      // Find the local numbers
      const unsigned n_position_type = nnodal_position_type();
      const unsigned nodal_dim = nodal_dimension();

      // Matrix structure to store all positional equations at a node
      DenseMatrix<int> Position_local_eqn_at_node(n_position_type, nodal_dim);

      // Map that store whether the node's equation numbers have already been
      // added to the local arrays
      std::map<Node*, bool> local_eqn_number_done;

      // Get number of dofs so far
      unsigned local_eqn_number = ndof();

      // A local queue to store the global equation numbers
      std::deque<unsigned long> global_eqn_number_queue;

      // Now loop over all the nodes again to find the master nodes
      // of any hanging nodes that have not yet been assigned
      for (unsigned n = 0; n < n_node; n++)
      {
        // POSITIONAL EQUATAIONS
        // If the node is a hanging node
        if (node_pt(n)->is_hanging())
        {
          // Find the local hang info object
          HangInfo* hang_info_pt = node_pt(n)->hanging_pt();
          // Find the number of master nodes
          unsigned n_master = hang_info_pt->nmaster();
          // Loop over the master nodes
          for (unsigned m = 0; m < n_master; m++)
          {
            // Get the m-th master node
            Node* Master_node_pt = hang_info_pt->master_node_pt(m);

            // If the local equation numbers associated with this master node
            // have not already been assigned, assign them
            if (local_eqn_number_done[Master_node_pt] == false)
            {
              // Now we need to test whether the master node is actually
              // a local node, in which case its local equation numbers
              // will already have been assigned and stored in
              // position_local_eqn(n,j,k)

              // Storage for the index of the local node
              // Initialised to n_node (beyond the possible indices of
              //"real" nodes)
              unsigned local_node_index = n_node;
              // Loop over the local nodes (again)
              for (unsigned n1 = 0; n1 < n_node; n1++)
              {
                // If the master node is a local node
                // get its index and break out of the loop
                if (Master_node_pt == node_pt(n1))
                {
                  local_node_index = n1;
                  break;
                }
              }

              // Now we test whether the node was found
              if (local_node_index < n_node)
              {
                // Loop over the number of position dofs
                for (unsigned j = 0; j < n_position_type; j++)
                {
                  // Loop over the dimension of each node
                  for (unsigned k = 0; k < nodal_dim; k++)
                  {
                    // Set the values in the node-based positional look-up
                    // scheme
                    Position_local_eqn_at_node(j, k) =
                      position_local_eqn(local_node_index, j, k);
                  }
                }
              }
              // Otherwise it's a new master node
              else
              {
                // Loop over the number of position dofs
                for (unsigned j = 0; j < n_position_type; j++)
                {
                  // Loop over the dimension of each node
                  for (unsigned k = 0; k < nodal_dim; k++)
                  {
                    // Get equation number (position_eqn_number)
                    // Note eqn_number is long !
                    long eqn_number = static_cast<SolidNode*>(Master_node_pt)
                                        ->position_eqn_number(j, k);
                    // If equation_number positive add to array
                    if (eqn_number >= 0)
                    {
                      // Add global equation number to the local queue
                      global_eqn_number_queue.push_back(eqn_number);
                      // Add pointer to the dof to the queue if required
                      if (store_local_dof_pt)
                      {
                        GeneralisedElement::Dof_pt_deque.push_back(
                          &(Master_node_pt->x_gen(j, k)));
                      }
                      // Add to pointer-based scheme
                      Position_local_eqn_at_node(j, k) = local_eqn_number;
                      // Increase the number of local variables
                      local_eqn_number++;
                    }
                    // Otherwise the value is pinned
                    else
                    {
                      Position_local_eqn_at_node(j, k) = Data::Is_pinned;
                    }
                  }
                }
              } // End of case when it's a new master node

              // Dofs included with this node have now been done
              local_eqn_number_done[Master_node_pt] = true;
              // Add to the pointer-based reference scheme
              Local_position_hang_eqn[Master_node_pt] =
                Position_local_eqn_at_node;
            }
          }
        }

      } // End of loop over nodes

      // Now add our global equations numbers to the internal element storage
      add_global_eqn_numbers(global_eqn_number_queue,
                             GeneralisedElement::Dof_pt_deque);
      // Clear the memory used in the deque
      if (store_local_dof_pt)
      {
        std::deque<double*>().swap(GeneralisedElement::Dof_pt_deque);
      }


    } // End of if nodes
  }

  //============================================================================
  /// This function calculates the entries of Jacobian matrix, used in
  /// the Newton method, associated with the elastic problem in which the
  /// nodal position is a variable. It does this using finite differences,
  /// rather than an analytical formulation, so can be done in total generality.
  /// Overload the standard case to include hanging node case
  //==========================================================================
  void RefineableSolidElement::fill_in_jacobian_from_solid_position_by_fd(
    Vector<double>& residuals, DenseMatrix<double>& jacobian)
  {
    // Find the number of nodes
    const unsigned n_node = nnode();

    // If there are no nodes, return straight away
    if (n_node == 0)
    {
      return;
    }

    // Call the update function to ensure that the element is in
    // a consistent state before finite differencing starts
    update_before_solid_position_fd();

    //  bool use_first_order_fd=false;

    // Find the number of positional dofs and nodal dimension
    const unsigned n_position_type = nnodal_position_type();
    const unsigned nodal_dim = nodal_dimension();

    // Find the number of dofs in the element
    const unsigned n_dof = ndof();

    // Create newres vector
    Vector<double> newres(n_dof); //, newres_minus(n_dof);

    // Used default value defined in GeneralisedElement
    const double fd_step = Default_fd_jacobian_step;

    // Integer storage for local unknowns
    int local_unknown = 0;

    // Loop over the nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      // Get the pointer to the node
      Node* const local_node_pt = node_pt(l);

      // If the node is not a hanging node
      if (local_node_pt->is_hanging() == false)
      {
        // Loop over position dofs
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Loop over dimension
          for (unsigned i = 0; i < nodal_dim; i++)
          {
            local_unknown = position_local_eqn(l, k, i);
            // If the variable is free
            if (local_unknown >= 0)
            {
              // Store a pointer to the (generalised) Eulerian nodal position
              double* const value_pt = &(local_node_pt->x_gen(k, i));

              // Save the old value of the (generalised) Eulerian nodal position
              const double old_var = *value_pt;

              // Increment the  (generalised) Eulerian nodal position
              *value_pt += fd_step;

              // Perform any auxialiary node updates
              local_node_pt->perform_auxiliary_node_update_fct();

              // Update any other dependent variables
              update_in_solid_position_fd(l);


              // Calculate the new residuals
              get_residuals(newres);

              //            if (use_first_order_fd)
              {
                // Do forward finite differences
                for (unsigned m = 0; m < n_dof; m++)
                {
                  // Stick the entry into the Jacobian matrix
                  jacobian(m, local_unknown) =
                    (newres[m] - residuals[m]) / fd_step;
                }
              }
              //             else
              //              {
              //               //Take backwards step for the  (generalised)
              //               Eulerian nodal
              //               // position
              //               node_pt(l)->x_gen(k,i) = old_var-fd_step;

              //               //Calculate the new residuals at backward
              //               position get_residuals(newres_minus);

              //               //Do central finite differences
              //               for(unsigned m=0;m<n_dof;m++)
              //                {
              //                 //Stick the entry into the Jacobian matrix
              //                 jacobian(m,local_unknown) =
              //                  (newres[m] - newres_minus[m])/(2.0*fd_step);
              //                }
              //              }

              // Reset the (generalised) Eulerian nodal position
              *value_pt = old_var;

              // Perform any auxialiary node updates
              local_node_pt->perform_auxiliary_node_update_fct();

              // Reset any other dependent variables
              reset_in_solid_position_fd(l);
            }
          }
        }
      }
      // Otherwise it's a hanging node
      else
      {
        // Find the local hanging object
        HangInfo* hang_info_pt = local_node_pt->hanging_pt();
        // Loop over the master nodes
        const unsigned n_master = hang_info_pt->nmaster();
        for (unsigned m = 0; m < n_master; m++)
        {
          // Get the pointer to the master node
          Node* const master_node_pt = hang_info_pt->master_node_pt(m);

          // Get the local equation numbers for the master node
          DenseMatrix<int> Position_local_eqn_at_node =
            Local_position_hang_eqn[master_node_pt];

          // Loop over position dofs
          for (unsigned k = 0; k < n_position_type; k++)
          {
            // Loop over dimension
            for (unsigned i = 0; i < nodal_dim; i++)
            {
              local_unknown = Position_local_eqn_at_node(k, i);
              // If the variable is free
              if (local_unknown >= 0)
              {
                // Store a pointer to the (generalised) Eulerian nodal position
                double* const value_pt = &(master_node_pt->x_gen(k, i));

                // Save the old value of the (generalised) Eulerian nodal
                // position
                const double old_var = *value_pt;

                // Increment the  (generalised) Eulerian nodal position
                *value_pt += fd_step;

                // Perform any auxialiary node updates
                master_node_pt->perform_auxiliary_node_update_fct();

                // Update any dependent variables
                update_in_solid_position_fd(l);

                // Calculate the new residuals
                get_residuals(newres);

                //              if (use_first_order_fd)
                {
                  // Do forward finite differences
                  for (unsigned m = 0; m < n_dof; m++)
                  {
                    // Stick the entry into the Jacobian matrix
                    jacobian(m, local_unknown) =
                      (newres[m] - residuals[m]) / fd_step;
                  }
                }
                //               else
                //                {
                //                 //Take backwards step for the  (generalised)
                //                 Eulerian nodal
                //                 // position
                //                 master_node_pt->x_gen(k,i) = old_var-fd_step;

                //                 //Calculate the new residuals at backward
                //                 position get_residuals(newres_minus);

                //                 //Do central finite differences
                //                 for(unsigned m=0;m<n_dof;m++)
                //                  {
                //                   //Stick the entry into the Jacobian matrix
                //                   jacobian(m,local_unknown) =
                //                    (newres[m] -
                //                    newres_minus[m])/(2.0*fd_step);
                //                  }
                //                }

                // Reset the (generalised) Eulerian nodal position
                *value_pt = old_var;

                // Perform any auxialiary node updates
                master_node_pt->perform_auxiliary_node_update_fct();

                // Reset any other dependent variables
                reset_in_solid_position_fd(l);
              }
            }
          }
        }
      } // End of hanging node case

    } // End of loop over nodes

    // End of finite difference loop
    // Final reset of any dependent data
    reset_after_solid_position_fd();
  }

} // namespace oomph
