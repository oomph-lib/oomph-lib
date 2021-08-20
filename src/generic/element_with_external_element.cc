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
// Functions for the ElementWithExternalElement class

#include "element_with_external_element.h"

namespace oomph
{
  //======================================================================
  /// \short Destructor, clean up any memory allocated for the equation
  /// numbering schemes and storage for the external elements
  //======================================================================
  ElementWithExternalElement::~ElementWithExternalElement()
  {
    // If storage has been allocated for external geometric data,
    // delete it
    if (External_interaction_geometric_data_local_eqn)
    {
      delete[] External_interaction_geometric_data_local_eqn;
    }
    if (External_interaction_geometric_data_index)
    {
      delete[] External_interaction_geometric_data_index;
    }
    if (External_interaction_geometric_data_pt)
    {
      delete[] External_interaction_geometric_data_pt;
    }

    // If storage has been allocated for external field data,
    // delete it
    if (External_interaction_field_data_local_eqn)
    {
      delete[] External_interaction_field_data_local_eqn;
    }
    if (External_interaction_field_data_index)
    {
      delete[] External_interaction_field_data_index;
    }
    if (External_interaction_field_data_pt)
    {
      delete[] External_interaction_field_data_pt;
    }

    // Flush the storage associated with the external elements
    this->flush_all_external_element_storage();
  }

  //========================================================================
  /// \short Initialise storage for source elements and their
  /// associated local coordinates for the specified number of interactions
  //========================================================================
  void ElementWithExternalElement::initialise_external_element_storage()
  {
    // Find the number of interactions
    const unsigned n_interaction = this->Ninteraction;
    // Find the number of integration points
    const unsigned n_intpt = this->integral_pt()->nweight();

    // Work out the required storage
    const unsigned n_external_element_storage = n_interaction * n_intpt;

    // If we have not allocated the correct amount of memory,
    // then do so
    if (Nexternal_element_storage != n_external_element_storage)
    {
      // Allocate storage for the pointers to the elements
      // Delete old storage, if needed
      if (External_element_pt)
      {
        delete[] External_element_pt;
      }
      // Allocate new memory
      External_element_pt = new FiniteElement*[n_external_element_storage];

      // Initialise all new pointers to zero
      for (unsigned i = 0; i < n_external_element_storage; i++)
      {
        External_element_pt[i] = 0;
      }

      // Alloacte storage for local coordinates
      // Delete old storage, if needed
      if (External_element_local_coord)
      {
        delete[] External_element_local_coord;
      }
      // Allocate the new memory
      External_element_local_coord =
        new Vector<double>[n_external_element_storage];

      // Finally record how much memory we have allocated
      Nexternal_element_storage = n_external_element_storage;
      Nintpt = n_intpt;
    }
  }

  //========================================================================
  /// \short Clear the storage for pointers to external  elements and their
  /// associated local coordinates.
  //========================================================================
  void ElementWithExternalElement::flush_all_external_element_storage()
  {
    // Delete the memory if it has been allocated
    if (External_element_pt)
    {
      delete[] External_element_pt;
      External_element_pt = 0;
    }

    if (External_element_local_coord)
    {
      delete[] External_element_local_coord;
      External_element_local_coord = 0;
    }

    // Reset the number of stored values to zero
    Nexternal_element_storage = 0;
  }


  //================================================================
  /// \short Function that must return all the data involved
  /// in the desired interactions from the external element
  /// Default is to call the identify_field_data_for_interaction()
  /// for each element.
  //================================================================
  void ElementWithExternalElement::
    identify_all_field_data_for_external_interaction(
      Vector<std::set<FiniteElement*>> const& external_elements_pt,
      std::set<std::pair<Data*, unsigned>>& paired_interaction_data)
  {
    // Loop over each interaction
    const unsigned n_interaction = this->ninteraction();
    for (unsigned i = 0; i < n_interaction; i++)
    {
      // Loop over each element in the set
      for (std::set<FiniteElement*>::const_iterator it =
             external_elements_pt[i].begin();
           it != external_elements_pt[i].end();
           it++)
      {
        (*it)->identify_field_data_for_interactions(paired_interaction_data);
      }
    } // End of loop over interactions
  }

  //=======================================================================
  /// \short Function that must return all geometric data involved
  /// in the desired interactions from the external element
  /// Default is to add all geometric data for each element for each
  /// interaction
  //=======================================================================
  void ElementWithExternalElement::
    identify_all_geometric_data_for_external_interaction(
      Vector<std::set<FiniteElement*>> const& external_elements_pt,
      std::set<Data*>& external_geometric_data_pt)
  {
    // Loop over each interaction
    const unsigned n_interaction = this->ninteraction();
    for (unsigned i = 0; i < n_interaction; i++)
    {
      // Loop over each element in the set
      for (std::set<FiniteElement*>::const_iterator it =
             external_elements_pt[i].begin();
           it != external_elements_pt[i].end();
           it++)
      {
        (*it)->identify_geometric_data(external_geometric_data_pt);
      }
    } // End of loop over interactions
  }

  /// \short Function to describe the local dofs of the element. The ostream
  /// specifies the output stream to which the description
  /// is written; the string stores the currently
  /// assembled output that is ultimately written to the
  /// output stream by Data::describe_dofs(...); it is typically
  /// built up incrementally as we descend through the
  /// call hierarchy of this function when called from
  /// Problem::describe_dofs(...)
  void ElementWithExternalElement::describe_local_dofs(
    std::ostream& out, const std::string& current_string) const
  {
    // Find the number of external field data
    const unsigned n_external_field_data = nexternal_interaction_field_data();
    // Now loop over the field data again to assign local equation numbers
    for (unsigned i = 0; i < n_external_field_data; i++)
    {
      std::stringstream conversion;
      conversion << " of External Interaction Field Data " << i
                 << current_string;
      std::string in(conversion.str());
      External_interaction_field_data_pt[i]->describe_dofs(out, in);
    }

    // Find the number of external geometric data
    unsigned n_external_geom_data = nexternal_interaction_geometric_data();

    // Now loop over the field data again assign local equation numbers
    for (unsigned i = 0; i < n_external_geom_data; i++)
    {
      std::stringstream conversion;
      conversion << " of External Interaction Geometric Data " << i
                 << current_string;
      std::string in(conversion.str());
      External_interaction_geometric_data_pt[i]->describe_dofs(out, in);
    }
    GeneralisedElement::describe_local_dofs(out, current_string);
  }

  //==========================================================================
  /// This function determines the all Data in external elemetns that
  /// affects the residuals of the element
  /// and adds their global equation numbers to the
  /// local-to-global look-up scheme. Note that we only include
  /// Data items into the element's External_interaction_data
  /// if they are not already
  /// included in the element's nodal positional Data, its internal
  /// or external Data.
  //==========================================================================
  void ElementWithExternalElement::
    assign_external_interaction_data_local_eqn_numbers(
      const bool& store_local_dof_pt)
  {
    // Reset number of stored field data to zero
    Nexternal_interaction_field_data = 0;
    // Clear all the internal field data storage, if it's been allocated
    if (External_interaction_field_data_pt)
    {
      delete[] External_interaction_field_data_pt;
      External_interaction_field_data_pt = 0;
    }
    if (External_interaction_field_data_index)
    {
      delete[] External_interaction_field_data_index;
      External_interaction_field_data_index = 0;
    }
    if (External_interaction_field_data_local_eqn)
    {
      delete[] External_interaction_field_data_local_eqn;
      External_interaction_field_data_local_eqn = 0;
    }

    // Reset number of stored geometric data to zero
    Nexternal_interaction_geometric_data = 0;
    // Clear all internal external data storage, if it's been allocated
    if (External_interaction_geometric_data_pt)
    {
      delete[] External_interaction_geometric_data_pt;
      External_interaction_geometric_data_pt = 0;
    }
    if (External_interaction_geometric_data_index)
    {
      delete[] External_interaction_geometric_data_index;
      External_interaction_geometric_data_index = 0;
    }
    if (External_interaction_geometric_data_local_eqn)
    {
      delete[] External_interaction_geometric_data_local_eqn;
      External_interaction_geometric_data_local_eqn = 0;
    }

    // Only bother with non-halo elements
#ifdef OOMPH_HAS_MPI
    if (!this->is_halo())
#endif
    {
      // If desired, determine the Data that affects the interactions but is not
      // already included in elements other generic Data.
      // The conditional test is done here, so that the vectors are cleared
      // and not re-filled if Add_external_interaction_data is false
      if (Add_external_interaction_data)
      {
        // Number of interactions
        const unsigned n_interaction = this->ninteraction();
        // Number of integration points
        const unsigned n_intpt = integral_pt()->nweight();

        // Sets of all (external) FiniteElements that affect the interactions
        // One set per interaction
        Vector<std::set<FiniteElement*>> external_interaction_elements_pt(
          n_interaction);

        // Loop over the interactions
        for (unsigned i = 0; i < n_interaction; i++)
        {
          // Loop over the integration points and the adjacent element at
          // each integration point to the set
          for (unsigned ipt = 0; ipt < n_intpt; ipt++)
          {
            // Add the element adjacent to the element into the set
            external_interaction_elements_pt[i].insert(
              external_element_pt(i, ipt));
          }
          // For safety erase any null pointers
          external_interaction_elements_pt[i].erase(0);
        }

        // Storage for a pairs of interaction data (pointer to Data and the
        // index of the value within this Data object) affecting the element.
        std::set<std::pair<Data*, unsigned>> paired_field_data;

        // Determine the field data that affects the external interactions
        // for all sets of external elements
        identify_all_field_data_for_external_interaction(
          external_interaction_elements_pt, paired_field_data);

        // It's just possible that some of the field data could be internal
        // or nodal data, so we should remove it if that's the case
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
            paired_field_data.erase(std::make_pair(dat_pt, i));
          }
        }

        // Loop over all the nodes
        const unsigned n_node = this->nnode();
        for (unsigned n = 0; n < n_node; n++)
        {
          // Find the node point
          Node* const nod_pt = this->node_pt(n);
          // Find the number of values stored at the node
          const unsigned n_value = nod_pt->nvalue();
          // Loop over values and erase all pairs from the set
          for (unsigned i = 0; i < n_value; i++)
          {
            paired_field_data.erase(std::make_pair(nod_pt, i));
          }
          SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
          if (solid_nod_pt != 0)
          {
            Data* pos_data_pt = solid_nod_pt->variable_position_pt();
            // Find the number of positional values stored at the node
            const unsigned n_value = pos_data_pt->nvalue();
            // Loop over values and erase all pairs from the set
            for (unsigned i = 0; i < n_value; i++)
            {
              paired_field_data.erase(std::make_pair(pos_data_pt, i));
            }
          }
        }

        // Now allocate storage for the external field data
        // associated indices and local equation numbers
        const unsigned n_external_interaction_field_data =
          paired_field_data.size();
        Nexternal_interaction_field_data = n_external_interaction_field_data;
        External_interaction_field_data_pt =
          new Data*[n_external_interaction_field_data];
        External_interaction_field_data_index =
          new unsigned[n_external_interaction_field_data];

        // Add the pairs of data to the field data vectors
        {
          unsigned count = 0;
          for (std::set<std::pair<Data*, unsigned>>::iterator it =
                 paired_field_data.begin();
               it != paired_field_data.end();
               it++)
          {
            External_interaction_field_data_pt[count] = it->first;
            External_interaction_field_data_index[count] = it->second;
            ++count;
          }
        }


        // Only bother to add geometric data if we're told to
        if (Add_external_geometric_data)
        {
          // Storage for a set of external geometric Data affecting the element
          std::set<Data*> external_geometric_data_pt;

          // Determine the geometric data that affects the external interactions
          // for all sets of external elements
          identify_all_geometric_data_for_external_interaction(
            external_interaction_elements_pt, external_geometric_data_pt);

          // Now loop over any geometric data of the Element itself
          // and erase them from the external_geometric_data_pt
          // because these data are actually intrinsic
          // data of this element and are counted and numbered elsewhere
          unsigned n_geom_data = ngeom_data();
          for (unsigned j = 0; j < n_geom_data; j++)
          {
            external_geometric_data_pt.erase(geom_data_pt(j));
          }

          // It is possible that the geometric data may have already been added
          // as external data. We should erase any common entries from the
          // External_interaction_data
          // but not touch the external data that has been set up by a
          //"knowledgeable" user
          unsigned n_external = nexternal_data();
          for (unsigned j = 0; j < n_external; j++)
          {
            external_geometric_data_pt.erase(external_data_pt(j));
          }

          // Loop over all the nodes of present element to avoid double counting
          const unsigned n_node = this->nnode();
          for (unsigned n = 0; n < n_node; n++)
          {
            Node* const nod_pt = this->node_pt(n);
            external_geometric_data_pt.erase(nod_pt);

            SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
            if (solid_nod_pt != 0)
            {
              Data* pos_data_pt = solid_nod_pt->variable_position_pt();
              //          std::ostringstream junk;
              //          junk << "Erasing ";
              //          unsigned nval=pos_data_pt->nvalue();
              //          for (unsigned i=0;i<nval;i++)
              //           {
              //            junk << pos_data_pt->eqn_number(i) << " ";
              //           }
              //          oomph_info << junk.str() << std::endl;
              external_geometric_data_pt.erase(pos_data_pt);
            }
          }


          // Next allocate storage for the geometric field data
          // Find out how many individual data we have
          unsigned n_external_interaction_geometric_data = 0;
          for (std::set<Data*>::iterator it =
                 external_geometric_data_pt.begin();
               it != external_geometric_data_pt.end();
               it++)
          {
            // Add the number of values stored in each geometric datum
            n_external_interaction_geometric_data += (*it)->nvalue();
          }

          // Now allocate storage
          Nexternal_interaction_geometric_data =
            n_external_interaction_geometric_data;
          External_interaction_geometric_data_pt =
            new Data*[n_external_interaction_geometric_data];
          External_interaction_geometric_data_index =
            new unsigned[n_external_interaction_geometric_data];

          // Now we can add all the geometric data to the geometric data vectors
          {
            unsigned count = 0;
            for (std::set<Data*>::iterator it =
                   external_geometric_data_pt.begin();
                 it != external_geometric_data_pt.end();
                 it++)
            {
              // Find the number of values stored in the geometric data
              unsigned n_value = (*it)->nvalue();
              // Loop over the values
              for (unsigned j = 0; j < n_value; j++)
              {
                // Add data to the external geometric data
                External_interaction_geometric_data_pt[count] = *it;
                External_interaction_geometric_data_index[count] = j;
                ++count;
              }
            }
          }
        }
      }

      // All external interaction data has now been specified

      // Find the number of external field data
      const unsigned n_external_field_data = nexternal_interaction_field_data();

      // If there are interaction data fill in the internal storage
      if (n_external_field_data > 0)
      {
        // Allocate storage for the local equation numbers associated with
        // external field data
        External_interaction_field_data_local_eqn =
          new int[n_external_field_data];

        // Find the number of local equations assigned so far
        unsigned local_eqn_number = ndof();

        // A local queue to store the global equation numbers
        std::deque<unsigned long> global_eqn_number_queue;

        // Now loop over the field data again to assign local equation numbers
        for (unsigned i = 0; i < n_external_field_data; i++)
        {
          // Get the GLOBAL equation number
          long eqn_number = External_interaction_field_data_pt[i]->eqn_number(
            External_interaction_field_data_index[i]);

          // If the GLOBAL equation number is positive (i.e. not pinned)
          if (eqn_number >= 0)
          {
            //        std::ostringstream junk;
            //        junk << "Adding global eqn " << eqn_number << " (";
            //        if (!(External_interaction_field_data_pt[i]->is_halo()))
            //         {
            //          junk << "not";
            //         }
            //        oomph_info << junk.str() << " halo" << std::endl;

            // Add the GLOBAL equation number to the local queue
            global_eqn_number_queue.push_back(eqn_number);
            // Add pointer to the dof to the queue if required
            if (store_local_dof_pt)
            {
              GeneralisedElement::Dof_pt_deque.push_back(
                External_interaction_field_data_pt[i]->value_pt(
                  External_interaction_field_data_index[i]));
            }

            // Add the local equation number to the local scheme
            External_interaction_field_data_local_eqn[i] = local_eqn_number;
            // Increase the local number
            local_eqn_number++;
          }
          else
          {
            // Set the local scheme to be pinned
            External_interaction_field_data_local_eqn[i] = Data::Is_pinned;
          }
        }
        // Now add our global equations numbers to the internal element storage
        add_global_eqn_numbers(global_eqn_number_queue,
                               GeneralisedElement::Dof_pt_deque);
        // Clear the memory used in the deque
        if (store_local_dof_pt)
        {
          std::deque<double*>().swap(GeneralisedElement::Dof_pt_deque);
        }
      }

      // Find the number of external geometric data
      unsigned n_external_geom_data = nexternal_interaction_geometric_data();

      // If there are external geometric data fill in the internal storage
      if (n_external_geom_data > 0)
      {
        // Allocate storage for the local equation numbers associated with
        // external geometric data
        External_interaction_geometric_data_local_eqn =
          new int[n_external_geom_data];

        // Find the number of local equations assigned so far
        unsigned local_eqn_number = ndof();

        // A local queue to store the global equation numbers
        std::deque<unsigned long> global_eqn_number_queue;

        // Now loop over the field data again assign local equation numbers
        for (unsigned i = 0; i < n_external_geom_data; i++)
        {
          // Get the GLOBAL equation number
          long eqn_number =
            External_interaction_geometric_data_pt[i]->eqn_number(
              External_interaction_geometric_data_index[i]);

          // If the GLOBAL equation number is positive (a free variable)
          if (eqn_number >= 0)
          {
            // Add the GLOBAL equation number to the local queue
            global_eqn_number_queue.push_back(eqn_number);
            // Add pointer to the dof to the queue if required
            if (store_local_dof_pt)
            {
              GeneralisedElement::Dof_pt_deque.push_back(
                External_interaction_geometric_data_pt[i]->value_pt(
                  External_interaction_geometric_data_index[i]));
            }

            // Add the local equation number to the local scheme
            External_interaction_geometric_data_local_eqn[i] = local_eqn_number;
            // Increase the local number
            local_eqn_number++;
          }
          else
          {
            // Set the local scheme to be pinned
            External_interaction_geometric_data_local_eqn[i] = Data::Is_pinned;
          }
        }
        // Now add our global equations numbers to the internal element storage
        add_global_eqn_numbers(global_eqn_number_queue,
                               GeneralisedElement::Dof_pt_deque);
        // Clear the memory used in the deque
        if (store_local_dof_pt)
        {
          std::deque<double*>().swap(GeneralisedElement::Dof_pt_deque);
        }
      }
    }
  }

  //============================================================================
  /// This function calculates the entries of Jacobian matrix, used in
  /// the Newton method, associated with the external interaction
  /// degrees of freedom for external fields.
  /// It does this using finite differences,
  /// rather than an analytical formulation, so can be done in total generality.
  //==========================================================================
  void ElementWithExternalElement::
    fill_in_jacobian_from_external_interaction_field_by_fd(
      Vector<double>& residuals, DenseMatrix<double>& jacobian)
  {
    // Locally cache the number of data
    const unsigned n_external_interaction_field_data =
      nexternal_interaction_field_data();

    // If there is no such data return
    if (n_external_interaction_field_data == 0)
    {
      return;
    }

    // Call the update function to ensure that the element is in a
    // consistent state before finite differencing
    update_before_external_interaction_field_fd();

    // Find the number of dofs in the element
    const unsigned n_dof = ndof();

    // Create newres vector
    Vector<double> newres(n_dof);

    // Integer storage for local unknown
    int local_unknown = 0;

    // Use the default finite difference step
    const double fd_step = Default_fd_jacobian_step;

    // Loop over the data
    for (unsigned i = 0; i < n_external_interaction_field_data; i++)
    {
      // Find the value of the local unknown
      local_unknown = External_interaction_field_data_local_eqn[i];
      // If it's not a boundary condition
      if (local_unknown >= 0)
      {
        // Store a pointer to the field value
        double* value_pt = External_interaction_field_data_pt[i]->value_pt(
          External_interaction_field_data_index[i]);

        // Save the old value of the field value
        double old_var = *value_pt;

        // Increment the value
        *value_pt += fd_step;

        // Now update any dependent variables
        update_in_external_interaction_field_fd(i);

        // Calculate the new residuals
        get_residuals(newres);

        // Do forward finite differences
        for (unsigned m = 0; m < n_dof; m++)
        {
          // Stick the entry into the Jacobian matrix
          jacobian(m, local_unknown) = (newres[m] - residuals[m]) / fd_step;
        }

        // Reset the variables
        *value_pt = old_var;

        // Reset any dependent variables
        reset_in_external_interaction_field_fd(i);
      }
    } // End of loop over external interaction data

    // End of finite difference loop
    // Final reset of any dependent data
    reset_after_external_interaction_field_fd();
  }


  //============================================================================
  /// This function calculates the entries of Jacobian matrix, used in
  /// the Newton method, associated with the external interaction
  /// degrees of freedom for external geometric data.
  /// It does this using finite differences,
  /// rather than an analytical formulation, so can be done in total generality.
  //==========================================================================
  void ElementWithExternalElement::
    fill_in_jacobian_from_external_interaction_geometric_by_fd(
      Vector<double>& residuals, DenseMatrix<double>& jacobian)
  {
    // Locally cache the number of data
    const unsigned n_external_interaction_geometric_data =
      nexternal_interaction_geometric_data();
    // If there is no such data return
    if (n_external_interaction_geometric_data == 0)
    {
      return;
    }

    // Call the update function to ensure that the element is in a
    // consistent state before finite differencing
    update_before_external_interaction_geometric_fd();

    // Find the number of dofs in the element
    const unsigned n_dof = ndof();

    // Create newres vector
    Vector<double> newres(n_dof);

    // Integer storage for local unknown
    int local_unknown = 0;

    // Use the default finite difference step
    const double fd_step = Default_fd_jacobian_step;

    // Loop over the data
    for (unsigned i = 0; i < n_external_interaction_geometric_data; i++)
    {
      // Find the value of the local unknown
      local_unknown = External_interaction_geometric_data_local_eqn[i];
      // If it's not a boundary condition
      if (local_unknown >= 0)
      {
        // Store a pointer to the geometric value
        double* value_pt = External_interaction_geometric_data_pt[i]->value_pt(
          External_interaction_geometric_data_index[i]);

        // Save the old value of the geometric value
        double old_var = *value_pt;

        // Increment the value
        *value_pt += fd_step;

        // Now update any dependent variables
        update_in_external_interaction_geometric_fd(i);

        // Calculate the new residuals
        get_residuals(newres);

        // Do forward finite differences
        for (unsigned m = 0; m < n_dof; m++)
        {
          // Stick the entry into the Jacobian matrix
          jacobian(m, local_unknown) = (newres[m] - residuals[m]) / fd_step;
        }

        // Reset the variables
        *value_pt = old_var;

        // Reset any dependent variables
        reset_in_external_interaction_geometric_fd(i);
      }
    } // End of loop over external interaction data

    // End of finite difference loop
    // Final reset of any dependent data
    reset_after_external_interaction_geometric_fd();
  }


  //============================================================================
  /// Output by plotting vector from integration point to
  /// corresponding point in external element for specified interaction
  /// index
  //==========================================================================
  void ElementWithExternalElement::output_external_elements(
    std::ostream& outfile, const unsigned& interaction_index)
  {
    // Dimension of element
    unsigned n_dim_el = dim();
    Vector<double> s(n_dim_el);

    // Vectors for coordintes
    unsigned n_dim = node_pt(0)->ndim();
    Vector<double> x(n_dim);
    Vector<double> x_ext(n_dim);

    // Loop over the integration points
    const unsigned n_intpt = this->integral_pt()->nweight();
    outfile << "ZONE I=" << n_intpt << std::endl;
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      for (unsigned i = 0; i < n_dim_el; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Eulerian coordinates of integration point
      interpolated_x(s, x);

      // Get pointer to external element
      FiniteElement* ext_el_pt = external_element_pt(interaction_index, ipt);
      // Get local coordinate in external element
      Vector<double> s_ext(
        external_element_local_coord(interaction_index, ipt));

      // Eulerian coordinates of point in external element
      ext_el_pt->interpolated_x(s_ext, x_ext);

      // Output coords of interation point
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << x[i] << " ";
      }
      // Write vector to point in external element
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << x_ext[i] - x[i] << " ";
      }
      outfile << std::endl;
    }
  }


} // namespace oomph
