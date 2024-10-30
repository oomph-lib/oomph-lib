// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
#ifndef MY_TRIANGLE_MESH_HEADER
#define MY_TRIANGLE_MESH_HEADER

#include "generic.h"
#include "meshes/triangle_mesh.h"

namespace oomph
{
  //============start_of_my_triangle_mesh_class===========================
  /// Upgrade TriangleMesh<ELEMENT> to be able to deep copy the mesh from
  /// another Mesh object.
  //======================================================================
  template<class ELEMENT>
  class MyTriangleMesh : public TriangleMesh<ELEMENT>
  {
  private:
    ErrorEstimator* Spatial_error_estimator_pt;

  public:
    /// Constructor with time stepper only
    MyTriangleMesh<ELEMENT>(
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Constructor with mesh parameters and/or time stepper
    MyTriangleMesh<ELEMENT>(
      TriangleMeshParameters& triangle_mesh_parameters,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Deep copy from passed mesh
    void deep_copy(Mesh* orig_mesh_pt);

    /// Access to spatial error estimator
    ErrorEstimator*& spatial_error_estimator_pt()
    {
      return Spatial_error_estimator_pt;
    }

    /// Access to spatial error estimator (const version
    ErrorEstimator* spatial_error_estimator_pt() const
    {
      return Spatial_error_estimator_pt;
    }
  };

  //======================================================================
  /// Constructor with no arguments.
  //======================================================================
  template<class ELEMENT>
  MyTriangleMesh<ELEMENT>::MyTriangleMesh(TimeStepper* time_stepper_pt)
    : TriangleMesh<ELEMENT>()
  {
    /// Set the time stepper to the default one as none is passed in.
    this->Time_stepper_pt = time_stepper_pt;
  }

  //======================================================================
  /// Constructor with mesh parameters and/or time stepper.
  /// Use the underlying TriangleMesh constructor.
  //======================================================================
  template<class ELEMENT>
  MyTriangleMesh<ELEMENT>::MyTriangleMesh(
    TriangleMeshParameters& triangle_mesh_parameters,
    TimeStepper* time_stepper_pt)
    : TriangleMesh<ELEMENT>(triangle_mesh_parameters, time_stepper_pt)
  {
  }

  //======================================================================
  /// Deep copy from passed mesh
  //======================================================================
  template<class ELEMENT>
  void MyTriangleMesh<ELEMENT>::deep_copy(Mesh* orig_mesh_pt)
  {
    unsigned nelem = orig_mesh_pt->nelement();
    this->Element_pt.resize(nelem);

    unsigned nnod = orig_mesh_pt->nnode();
    this->Node_pt.resize(nnod);

    // Create a map storing the node_id of the mesh used to update the
    // node position in the update_triangulateio function
    std::map<Node*, unsigned> old_global_number;

    // Store the TriangulateIO node id
    for (unsigned inod = 0; inod < nnod; inod++)
    {
      Node* old_node_pt = orig_mesh_pt->node_pt(inod);
      old_global_number[old_node_pt] = inod;
    }

    // Initialize the old node id vector
    this->Oomph_vertex_nodes_id.resize(nnod);

    unsigned nbound = orig_mesh_pt->nboundary();
    this->set_nboundary(nbound);

    for (unsigned e = 0; e < nelem; e++)
    {
      this->Element_pt[e] = new ELEMENT;
    }

    /// Number of nodes in each element. Assuming the same for all elements.
    unsigned nnod_el = orig_mesh_pt->finite_element_pt(0)->nnode();

    // Setup map to check the (pseudo-)global node number
    // Nodes whose number is zero haven't been copied across
    // into the mesh yet.
    std::map<Node*, unsigned> global_number;
    unsigned global_count = 0;

    /// Loop over elements in original mesh, visit their nodes
    for (unsigned e = 0; e < nelem; e++)
    {
      /// Loop over all nodes in element
      for (unsigned j = 0; j < nnod_el; j++)
      {
        // Pointer to node in the scaffold mesh
        Node* orig_node_pt = orig_mesh_pt->finite_element_pt(e)->node_pt(j);

        // Get the (pseudo-)global node number in scaffold mesh
        // (It's zero [=default] if not visited this one yet)
        unsigned j_global = global_number[orig_node_pt];

        // Haven't done this one yet
        if (j_global == 0)
        {
          // Find and store the node_id in the old nodes map
          this->Oomph_vertex_nodes_id[global_count] =
            old_global_number[orig_node_pt];

          // Get pointer to set of mesh boundaries that this
          // scaffold node occupies; NULL if the node is not on any boundary
          std::set<unsigned>* boundaries_pt;
          orig_node_pt->get_boundaries_pt(boundaries_pt);

          // Storage for the new node
          Node* new_node_pt = 0;

          // Is it on boundaries
          if (boundaries_pt != 0)
          {
            // Create new boundary node
            new_node_pt = this->finite_element_pt(e)->construct_boundary_node(
              j, this->Time_stepper_pt);

            // Add to boundaries
            for (std::set<unsigned>::iterator it = boundaries_pt->begin();
                 it != boundaries_pt->end();
                 ++it)
            {
              this->add_boundary_node(*it, new_node_pt);
            }
          }
          // build normal node
          else
          {
            /// Create new node, using the NEW element's construct_node
            /// member function
            new_node_pt = this->finite_element_pt(e)->construct_node(
              j, this->Time_stepper_pt);
          }

          // Give it a number (not necessarily the global node
          // number in the scaffold mesh -- we just need something
          // to keep track...)
          global_count++;
          global_number[orig_node_pt] = global_count;

          // Copy new node, created using the NEW element's construct_node
          // function into global storage, using the same global
          // node number that we've just associated with the
          // corresponding node in the scaffold mesh
          this->Node_pt[old_global_number[orig_node_pt]] = new_node_pt;

          // Assign coordinates
          for (unsigned i = 0; i < this->finite_element_pt(e)->dim(); i++)
          {
            new_node_pt->x(i) = orig_node_pt->x(i);
          }
        }
        // This one has already been done: Copy accross
        else
        {
          this->finite_element_pt(e)->node_pt(j) =
            this->Node_pt[old_global_number[orig_node_pt]];
        }
      }
    }

    // ============================================================
    /*/// Construct the elements

    /// Number of elements
    unsigned nelem = orig_mesh_pt->nelement();
    this->Element_pt.resize(nelem);
    for (unsigned e = 0; e < nelem; e++)
    {
      this->Element_pt[e] = new ELEMENT;
    }

    // Number of boundaries
    unsigned nbound = orig_mesh_pt->nboundary();
    this->set_nboundary(nbound);

    /// Construct the nodes

    /// Number of nodes
    unsigned nnod = orig_mesh_pt->nnode();
    this->Node_pt.resize(nnod);

    /// In the first instance build all nodes from within all the elements

    /// Number of nodes in each element. Assuming the same for all elements.
    unsigned nnod_el = orig_mesh_pt->finite_element_pt(0)->nnode();

    /// Loop over elements in original mesh, visit their nodes
    for (unsigned e = 0; e < nelem; e++)
    {
      /// Loop over all nodes in element
      for (unsigned j = 0; j < nnod_el; j++)
      {
        /// Create new node, using the NEW element's construct_node
        /// member function
        this->finite_element_pt(e)->construct_node(j, this->Time_stepper_pt);
      }
    }

    /// Setup nodes.
    /// Upgrade boundary nodes and assign to boundaries.
    /// Delete multiple nodes at the same point.

    /// A map from each node pointer to its assigned global count number.
    std::map<Node*, unsigned> global_number;

    /// Global counter to loop through each unique node
    unsigned global_count = 0;

    /// Loop over elements in original mesh, visit their nodes
    for (unsigned e = 0; e < nelem; e++)
    {
      /// Loop over all nodes in element
      for (unsigned j = 0; j < nnod_el; j++)
      {
        /// Pointer to current node in the original mesh
        Node* orig_node_pt = orig_mesh_pt->finite_element_pt(e)->node_pt(j);

        /// Get the (pseudo-)global node number in original mesh
        /// (It's zero [=default] if not visited this one yet)

        /// The (pseudo-)global node number in original mesh
        unsigned j_global = global_number[orig_node_pt];

        /// Haven't setup this node yet
        if (j_global == 0)
        {
          /// Give it a number (not necessarily the global node
          /// number in the original mesh -- we just need something
          /// to keep track...)
          global_count++;
          global_number[orig_node_pt] = global_count;

          /// Copy new node, created using the NEW element's construct_node
          /// function into global storage, using the same global
          /// node number that we've just associated with the
          /// corresponding node in the original mesh
          this->Node_pt[global_count - 1] =
            this->finite_element_pt(e)->node_pt(j);

          /// Assign coordinates
          this->Node_pt[global_count - 1]->x(0) = orig_node_pt->x(0);
          this->Node_pt[global_count - 1]->x(1) = orig_node_pt->x(1);


          /// Get pointer to set of mesh boundaries that this
          /// original node occupies; NULL if the node is not on any boundary
          std::set<unsigned>* boundaries_pt;
          orig_node_pt->get_boundaries_pt(boundaries_pt);

          /// Loop over the mesh boundaries that the node in the original mesh
          /// occupies and assign new node to the same ones.
          if (boundaries_pt != 0)
          {
            /// Convert it to a boundary node
            this->convert_to_boundary_node(this->Node_pt[global_count - 1]);
            /// Add the node to the boundaries
            for (std::set<unsigned>::iterator it = boundaries_pt->begin();
                 it != boundaries_pt->end();
                 ++it)
            {
              this->add_boundary_node(*it, this->Node_pt[global_count - 1]);
            }
          }
        }
        /// This one has already been setup/visited, therefore it is not
    unique,
        /// so we delete it
        else
        {
          /// Delete this node
          delete this->finite_element_pt(e)->node_pt(j);

          /// Overwrite the element's pointer to local node by
          /// pointer to the already existing node -- identified by
          /// the number kept in the map
          this->finite_element_pt(e)->node_pt(j) = this->Node_pt[j_global -
    1];
        }
      }
    }
*/

    // ============================================================

    // Wipe/allocate storage for arrays
    this->Boundary_element_pt.clear();
    this->Face_index_at_boundary.clear();
    this->Boundary_element_pt.resize(nbound);
    this->Face_index_at_boundary.resize(nbound);

    Vector<GeneralisedElement*> orig_element_pt = orig_mesh_pt->element_pt();

    // Loop over the boundaries
    for (unsigned b = 0; b < nbound; b++)
    {
      unsigned n_bound_ele = orig_mesh_pt->nboundary_element(b);
      for (unsigned n = 0; n < n_bound_ele; n++)
      {
        unsigned el_index =
          std::distance(orig_element_pt.begin(),
                        std::find(orig_element_pt.begin(),
                                  orig_element_pt.end(),
                                  orig_mesh_pt->boundary_element_pt(b, n)));
        this->Boundary_element_pt[b].push_back(
          dynamic_cast<FiniteElement*>(this->element_pt(el_index)));
        this->Face_index_at_boundary[b].push_back(
          orig_mesh_pt->face_index_at_boundary(b, n));
      }
    }
    this->Lookup_for_elements_next_boundary_is_setup = true;

    // ============================================================

    // this->setup_boundary_element_info();
  }
} // namespace oomph
#endif
