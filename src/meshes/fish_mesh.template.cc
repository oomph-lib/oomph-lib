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
#ifndef OOMPH_FISH_MESH_TEMPLATE_CC
#define OOMPH_FISH_MESH_TEMPLATE_CC

#include "fish_mesh.template.h"


namespace oomph
{
  //=================================================================
  /// Constructor: Pass pointer to timestepper.
  /// (defaults to (Steady) default timestepper defined in Mesh)
  //=================================================================
  template<class ELEMENT>
  FishMesh<ELEMENT>::FishMesh(TimeStepper* time_stepper_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Fish back is a circle of radius 1, centred at (0.5,0.0)
    double x_c = 0.5;
    double y_c = 0.0;
    double r_back = 1.0;
    Back_pt = new Circle(x_c, y_c, r_back);

    // I've created the fishback -- I need to kill it.
    Must_kill_fish_back = true;

    // Now build the mesh
    build_mesh(time_stepper_pt);
  }


  //=================================================================
  /// Constructor: Pass pointer GeomObject that defines
  /// the fish's back and pointer to timestepper.
  /// (defaults to (Steady) default timestepper defined in Mesh)
  //=================================================================
  template<class ELEMENT>
  FishMesh<ELEMENT>::FishMesh(GeomObject* back_pt, TimeStepper* time_stepper_pt)
    : Back_pt(back_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Back_pt has already been set....
    Must_kill_fish_back = false;

    // Now build the mesh
    build_mesh(time_stepper_pt);
  }


  //============================start_build_mesh=====================
  /// Build the mesh, using the geometric object that
  /// defines the fish's back.
  //=================================================================
  template<class ELEMENT>
  void FishMesh<ELEMENT>::build_mesh(TimeStepper* time_stepper_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Build domain: Pass the pointer to the geometric object that
    // represents the fish's back (pointed to by the FishMesh's
    // private data member, Back_pt, and the values of the
    // Lagrangian coordinate along this object at the "tail"
    // and the "nose":
    double xi_lo = 2.6;
    double xi_hi = 0.4;
    Domain_pt = new FishDomain(Back_pt, xi_lo, xi_hi);

    // Plot the domain? Keep this here in case we need to doc it
    bool plot_it = false;
    if (plot_it)
    {
      // Output the domain
      std::ofstream some_file;

      // Number of plot points in each coordinate direction.
      unsigned npts = 10;

      // Output domain
      some_file.open("fish_domain.dat");
      Domain_pt->output(some_file, npts);
      Domain_pt->output_macro_element_boundaries(some_file, npts);
      some_file.close();
    }

    // Set the number of boundaries
    set_nboundary(7);


    // We will store boundary coordinates on the curvilinear boundaries
    //(boundaries  0 and 4) along the fish's belly and its back.
    Boundary_coordinate_exists[0] = true;
    Boundary_coordinate_exists[4] = true;

    // Allocate the storage for the elements
    unsigned nelem = 4;
    Element_pt.resize(nelem);

    // Create  dummy element so we can determine the number of nodes
    ELEMENT* dummy_el_pt = new ELEMENT;

    // Read out the number of linear points in the element
    unsigned n_node_1d = dummy_el_pt->nnode_1d();

    // Kill the element
    delete dummy_el_pt;

    // Can now allocate the store for the nodes
    unsigned nnodes_total =
      (2 * (n_node_1d - 1) + 1) * (2 * (n_node_1d - 1) + 1);
    Node_pt.resize(nnodes_total);


    Vector<double> s(2), s_fraction(2);
    Vector<double> r(2);


    // Create elements and all nodes in element
    //-----------------------------------------
    // (ignore repetitions for now -- we'll clean them up later)
    //----------------------------------------------------------
    for (unsigned e = 0; e < nelem; e++)
    {
      // Create element
      Element_pt[e] = new ELEMENT;

      // Loop over rows in y/s_1-direction
      for (unsigned i1 = 0; i1 < n_node_1d; i1++)
      {
        // Loop over rows in x/s_0-direction
        for (unsigned i0 = 0; i0 < n_node_1d; i0++)
        {
          // Local node number
          unsigned j_local = i0 + i1 * n_node_1d;

          // Create the node and store pointer to it
          Node* node_pt =
            finite_element_pt(e)->construct_node(j_local, time_stepper_pt);

          // Work out the node's coordinates in the finite element's local
          // coordinate system:
          finite_element_pt(e)->local_fraction_of_node(j_local, s_fraction);

          s[0] = -1.0 + 2.0 * s_fraction[0];
          s[1] = -1.0 + 2.0 * s_fraction[1];

          // Get the global position of the node from macro element mapping
          Domain_pt->macro_element_pt(e)->macro_map(s, r);

          // Set the nodal position
          node_pt->x(0) = r[0];
          node_pt->x(1) = r[1];
        }
      }
    } // end of loop over elements


    // Kill repeated nodes and replace by pointers to nodes in appropriate
    //---------------------------------------------------------------------
    // neighbour. Also add node pointers to boundary arrays.
    //------------------------------------------------------

    // Initialise number of global nodes
    unsigned node_count = 0;

    // Check for error in node killing
    bool stopit = false;


    // Max tolerance for error in node killing
    double Max_tol_in_node_killing = 1.0e-12;

    // First element: Lower body: all nodes survive
    //---------------------------------------------
    unsigned e = 0;

    // Loop over rows in y/s_1-direction
    for (unsigned i1 = 0; i1 < n_node_1d; i1++)
    {
      // Loop over rows in x/s_0-direction
      for (unsigned i0 = 0; i0 < n_node_1d; i0++)
      {
        // Local node number
        unsigned j_local = i0 + i1 * n_node_1d;

        // No duplicate node: Copy new node across to mesh
        Node_pt[node_count] = finite_element_pt(e)->node_pt(j_local);

        // Set boundaries:

        // If we're on the boundary need to convert the node
        // into a boundary node
        if ((i0 == 0) || (i1 == 0))
        {
          this->convert_to_boundary_node(Node_pt[node_count]);
        }

        // Lower jaw: boundary 6
        if (i0 == 0)
        {
          add_boundary_node(6, Node_pt[node_count]);
        }

        // Belly: boundary 0
        if (i1 == 0)
        {
          add_boundary_node(0, Node_pt[node_count]);

          // Set the boundary coordinate
          Vector<double> zeta(1);
          zeta[0] =
            xi_lo + (xi_hi - xi_lo) * double(i0) / double(n_node_1d - 1);
          Node_pt[node_count]->set_coordinates_on_boundary(0, zeta);
        }

        // Increment node counter
        node_count++;
      }
    }


    // Lower fin: Western row of nodes is duplicate from element 0
    //------------------------------------------------------------
    e = 1;

    // Loop over rows in y/s_1-direction
    for (unsigned i1 = 0; i1 < n_node_1d; i1++)
    {
      // Loop over rows in x/s_0-direction
      for (unsigned i0 = 0; i0 < n_node_1d; i0++)
      {
        // Local node number
        unsigned j_local = i0 + i1 * n_node_1d;

        // Has the node been killed?
        bool killed = false;

        // First vertical row of nodes in s_1 direction get killed
        // and re-directed to nodes in element 0
        if (i0 == 0)
        {
          // Neighbour element
          unsigned e_neigh = 0;

          // Node in neighbour element
          unsigned i0_neigh = n_node_1d - 1;
          unsigned i1_neigh = i1;
          unsigned j_local_neigh = i0_neigh + i1_neigh * n_node_1d;


          // Check:
          for (unsigned i = 0; i < 2; i++)
          {
            double error = std::fabs(
              finite_element_pt(e)->node_pt(j_local)->x(i) -
              finite_element_pt(e_neigh)->node_pt(j_local_neigh)->x(i));
            if (error > Max_tol_in_node_killing)
            {
              oomph_info << "Error in node killing for i " << i << " " << error
                         << std::endl;
              stopit = true;
            }
          }

          // Kill node
          delete finite_element_pt(e)->node_pt(j_local);
          killed = true;

          // Set pointer to neighbour:
          finite_element_pt(e)->node_pt(j_local) =
            finite_element_pt(e_neigh)->node_pt(j_local_neigh);
        }


        // No duplicate node: Copy across to mesh
        if (!killed)
        {
          // Copy the node across
          Node_pt[node_count] = finite_element_pt(e)->node_pt(j_local);

          // If we're on a boundary turn the node into
          // a boundary node
          if ((i1 == 0) || (i0 == n_node_1d - 1))
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
          }

          // Increment node counter
          node_count++;
        }

        // Set boundaries:

        // Bottom of tail: boundary 1
        if (i1 == 0)
        {
          add_boundary_node(1, finite_element_pt(e)->node_pt(j_local));
        }

        // Vertical bit of tail: boundary 2
        if (i0 == n_node_1d - 1)
        {
          add_boundary_node(2, finite_element_pt(e)->node_pt(j_local));
        }
      }
    }


    // Upper body: Southern row of nodes is duplicate from element 0
    //--------------------------------------------------------------
    e = 2;

    // Loop over rows in y/s_1-direction
    for (unsigned i1 = 0; i1 < n_node_1d; i1++)
    {
      // Loop over rows in x/s_0-direction
      for (unsigned i0 = 0; i0 < n_node_1d; i0++)
      {
        // Local node number
        unsigned j_local = i0 + i1 * n_node_1d;

        // Has the node been killed?
        bool killed = false;

        // First horizontal row of nodes in s_0 direction get killed
        // and re-directed to nodes in element 0
        if (i1 == 0)
        {
          // Neighbour element
          unsigned e_neigh = 0;

          // Node in neighbour element
          unsigned i0_neigh = i0;
          unsigned i1_neigh = n_node_1d - 1;
          unsigned j_local_neigh = i0_neigh + i1_neigh * n_node_1d;

          // Check:
          for (unsigned i = 0; i < 2; i++)
          {
            double error = std::fabs(
              finite_element_pt(e)->node_pt(j_local)->x(i) -
              finite_element_pt(e_neigh)->node_pt(j_local_neigh)->x(i));
            if (error > Max_tol_in_node_killing)
            {
              oomph_info << "Error in node killing for i " << i << " " << error
                         << std::endl;
              stopit = true;
            }
          }

          // Kill node
          delete finite_element_pt(e)->node_pt(j_local);
          killed = true;

          // Set pointer to neighbour:
          finite_element_pt(e)->node_pt(j_local) =
            finite_element_pt(e_neigh)->node_pt(j_local_neigh);
        }

        // No duplicate node: Copy across to mesh
        if (!killed)
        {
          // Copy the old node across to the mesh
          Node_pt[node_count] = finite_element_pt(e)->node_pt(j_local);

          // If we're on a boundary, convert the node into a boundary
          // node. This will automatically update the entry in the mesh
          if ((i1 == n_node_1d - 1) || (i0 == 0))
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
          }

          // Increment node counter
          node_count++;
        }

        // Set boundaries:

        // Back: boundary 4
        if (i1 == n_node_1d - 1)
        {
          add_boundary_node(4, finite_element_pt(e)->node_pt(j_local));

          // Set the boundary coordinate
          Vector<double> zeta(1);
          zeta[0] =
            xi_lo + (xi_hi - xi_lo) * double(i0) / double(n_node_1d - 1);
          finite_element_pt(e)->node_pt(j_local)->set_coordinates_on_boundary(
            4, zeta);
        }

        // Upper jaw: boundary 5
        if (i0 == 0)
        {
          add_boundary_node(5, finite_element_pt(e)->node_pt(j_local));
        }
      }
    }


    // Upper fin: Western/southern row of nodes is duplicate from element 2/1
    //-----------------------------------------------------------------------
    e = 3;

    // Loop over rows in y/s_1-direction
    for (unsigned i1 = 0; i1 < n_node_1d; i1++)
    {
      // Loop over rows in x/s_0-direction
      for (unsigned i0 = 0; i0 < n_node_1d; i0++)
      {
        // Local node number
        unsigned j_local = i0 + i1 * n_node_1d;

        // Has the node been killed?
        bool killed = false;

        // First vertical row of nodes in s_1 direction get killed
        // and re-directed to nodes in element 2
        if (i0 == 0)
        {
          // Neighbour element
          unsigned e_neigh = 2;

          // Node in neighbour element
          unsigned i0_neigh = n_node_1d - 1;
          unsigned i1_neigh = i1;
          unsigned j_local_neigh = i0_neigh + i1_neigh * n_node_1d;


          // Check:
          for (unsigned i = 0; i < 2; i++)
          {
            double error = std::fabs(
              finite_element_pt(e)->node_pt(j_local)->x(i) -
              finite_element_pt(e_neigh)->node_pt(j_local_neigh)->x(i));
            if (error > Max_tol_in_node_killing)
            {
              oomph_info << "Error in node killing for i " << i << " " << error
                         << std::endl;
              stopit = true;
            }
          }

          // Kill node
          delete finite_element_pt(e)->node_pt(j_local);
          killed = true;

          // Set pointer to neighbour:
          finite_element_pt(e)->node_pt(j_local) =
            finite_element_pt(e_neigh)->node_pt(j_local_neigh);
        }


        // First horizontal row of nodes in s_0 direction (apart from
        // first node get killed and re-directed to nodes in element 1
        if ((i0 != 0) && (i1 == 0))
        {
          // Neighbour element
          unsigned e_neigh = 1;

          // Node in neighbour element
          unsigned i0_neigh = i0;
          unsigned i1_neigh = n_node_1d - 1;
          unsigned j_local_neigh = i0_neigh + i1_neigh * n_node_1d;


          // Check:
          for (unsigned i = 0; i < 2; i++)
          {
            double error = std::fabs(
              finite_element_pt(e)->node_pt(j_local)->x(i) -
              finite_element_pt(e_neigh)->node_pt(j_local_neigh)->x(i));
            if (error > Max_tol_in_node_killing)
            {
              oomph_info << "Error in node killing for i " << i << " " << error
                         << std::endl;
              stopit = true;
            }
          }

          // Kill node
          delete finite_element_pt(e)->node_pt(j_local);
          killed = true;

          // Set pointer to neighbour:
          finite_element_pt(e)->node_pt(j_local) =
            finite_element_pt(e_neigh)->node_pt(j_local_neigh);
        }


        // No duplicate node: Copy across to mesh
        if (!killed)
        {
          // Now copy the node across
          Node_pt[node_count] = finite_element_pt(e)->node_pt(j_local);

          // If we're on the boundary, convert the node into
          // a boundary node
          if ((i1 == n_node_1d - 1) || (i0 == n_node_1d - 1))
          {
            this->convert_to_boundary_node(Node_pt[node_count]);
          }

          // Increment node counter
          node_count++;
        }

        // Set boundaries:

        // To of tail: boundary 3
        if (i1 == n_node_1d - 1)
        {
          add_boundary_node(3, finite_element_pt(e)->node_pt(j_local));
        }


        // Vertical bit of tail: boundary 2
        if (i0 == n_node_1d - 1)
        {
          add_boundary_node(2, finite_element_pt(e)->node_pt(j_local));
        }
      }
    }

    // Terminate if there's been an error
    if (stopit)
    {
      std::ostringstream error_message;
      error_message << "Error occured in node killing!\n";
      error_message
        << "Max. permitted difference in position of the two nodes\n "
        << "that get 'merged' : " << Max_tol_in_node_killing << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


    // Loop over all elements and set macro element pointer
    unsigned n_element = this->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      // Get pointer to element
      FiniteElement* el_pt = this->finite_element_pt(e);

      // Set pointer to macro element to enable MacroElement-based
      // remesh. Also enables the curvlinear boundaries
      // of the mesh/domain get picked up during adaptive
      // mesh refinement in derived classes.
      el_pt->set_macro_elem_pt(this->Domain_pt->macro_element_pt(e));
    }


    // Setup boundary element lookup schemes
    setup_boundary_element_info();


    // Check the boundary coordinates
#ifdef PARANOID
    {
      Vector<double> zeta(1);
      Vector<double> r(2);
      bool stopit = false;
      unsigned num_bound = nboundary();
      for (unsigned ibound = 0; ibound < num_bound; ibound++)
      {
        if (boundary_coordinate_exists(ibound))
        {
          unsigned num_nod = nboundary_node(ibound);
          for (unsigned inod = 0; inod < num_nod; inod++)
          {
            // Get the boundary coordinate
            boundary_node_pt(ibound, inod)
              ->get_coordinates_on_boundary(ibound, zeta);

            // Get position from wall object
            Back_pt->position(zeta, r);

            // Flip it
            if (ibound == 0) r[1] = -r[1];

            // Check:
            for (unsigned i = 0; i < 2; i++)
            {
              double error =
                std::fabs(r[i] - boundary_node_pt(ibound, inod)->x(i));
              if (error > Max_tol_in_node_killing)
              {
                oomph_info << "Error in boundary coordinate for direction " << i
                           << " on boundary " << ibound << ":" << error
                           << std::endl;

                oomph_info << "x: " << r[0] << " "
                           << boundary_node_pt(ibound, inod)->x(0) << std::endl;

                oomph_info << "y: " << r[1] << " "
                           << boundary_node_pt(ibound, inod)->x(1) << std::endl
                           << std::endl;
                stopit = true;
              }
            }
          }
        }
      }

      // Terminate if there's been an error
      if (stopit)
      {
        std::ostringstream error_message;
        error_message << "Error occured in boundary coordinate setup!\n";
        error_message << "Max. tolerance: " << Max_tol_in_node_killing
                      << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif
  }


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //=========start_setup_adaptivity=========================================
  ///  Setup all the information that's required for spatial adaptivity:
  /// Build quadtree forest.
  //========================================================================
  template<class ELEMENT>
  void RefineableFishMesh<ELEMENT>::setup_adaptivity()
  {
    // Setup quadtree forest
    this->setup_quadtree_forest();

  } // end of setup_adaptivity


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // AlgebraicElement fish-shaped mesh
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //======================================================================
  /// \short Setup algebraic update operation. Nodes are "suspended"
  /// from the fish's back and the upper edge of the fin. Nodes
  /// in the lower half are placed symmetrically.
  //======================================================================
  template<class ELEMENT>
  void AlgebraicFishMesh<ELEMENT>::setup_algebraic_node_update()
  {
#ifdef PARANOID
    /// Pointer to algebraic element in lower body
    AlgebraicElementBase* lower_body_pt =
      dynamic_cast<AlgebraicElementBase*>(Mesh::element_pt(0));

    if (lower_body_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "Element in AlgebraicFishMesh must be\n"
                    << "derived from AlgebraicElementBase\n"
                    << "but it is of type: "
                    << typeid(Mesh::element_pt(0)).name() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Read out the number of linear points in the element
    unsigned n_p =
      dynamic_cast<ELEMENT*>(FishMesh<ELEMENT>::Mesh::finite_element_pt(0))
        ->nnode_1d();

    // Element 0: Lower body
    //----------------------
    {
      unsigned ielem = 0;
      FiniteElement* el_pt = Mesh::finite_element_pt(ielem);

      // Loop over rows in y/s_1-direction
      for (unsigned i1 = 0; i1 < n_p; i1++)
      {
        // Loop over rows in x/s_0-direction
        for (unsigned i0 = 0; i0 < n_p; i0++)
        {
          // Local node number
          unsigned jnod = i0 + i1 * n_p;

          // One geometric object is involved in update operation
          Vector<GeomObject*> geom_object_pt(1);
          geom_object_pt[0] = this->Back_pt;

          // The update function requires three parameters:
          Vector<double> ref_value(3);

          // First reference value: fractional x-position
          ref_value[0] = double(i0) / double(n_p - 1);

          // Second reference value: fractional position along
          // straight line from position on horizontal symmetry line to
          // point on fish back
          ref_value[1] = 1.0 - double(i1) / double(n_p - 1);

          // Third reference value: Sign (are we above or below the
          // symmetry line?)
          ref_value[2] = -1.0;

          // Setup algebraic update for node: Pass update information
          dynamic_cast<AlgebraicNode*>(el_pt->node_pt(jnod))
            ->add_node_update_info(this->Lower_body, // enumerated ID
                                   this, // mesh
                                   geom_object_pt, // vector of geom objects
                                   ref_value); // vector of ref. values
        }
      }
    }


    // Element 1: Lower fin
    //---------------------
    {
      unsigned ielem = 1;
      FiniteElement* el_pt = Mesh::finite_element_pt(ielem);

      // Loop over rows in y/s_1-direction
      for (unsigned i1 = 0; i1 < n_p; i1++)
      {
        // Loop over rows in x/s_0-direction
        for (unsigned i0 = 0; i0 < n_p; i0++)
        {
          // Local node number
          unsigned jnod = i0 + i1 * n_p;

          // One geometric object is involved in update operation
          Vector<GeomObject*> geom_object_pt(1);
          geom_object_pt[0] = this->Back_pt;

          // The update function requires three parameters:
          Vector<double> ref_value(3);

          // First reference value: fractional x-position
          ref_value[0] = double(i0) / double(n_p - 1);

          // Second reference value: fractional position along
          // straight line from position on horizontal symmetry line to
          // point on fish back
          ref_value[1] = 1.0 - double(i1) / double(n_p - 1);

          // Third reference value: Sign (are we above or below the
          // symmetry line?)
          ref_value[2] = -1.0;

          // Setup algebraic update for node: Pass update information
          dynamic_cast<AlgebraicNode*>(el_pt->node_pt(jnod))
            ->add_node_update_info(this->Lower_fin, // enumerated ID
                                   this, // mesh
                                   geom_object_pt, // vector of geom objects
                                   ref_value); // vector of ref. values
        }
      }
    }


    // Element 2: Upper body
    //----------------------
    {
      unsigned ielem = 2;
      FiniteElement* el_pt = Mesh::finite_element_pt(ielem);

      // Loop over rows in y/s_1-direction
      for (unsigned i1 = 0; i1 < n_p; i1++)
      {
        // Loop over rows in x/s_0-direction
        for (unsigned i0 = 0; i0 < n_p; i0++)
        {
          // Local node number
          unsigned jnod = i0 + i1 * n_p;

          // One geometric object is involved in update operation
          Vector<GeomObject*> geom_object_pt(1);
          geom_object_pt[0] = this->Back_pt;

          // The update function requires three parameters:
          Vector<double> ref_value(3);

          // First reference value: fractional x-position
          ref_value[0] = double(i0) / double(n_p - 1);

          // Second reference value: fractional position along
          // straight line from position on horizontal symmetry line to
          // point on fish back
          ref_value[1] = double(i1) / double(n_p - 1);

          // Third reference value: Sign (are we above or below the
          // symmetry line?)
          ref_value[2] = 1.0;

          // Setup algebraic update for node: Pass update information
          dynamic_cast<AlgebraicNode*>(el_pt->node_pt(jnod))
            ->add_node_update_info(this->Upper_body, // enumerated ID
                                   this, // mesh
                                   geom_object_pt, // vector of geom objects
                                   ref_value); // vector of ref. values
        }
      }
    }


    // Element 3: Upper fin
    //---------------------
    {
      unsigned ielem = 3;
      FiniteElement* el_pt = Mesh::finite_element_pt(ielem);

      // Loop over rows in y/s_1-direction
      for (unsigned i1 = 0; i1 < n_p; i1++)
      {
        // Loop over rows in x/s_0-direction
        for (unsigned i0 = 0; i0 < n_p; i0++)
        {
          // Local node number
          unsigned jnod = i0 + i1 * n_p;

          // One geometric object is involved in update operation
          Vector<GeomObject*> geom_object_pt(1);
          geom_object_pt[0] = this->Back_pt;

          // The update function requires three parameters:
          Vector<double> ref_value(3);

          // First reference value: fractional x-position
          ref_value[0] = double(i0) / double(n_p - 1);

          // Second reference value: fractional position along
          // straight line from position on horizontal symmetry line to
          // point on fish back
          ref_value[1] = double(i1) / double(n_p - 1);

          // Third reference value: Sign (are we above or below the
          // symmetry line?)
          ref_value[2] = 1.0;

          // Setup algebraic update for node: Pass update information
          dynamic_cast<AlgebraicNode*>(el_pt->node_pt(jnod))
            ->add_node_update_info(this->Upper_fin, // enumerated ID
                                   this, // mesh
                                   geom_object_pt, // vector of geom objects
                                   ref_value); // vector of ref. values
        }
      }
    }
  }


  //======================================================================
  /// \short Algebraic update function: Update in (upper or lower) body
  /// according to wall shape at time level t (t=0: present; t>0: previous)
  //======================================================================
  template<class ELEMENT>
  void AlgebraicFishMesh<ELEMENT>::node_update_in_body(const unsigned& t,
                                                       AlgebraicNode*& node_pt)
  {
    // Pointer to geometric object that represents the fish back:
    GeomObject* back_pt = node_pt->geom_object_pt(unsigned(0));

    // Fixed reference value: x-position of mouth
    double x_mouth = this->Domain_pt->x_mouth();

    // Fixed reference value: Lagrangian coordinate of point
    // over mouth
    double zeta_mouth = this->Domain_pt->xi_nose();

    // Fixed reference value: Lagrangian coordinate of point
    // near tail
    double zeta_near_tail = this->Domain_pt->xi_tail();

    // First reference value: fractional x-position
    double fract_x = node_pt->ref_value(unsigned(0));

    // Second reference value: fractional position along
    // straight line from position on horizontal symmetry line to
    // point on fish back
    double fract_y = node_pt->ref_value(1);

    // Third reference value: Sign (are we above or below the
    // symmetry line?)
    double sign = node_pt->ref_value(2);

    // Get position on fish back
    Vector<double> zeta(back_pt->nlagrangian());
    zeta[0] = zeta_mouth + fract_x * (zeta_near_tail - zeta_mouth);
    Vector<double> r_back(back_pt->ndim());
    back_pt->position(t, zeta, r_back);

    // Get position of point on fish back near tail
    zeta[0] = zeta_near_tail;
    Vector<double> r_near_tail(back_pt->ndim());
    back_pt->position(t, zeta, r_near_tail);

    // Get position on symmetry line
    Vector<double> r_sym(2);
    r_sym[0] = x_mouth + fract_x * (r_near_tail[0] - x_mouth);
    r_sym[1] = 0.0;

    // Assign new nodal coordinate
    node_pt->x(t, 0) = r_sym[0] + fract_y * (r_back[0] - r_sym[0]);
    node_pt->x(t, 1) = sign * (r_sym[1] + fract_y * (r_back[1] - r_sym[1]));
  }


  //======================================================================
  /// \short Algebraic update function: Update in (upper or lower) fin
  /// according to wall shape at time level t (t=0: present; t>0: previous)
  //======================================================================
  template<class ELEMENT>
  void AlgebraicFishMesh<ELEMENT>::node_update_in_fin(const unsigned& t,
                                                      AlgebraicNode*& node_pt)
  {
    // Pointer to geometric object that represents the fish back:
    GeomObject* back_pt = node_pt->geom_object_pt(unsigned(0));

    // Fixed reference value: End coordinate on fish back
    double zeta_wall = this->Domain_pt->xi_tail();

    // Fixed reference value: x-position of end of tail
    double x_tail = this->Domain_pt->x_fin();

    // Fixed reference value: y-position of end of tail
    double y_tail = this->Domain_pt->y_fin();

    // First reference value: fractional position in x-direction
    double fract_x = node_pt->ref_value(unsigned(0));

    // Second reference value: fractional position along
    // vertical line from position on horizontal symmetry line to
    // point on upper end of tail
    double fract_y = node_pt->ref_value(1);

    // Third reference value: Sign (are we above or below the
    // symmetry line?)
    double sign = node_pt->ref_value(2);

    // Get position on fish back
    Vector<double> zeta(back_pt->nlagrangian());
    zeta[0] = zeta_wall;
    Vector<double> r_back(back_pt->ndim());
    back_pt->position(t, zeta, r_back);

    // y-position on top edge of fin:
    double y_fin_edge = r_back[1] + fract_x * (y_tail - r_back[1]);

    // Assign new nodal coordinate
    node_pt->x(t, 0) = r_back[0] + fract_x * (x_tail - r_back[0]);
    node_pt->x(t, 1) = sign * fract_y * y_fin_edge;
  }

} // namespace oomph
#endif
