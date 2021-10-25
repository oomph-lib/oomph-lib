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
#ifndef OOMPH_TWO_LAYER_SPINE_MESH_TEMPLATE_CC
#define OOMPH_TWO_LAYER_SPINE_MESH_TEMPLATE_CC

#include "two_layer_spine_mesh.template.h"
#include "rectangular_quadmesh.template.cc"


namespace oomph
{
  //===========================================================================
  /// Constuctor for spine 2D mesh: Pass number of elements in x-direction,
  /// number of elements in y-direction in bottom and top layer, respectively,
  /// axial length and height of top and bottom layers, and pointer to
  /// timestepper (defaults to Static timestepper).
  ///
  /// The mesh contains two layers of elements (of type ELEMENT;
  /// e.g  SpineElement<QCrouzeixRaviartElement<2>)
  /// and an interfacial layer of corresponding Spine interface elements
  /// of type INTERFACE_ELEMENT, e.g.
  /// SpineLineFluidInterfaceElement<ELEMENT> for 2D planar
  /// problems.
  //===========================================================================
  template<class ELEMENT>
  TwoLayerSpineMesh<ELEMENT>::TwoLayerSpineMesh(const unsigned& nx,
                                                const unsigned& ny1,
                                                const unsigned& ny2,
                                                const double& lx,
                                                const double& h1,
                                                const double& h2,
                                                TimeStepper* time_stepper_pt)
    : RectangularQuadMesh<ELEMENT>(
        nx, ny1 + ny2, 0.0, lx, 0.0, h1 + h2, false, false, time_stepper_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Mesh can only be built with spine elements
    MeshChecker::assert_geometric_element<SpineFiniteElement, ELEMENT>(2);

    // We've called the "generic" constructor for the RectangularQuadMesh
    // which doesn't do much...
    // Now set up the parameters that characterise the mesh geometry
    // then call the build function

    // Number of elements in bottom and top layers
    Ny1 = ny1;
    Ny2 = ny2;

    // Set height of upper and lower layers
    H1 = h1;
    H2 = h2;

    // Now build the mesh:
    build_two_layer_mesh(time_stepper_pt);
  }


  //===========================================================================
  /// Constuctor for spine 2D mesh: Pass number of elements in x-direction,
  /// number of elements in y-direction in bottom and top layer, respectively,
  /// axial length and height of top and bottom layers, a boolean
  /// flag to make the mesh periodic in the x-direction, and pointer to
  /// timestepper (defaults to Static timestepper).
  ///
  /// The mesh contains two layers of elements (of type ELEMENT;
  /// e.g  SpineElement<QCrouzeixRaviartElement<2>)
  /// and an interfacial layer of corresponding Spine interface elements
  /// of type INTERFACE_ELEMENT, e.g.
  /// SpineLineFluidInterfaceElement<ELEMENT> for 2D planar
  /// problems.
  //===========================================================================
  template<class ELEMENT>
  TwoLayerSpineMesh<ELEMENT>::TwoLayerSpineMesh(const unsigned& nx,
                                                const unsigned& ny1,
                                                const unsigned& ny2,
                                                const double& lx,
                                                const double& h1,
                                                const double& h2,
                                                const bool& periodic_in_x,
                                                TimeStepper* time_stepper_pt)
    : RectangularQuadMesh<ELEMENT>(nx,
                                   ny1 + ny2,
                                   0.0,
                                   lx,
                                   0.0,
                                   h1 + h2,
                                   periodic_in_x,
                                   false,
                                   time_stepper_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Mesh can only be built with spine elements
    MeshChecker::assert_geometric_element<SpineFiniteElement, ELEMENT>(2);

    // We've called the "generic" constructor for the RectangularQuadMesh
    // which doesn't do much...
    // Now set up the parameters that characterise the mesh geometry
    // then call the constructor

    // Number of elements in bottom and top layers
    Ny1 = ny1;
    Ny2 = ny2;

    // Set height of upper and lower layers
    H1 = h1;
    H2 = h2;

    // Now build the mesh:
    build_two_layer_mesh(time_stepper_pt);
  }


  //===========================================================================
  /// Constuctor for spine 2D mesh: Pass number of elements in x-direction,
  /// number of elements in y-direction in bottom and top layer, respectively,
  /// axial length and height of top and bottom layers, a boolean
  /// flag to make the mesh periodic in the x-direction, a boolean flag to
  /// specify whether or not to call the "build_two_layer_mesh" function,
  /// and pointer to timestepper (defaults to Static timestepper).
  ///
  /// The mesh contains two layers of elements (of type ELEMENT;
  /// e.g  SpineElement<QCrouzeixRaviartElement<2>)
  /// and an interfacial layer of corresponding Spine interface elements
  /// of type INTERFACE_ELEMENT, e.g.
  /// SpineLineFluidInterfaceElement<ELEMENT> for 2D planar
  /// problems.
  //===========================================================================
  template<class ELEMENT>
  TwoLayerSpineMesh<ELEMENT>::TwoLayerSpineMesh(const unsigned& nx,
                                                const unsigned& ny1,
                                                const unsigned& ny2,
                                                const double& lx,
                                                const double& h1,
                                                const double& h2,
                                                const bool& periodic_in_x,
                                                const bool& build_mesh,
                                                TimeStepper* time_stepper_pt)
    : RectangularQuadMesh<ELEMENT>(nx,
                                   ny1 + ny2,
                                   0.0,
                                   lx,
                                   0.0,
                                   h1 + h2,
                                   periodic_in_x,
                                   false,
                                   time_stepper_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Mesh can only be built with spine elements
    MeshChecker::assert_geometric_element<SpineFiniteElement, ELEMENT>(2);

    // We've called the "generic" constructor for the RectangularQuadMesh
    // which doesn't do much...
    // Now set up the parameters that characterise the mesh geometry
    // then call the constructor

    // Number of elements in bottom and top layers
    Ny1 = ny1;
    Ny2 = ny2;

    // Set height of upper and lower layers
    H1 = h1;
    H2 = h2;

    // Only build the mesh here if build_mesh=true
    // This is useful when calling this constructor from a derived class
    // (such as Axisym2x6TwoLayerSpineMesh) where the mesh building
    // needs to be called from *its* constructor and this constructor is
    // only used to pass arguments to the RectangularQuadMesh constructor.
    if (build_mesh)
    {
      build_two_layer_mesh(time_stepper_pt);
    }
  }


  //==================================================================
  /// The spacing function for the x co-ordinate, which is the
  /// same as the default function.
  //==================================================================
  template<class ELEMENT>
  double TwoLayerSpineMesh<ELEMENT>::x_spacing_function(unsigned xelement,
                                                        unsigned xnode,
                                                        unsigned yelement,
                                                        unsigned ynode)
  {
    // Calculate the values of equal increments in nodal values
    double xstep = (this->Xmax - this->Xmin) / ((this->Np - 1) * this->Nx);
    // Return the appropriate value
    return (this->Xmin + xstep * ((this->Np - 1) * xelement + xnode));
  }

  //==================================================================
  /// The spacing function for the y co-ordinates, which splits
  /// the region into two regions (1 and 2), according to the
  /// heights H1 and H2, with Ny1 and Ny2 elements respectively.
  template<class ELEMENT>
  double TwoLayerSpineMesh<ELEMENT>::y_spacing_function(unsigned xelement,
                                                        unsigned xnode,
                                                        unsigned yelement,
                                                        unsigned ynode)
  {
    // Set up some spacing parameters
    // The lower region a starts at Ymin
    double Ymin = RectangularQuadMesh<ELEMENT>::Ymin;
    // The upper region a ends at Ymax
    double Ymax = RectangularQuadMesh<ELEMENT>::Ymax;
    // Number of nodes per element
    unsigned n_p = RectangularQuadMesh<ELEMENT>::Np;
    // The lower region starts at Ymin
    double y1init = Ymin;
    // The upper region starts at H1 - Ymin
    double y2init = H1 - Ymin;
    // Calculate the space between each node in each region,
    // Assumming uniform spacing
    // Lower region has a length (H1-Ymin)
    double y1step = (H1 - Ymin) / ((n_p - 1) * Ny1);
    // Upper region has a length (Ymax-H1)
    double y2step = (Ymax - H1) / ((n_p - 1) * Ny2);

    // Now return the actual node position, it's different in the two
    // regions, of course
    if (yelement < Ny1)
    {
      return (y1init + y1step * ((n_p - 1) * yelement + ynode));
    }
    else
    {
      return (y2init + y2step * ((n_p - 1) * (yelement - Ny1) + ynode));
    }
  }

  //===========================================================================
  /// Helper function that actually builds the two-layer spine mesh
  /// based on the parameters set in the various constructors
  //===========================================================================
  template<class ELEMENT>
  void TwoLayerSpineMesh<ELEMENT>::build_two_layer_mesh(
    TimeStepper* time_stepper_pt)
  {
    // Build the underlying quad mesh:
    RectangularQuadMesh<ELEMENT>::build_mesh(time_stepper_pt);

    // Reset the number of boundaries
    set_nboundary(7);

    // Set up the pointers to elements in the upper and lower fluid
    // Calculate numbers of nodes in upper and lower regions
    unsigned long n_lower = this->Nx * Ny1;
    unsigned long n_upper = this->Nx * Ny2;
    // Loop over lower elements and push back
    Lower_layer_element_pt.reserve(n_lower);
    for (unsigned e = 0; e < n_lower; e++)
    {
      Lower_layer_element_pt.push_back(this->finite_element_pt(e));
    }
    // Loop over upper elements and push back
    Upper_layer_element_pt.reserve(n_upper);
    for (unsigned e = n_lower; e < (n_lower + n_upper); e++)
    {
      Upper_layer_element_pt.push_back(this->finite_element_pt(e));
    }

    // Set the elements adjacent to the interface
    Interface_lower_boundary_element_pt.resize(this->Nx);
    Interface_upper_boundary_element_pt.resize(this->Nx);
    {
      unsigned count_lower = this->Nx * (this->Ny1 - 1);
      unsigned count_upper = count_lower + this->Nx;
      for (unsigned e = 0; e < this->Nx; e++)
      {
        Interface_lower_boundary_element_pt[e] =
          this->finite_element_pt(count_lower);
        ++count_lower;
        Interface_upper_boundary_element_pt[e] =
          this->finite_element_pt(count_upper);
        ++count_upper;
      }
    }

    // Relabel the boundary nodes
    // Storage for the boundary coordinates that will be transferred directly
    Vector<double> b_coord;
    {
      // Store the interface level
      const double y_interface = H1 - this->Ymin;

      // Nodes on boundary 3 are now on boundaries 4 and 5
      unsigned n_boundary_node = this->nboundary_node(3);
      // Loop over the nodes remove them from boundary 3
      // and add them to boundary 4 or 5 depending on their y coordinate
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        // Cache pointer to the node
        Node* const nod_pt = this->boundary_node_pt(3, n);
        // Get the boundary coordinates if set
        if (boundary_coordinate_exists(3))
        {
          b_coord.resize(nod_pt->ncoordinates_on_boundary(3));
          nod_pt->get_coordinates_on_boundary(3, b_coord);
          // Indicate that the boundary coordinates are to be set on the
          // new boundaries
          this->Boundary_coordinate_exists[4] = true;
          this->Boundary_coordinate_exists[5] = true;
        }

        // Now remove the node from the boundary
        nod_pt->remove_from_boundary(3);

        // Find the height of the node
        double y = nod_pt->x(1);
        // If it's above or equal to the interface, it's on boundary 4
        if (y >= y_interface)
        {
          this->add_boundary_node(4, nod_pt);
          // Add the boundary coordinate if it has been set up
          if (this->Boundary_coordinate_exists[4])
          {
            nod_pt->set_coordinates_on_boundary(4, b_coord);
          }
        }
        // otherwise it's on boundary 5
        if (y <= y_interface)
        {
          this->add_boundary_node(5, nod_pt);
          // Add the boundary coordinate if it has been set up
          if (this->Boundary_coordinate_exists[5])
          {
            nod_pt->set_coordinates_on_boundary(5, b_coord);
          }
        }
      }

      // Now clear the boundary node information from boundary 3
      this->Boundary_node_pt[3].clear();

      // Relabel the nodes on boundary 2 to be on boundary 3
      n_boundary_node = this->nboundary_node(2);
      // Loop over the nodes remove them from boundary 2
      // and add them to boundary 3
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        // Cache pointer to the node
        Node* const nod_pt = this->boundary_node_pt(2, n);
        // Get the boundary coordinates if set
        if (this->boundary_coordinate_exists(2))
        {
          b_coord.resize(nod_pt->ncoordinates_on_boundary(2));
          nod_pt->get_coordinates_on_boundary(2, b_coord);
          this->Boundary_coordinate_exists[3] = true;
        }

        // Now remove the node from the boundary 2
        nod_pt->remove_from_boundary(2);
        // and add to boundary 3
        this->add_boundary_node(3, nod_pt);
        if (this->Boundary_coordinate_exists[3])
        {
          nod_pt->set_coordinates_on_boundary(3, b_coord);
        }
      }

      // Clear the information from boundary 2
      this->Boundary_node_pt[2].clear();

      std::list<Node*> nodes_to_be_removed;

      // Nodes on boundary 1 are now on boundaries 1 and 2
      n_boundary_node = this->nboundary_node(1);
      // Loop over the nodes remove some of them from boundary 1
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        // Cache pointer to the node
        Node* const nod_pt = this->boundary_node_pt(1, n);

        // Find the height of the node
        double y = nod_pt->x(1);
        // If it's above or equal to the interface it's on boundary 2
        if (y >= y_interface)
        {
          // Get the boundary coordinates if set
          if (this->boundary_coordinate_exists(1))
          {
            b_coord.resize(nod_pt->ncoordinates_on_boundary(1));
            nod_pt->get_coordinates_on_boundary(1, b_coord);
            this->Boundary_coordinate_exists[2] = true;
          }

          // Now remove the node from the boundary 1 if above interace
          if (y > y_interface)
          {
            nodes_to_be_removed.push_back(nod_pt);
          }
          // Always add to boundary 2
          this->add_boundary_node(2, nod_pt);
          // Add the boundary coordinate if it has been set up
          if (this->Boundary_coordinate_exists[2])
          {
            nod_pt->set_coordinates_on_boundary(2, b_coord);
          }
        }
      }

      // Loop over the nodes that are to be removed from boundary 1 and remove
      // them
      for (std::list<Node*>::iterator it = nodes_to_be_removed.begin();
           it != nodes_to_be_removed.end();
           ++it)
      {
        this->remove_boundary_node(1, *it);
      }
      nodes_to_be_removed.clear();
    }
    // Allocate memory for the spines and fractions along spines
    //---------------------------------------------------------

    // Read out number of linear points in the element
    unsigned n_p = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();


    // Add the nodes on the interface to the boundary 6
    // Storage for boundary coordinates (x-coordinate)
    b_coord.resize(1);
    this->Boundary_coordinate_exists[6];
    // Starting index of the nodes
    unsigned n_start = 0;
    for (unsigned e = 0; e < this->Nx; e++)
    {
      // If we are past the
      if (e > 0)
      {
        n_start = 1;
      }
      // Get the layer of elements just above the interface
      FiniteElement* el_pt = this->finite_element_pt(this->Nx * this->Ny1 + e);
      // The first n_p nodes lie on the boundary
      for (unsigned n = n_start; n < n_p; n++)
      {
        Node* nod_pt = el_pt->node_pt(n);
        this->convert_to_boundary_node(nod_pt);
        this->add_boundary_node(6, nod_pt);
        b_coord[0] = nod_pt->x(0);
        nod_pt->set_coordinates_on_boundary(6, b_coord);
      }
    }

    // Allocate store for the spines:
    if (this->Xperiodic)
    {
      Spine_pt.reserve((n_p - 1) * this->Nx);
    }
    else
    {
      Spine_pt.reserve((n_p - 1) * this->Nx + 1);
    }


    // FIRST SPINE
    //-----------

    // Element 0
    // Node 0
    // Assign the new spine of height of the lower layer
    Spine* new_spine_pt = new Spine(H1);
    Spine_pt.push_back(new_spine_pt);

    // Get pointer to node
    SpineNode* nod_pt = element_node_pt(0, 0);

    // Set the pointer to the spine
    nod_pt->spine_pt() = new_spine_pt;
    // Set the fraction
    nod_pt->fraction() = 0.0;
    // Pointer to the mesh that implements the update fct
    nod_pt->spine_mesh_pt() = this;
    // ID of update function within this mesh: 0 = lower; 1 = upper
    nod_pt->node_update_fct_id() = 0;

    // Loop vertically along the spine
    // Loop over the elements in fluid 1
    for (unsigned long i = 0; i < Ny1; i++)
    {
      // Loop over the vertical nodes, apart from the first
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // Get pointer to node
        SpineNode* nod_pt = element_node_pt(i * this->Nx, l1 * n_p);
        // Set the pointer to the spine
        nod_pt->spine_pt() = new_spine_pt;
        // Set the fraction
        nod_pt->fraction() = (nod_pt->x(1) - this->Ymin) / (H1);
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this;
        // ID of update function within this mesh: 0 = lower; 1 = upper
        nod_pt->node_update_fct_id() = 0;
      }
    }

    // Loop over the elements in fluid 2
    for (unsigned long i = 0; i < Ny2; i++)
    {
      // Loop over vertical nodes, apart from the first
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // Get pointer to node
        SpineNode* nod_pt = element_node_pt((Ny1 + i) * this->Nx, l1 * n_p);

        // Set the pointer to the spine
        nod_pt->spine_pt() = new_spine_pt;
        // Set the fraction
        nod_pt->fraction() = (nod_pt->x(1) - (this->Ymin + H1)) / (H2);
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this;
        // ID of update function within this mesh: 0 = lower; 1 = upper
        nod_pt->node_update_fct_id() = 1;
      }
    }


    // LOOP OVER OTHER SPINES
    //----------------------

    // Now loop over the elements horizontally
    for (unsigned long j = 0; j < this->Nx; j++)
    {
      // Loop over the nodes in the elements horizontally, ignoring
      // the first column

      // Last spine needs special treatment in x-periodic meshes:
      unsigned n_pmax = n_p;
      if ((this->Xperiodic) && (j == this->Nx - 1)) n_pmax = n_p - 1;

      for (unsigned l2 = 1; l2 < n_pmax; l2++)
      {
        // Assign the new spine with length the height of the lower layer
        new_spine_pt = new Spine(H1);
        Spine_pt.push_back(new_spine_pt);

        // Get pointer to node
        SpineNode* nod_pt = element_node_pt(j, l2);

        // Set the pointer to the spine
        nod_pt->spine_pt() = new_spine_pt;
        // Set the fraction
        nod_pt->fraction() = 0.0;
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this;
        // ID of update function within this mesh: 0 = lower; 1 = upper
        nod_pt->node_update_fct_id() = 0;

        // Loop vertically along the spine
        // Loop over the elements in fluid 1
        for (unsigned long i = 0; i < Ny1; i++)
        {
          // Loop over the vertical nodes, apart from the first
          for (unsigned l1 = 1; l1 < n_p; l1++)
          {
            // Get pointer to node
            SpineNode* nod_pt =
              element_node_pt(i * this->Nx + j, l1 * n_p + l2);
            // Set the pointer to the spine
            nod_pt->spine_pt() = new_spine_pt;
            // Set the fraction
            nod_pt->fraction() = (nod_pt->x(1) - this->Ymin) / H1;
            // Pointer to the mesh that implements the update fct
            nod_pt->spine_mesh_pt() = this;
            // ID of update function within this mesh: 0 = lower; 1 = upper
            nod_pt->node_update_fct_id() = 0;
          }
        }

        // Loop over the elements in fluid 2
        for (unsigned long i = 0; i < Ny2; i++)
        {
          // Loop over vertical nodes, apart from the first
          for (unsigned l1 = 1; l1 < n_p; l1++)
          {
            // Get pointer to node
            SpineNode* nod_pt =
              element_node_pt((Ny1 + i) * this->Nx + j, l1 * n_p + l2);

            // Set the pointer to the spine
            nod_pt->spine_pt() = new_spine_pt;
            // Set the fraction
            nod_pt->fraction() = (nod_pt->x(1) - (this->Ymin + H1)) / H2;
            // Pointer to the mesh that implements the update fct
            nod_pt->spine_mesh_pt() = this;
            // ID of update function within this mesh: 0 = lower; 1 = upper
            nod_pt->node_update_fct_id() = 1;
          }
        }
      }
    }


    // Last spine needs special treatment for periodic meshes
    // because it's the same as the first one...
    if (this->Xperiodic)
    {
      // Last spine is the same as first one...
      Spine* final_spine_pt = Spine_pt[0];

      // Get pointer to node
      SpineNode* nod_pt = element_node_pt((this->Nx - 1), (n_p - 1));

      // Set the pointer to the spine
      nod_pt->spine_pt() = final_spine_pt;
      // Set the fraction to be the same as for the nodes on the first row
      nod_pt->fraction() = element_node_pt(0, 0)->fraction();
      // Pointer to the mesh that implements the update fct
      nod_pt->spine_mesh_pt() = element_node_pt(0, 0)->spine_mesh_pt();
      // ID of update function within this mesh: 0 = lower; 1 = upper
      nod_pt->node_update_fct_id() =
        element_node_pt(0, 0)->node_update_fct_id();

      // Now loop vertically along the spine
      for (unsigned i = 0; i < (Ny1 + Ny2); i++)
      {
        // Loop over the vertical nodes, apart from the first
        for (unsigned l1 = 1; l1 < n_p; l1++)
        {
          // Get pointer to node
          SpineNode* nod_pt = element_node_pt(i * this->Nx + (this->Nx - 1),
                                              l1 * n_p + (n_p - 1));

          // Set the pointer to the spine
          nod_pt->spine_pt() = final_spine_pt;
          // Set the fraction to be the same as in first row
          nod_pt->fraction() =
            element_node_pt(i * this->Nx, l1 * n_p)->fraction();
          // ID of update function within this mesh: 0 = lower; 1 = upper
          nod_pt->node_update_fct_id() =
            element_node_pt(i * this->Nx, l1 * n_p)->node_update_fct_id();
          // Pointer to the mesh that implements the update fct
          nod_pt->spine_mesh_pt() =
            element_node_pt(i * this->Nx, l1 * n_p)->spine_mesh_pt();
        }
      }
    }


    // Assign the 1D Line elements
    //---------------------------

    // Get the present number of elements
    /*unsigned long element_count = Element_pt.size();

    //Loop over the horizontal elements
    for(unsigned i=0;i<this->Nx;i++)
     {
     //Construct a new 1D line element on the face on which the local
     //coordinate 1 is fixed at its max. value (1) -- Face 2
      FiniteElement *interface_element_element_pt =
       new INTERFACE_ELEMENT(finite_element_pt(this->Nx*(Ny1-1)+i),2);

      //Push it back onto the stack
      Element_pt.push_back(interface_element_element_pt);

      //Push it back onto the stack of interface elements
      Interface_element_pt.push_back(interface_element_element_pt);

      element_count++;
      }*/

    // Setup the boundary information
    this->setup_boundary_element_info();
  }


  //======================================================================
  /// Reorder the elements, so we loop over them vertically first
  /// (advantageous in "wide" domains if a frontal solver is used).
  //======================================================================
  /*template <class ELEMENT, >
  void TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>::element_reorder()
  {

   unsigned Nx = this->Nx;
   //Find out how many elements there are
   unsigned long Nelement = nelement();
   //Find out how many fluid elements there are
   unsigned long Nfluid = Nx*(Ny1+Ny2);
   //Create a dummy array of elements
   Vector<FiniteElement *> dummy;

   //Loop over the elements in horizontal order
   for(unsigned long j=0;j<Nx;j++)
    {
     //Loop over the elements in lower layer vertically
     for(unsigned long i=0;i<Ny1;i++)
      {
       //Push back onto the new stack
       dummy.push_back(finite_element_pt(Nx*i + j));
      }

     //Push back the line element onto the stack
     dummy.push_back(finite_element_pt(Nfluid+j));

     //Loop over the elements in upper layer vertically
     for(unsigned long i=Ny1;i<(Ny1+Ny2);i++)
      {
       //Push back onto the new stack
       dummy.push_back(finite_element_pt(Nx*i + j));
      }
    }

   //Now copy the reordered elements into the element_pt
   for(unsigned long e=0;e<Nelement;e++)
    {
     Element_pt[e] = dummy[e];
    }

    }*/

} // namespace oomph
#endif
