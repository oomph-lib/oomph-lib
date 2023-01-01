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
#ifndef OOMPH_HORIZONTAL_SINGLE_LAYER_SPINE_MESH_TEMPLATE_CC
#define OOMPH_HORIZONTAL_SINGLE_LAYER_SPINE_MESH_TEMPLATE_CC

#include "horizontal_single_layer_spine_mesh.template.h"
#include "rectangular_quadmesh.template.cc"


namespace oomph
{
  //===========================================================================
  /// Constructor for spine 2D mesh: Pass number of elements in x-direction,
  /// number of elements in y-direction, axial length and height of layer,
  /// and pointer to timestepper (defaults to Static timestepper).
  ///
  /// The mesh contains a layer of spinified fluid elements (of type ELEMENT;
  /// e.g  SpineElement<QCrouzeixRaviartElement<2>)
  /// and information about how the internal nodes positions are affected
  /// by changes in spine length. Additional equations that determine the
  /// spine heights must be specified in order to use this mesh.
  //===========================================================================
  template<class ELEMENT>
  HorizontalSingleLayerSpineMesh<ELEMENT>::HorizontalSingleLayerSpineMesh(
    const unsigned& nx,
    const unsigned& ny,
    const double& lx,
    const double& h,
    TimeStepper* time_stepper_pt)
    : RectangularQuadMesh<ELEMENT>(
        nx, ny, 0.0, lx, 0.0, h, false, false, time_stepper_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Mesh can only be built with spine elements
    MeshChecker::assert_geometric_element<SpineFiniteElement, ELEMENT>(2);

    // We've called the "generic" constructor for the RectangularQuadMesh
    // which doesn't do much...

    // Now build the mesh:
    build_horizontal_single_layer_mesh(time_stepper_pt);
  }

  //===========================================================================
  /// Helper function that actually builds the single-layer spine mesh
  /// based on the parameters set in the various constructors
  //===========================================================================
  template<class ELEMENT>
  void HorizontalSingleLayerSpineMesh<
    ELEMENT>::build_horizontal_single_layer_mesh(TimeStepper* time_stepper_pt)
  {
    // Build the underlying quad mesh:
    RectangularQuadMesh<ELEMENT>::build_mesh(time_stepper_pt);

    // Read out the number of elements in the x-direction
    unsigned n_x = this->Nx;
    unsigned n_y = this->Ny;

    // Allocate memory for the spines and fractions along spines
    //---------------------------------------------------------

    // Read out number of linear points in the element
    unsigned n_p = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();
    Spine_pt.reserve((n_p - 1) * n_y + 1);

    // FIRST SPINE
    // -----------

    // Element 0
    // Node 0
    // Assign the new spine with unit length
    Spine* new_spine_pt = new Spine(1.0);
    Spine_pt.push_back(new_spine_pt);


    // Get pointer to node
    SpineNode* nod_pt = element_node_pt(0, 0);
    // Set the pointer to the spine
    nod_pt->spine_pt() = new_spine_pt;
    // Set the fraction
    nod_pt->fraction() = 0.0;
    // Pointer to the mesh that implements the update fct
    nod_pt->spine_mesh_pt() = this;

    // Loop HORIZONTAL along the spine
    // Loop over the elements
    for (unsigned long i = 0; i < n_x; i++)
    {
      // Loop over the HORIZONTAL nodes, apart from the first
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // Get pointer to node
        // SpineNode* nod_pt=element_node_pt(i*n_y,l1*n_p);

        // Get pointer to node(without reoder)
        SpineNode* nod_pt = element_node_pt(i, l1);
        // Set the pointer to the spine
        nod_pt->spine_pt() = new_spine_pt;
        // Set the fraction
        nod_pt->fraction() =
          (double(i) + double(l1) / double(n_p - 1)) / double(n_x);
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this;
      }
    }


    // LOOP OVER OTHER SPINES
    // ----------------------

    // Now loop over the elements VERTICALLY
    for (unsigned long j = 0; j < n_y; j++)
    {
      // Loop over the nodes in the elements horizontally, ignoring
      // the first row

      // Last spine needs special treatment in x-periodic meshes:
      unsigned n_pmax = n_p;

      for (unsigned l2 = 1; l2 < n_pmax; l2++)
      {
        // Assign the new spine with unit height
        new_spine_pt = new Spine(1.0);
        Spine_pt.push_back(new_spine_pt);

        // Get the node
        // SpineNode* nod_pt=element_node_pt(j,l2);

        // Get the node (without reorder)
        SpineNode* nod_pt = element_node_pt(j * n_x, l2 * n_p);
        // Set the pointer to spine
        nod_pt->spine_pt() = new_spine_pt;
        // Set the fraction
        nod_pt->fraction() = 0.0;
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this;

        // Loop HORIZONTALLY along the spine
        // Loop over the elements
        for (unsigned long i = 0; i < n_x; i++)
        {
          // Loop over the HORIZONTAL nodes, apart from the first
          for (unsigned l1 = 1; l1 < n_p; l1++)
          {
            // Get the node
            // SpineNode* nod_pt=element_node_pt(i*n_y+j,l1*n_p+l2);

            // Get the node (without reorder)
            SpineNode* nod_pt = element_node_pt(j * n_x + i, l2 * n_p + l1);
            // Set the pointer to the spine
            nod_pt->spine_pt() = new_spine_pt;
            // Set the fraction
            nod_pt->fraction() =
              (double(i) + double(l1) / double(n_p - 1)) / double(n_x);
            // Pointer to the mesh that implements the update fct
            nod_pt->spine_mesh_pt() = this;
          }
        }
      }
    }
  }


} // namespace oomph
#endif
