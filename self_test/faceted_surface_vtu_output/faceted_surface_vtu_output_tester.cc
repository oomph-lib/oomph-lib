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

// Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/triangle_mesh.h"

// Include definitions of faceted surfaces
#include "tetmesh_faceted_surfaces.h"

// Include the definition of SphericalTetMeshFacetedSurface
#include "spherical_tetmesh_faceted_surface.h"

using namespace oomph;

int main()
{
  // Introduce a dummy ID since we're not concerned with boundary IDs.
  unsigned dummy_id = 1;


  // CubicTetMeshFacetedSurface
  //----------------------------
  // Create a cubic faceted surface
  double box_half_width = 1.0;
  double box_half_length = 2.0;
  CubicTetMeshFacetedSurface cubic_faceted_surface(
    box_half_width, box_half_length, dummy_id);
  // Output into the Paraview .vtu format
  cubic_faceted_surface.output_paraview("cubic_faceted.vtu");


  // DiskTetMeshFacetedSurface
  //---------------------------
  // Create a faceted disk with an oscillating displacement
  double epsilon = 0.1; // Set the amplitude
  unsigned n = 4; // Set the mode
  double z_offset = 0.0; // Set the offset in the z coordinate
  // Set the number of lines discretising half of the disk's arc
  unsigned half_nsegment = 10;

  // Create the Geometric object
  WarpedCircularDisk disk(epsilon, n, z_offset);

  // Create the faceted surface
  DiskTetMeshFacetedSurface disk_faceted_surface(
    &disk, half_nsegment, dummy_id, dummy_id);
  // Output into the Paraview .vtu format
  disk_faceted_surface.output_paraview("disk_faceted.vtu");


  // SphericalTetMeshFacetedSurface
  //--------------------------------
  // Create a spherical faceted surface with a fixed number of discrete faces
  SphericalTetMeshFacetedSurface spherical_faceted_surface;
  // Output into the Paraview .vtu format
  spherical_faceted_surface.output_paraview("spherical_faceted.vtu");


  // RectangularTetMeshFacetedSurface
  //----------------------------------
  // Create a simple rectangular faceted surface consisting of a single facet.
  double half_x_width = 1.0;
  double half_y_length = 2.0;
  Vector<double> offset(3, 0.0); // Set the rectangle at the origin
  RectangularTetMeshFacetedSurface rectangular_faceted_surface(
    half_x_width, half_y_length, offset, dummy_id);
  // Output into the Paraview .vtu format
  rectangular_faceted_surface.output_paraview("rectangular_faceted.vtu");


  // DiskWithTwoLayersTetMeshFacetedSurface
  //----------------------------------------
  // Create a disk with two layers, looks like a coin
  DiskWithTwoLayersTetMeshFacetedSurface disk_with_two_layers_faceted_surface(
    &disk, half_nsegment, dummy_id, dummy_id, dummy_id, dummy_id, dummy_id);
  // Output into the Paraview .vtu format
  disk_with_two_layers_faceted_surface.output_paraview(
    "disk_with_two_layers_faceted.vtu");
}