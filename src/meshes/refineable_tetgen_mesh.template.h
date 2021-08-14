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

#ifndef OOMPH_REFINEABLE_TETGEN_MESH_HEADER
#define OOMPH_REFINEABLE_TETGEN_MESH_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


#include "../generic/tetgen_scaffold_mesh.h"
#include "../generic/tet_mesh.h"
#include "../generic/refineable_mesh.h"
#include "tetgen_mesh.template.h"

namespace oomph
{
  //=========================================================================
  // Unstructured refineable TetgenMesh
  //=========================================================================
  template<class ELEMENT>
  class RefineableTetgenMesh : public virtual TetgenMesh<ELEMENT>,
                               public virtual RefineableTetMeshBase
  {
  public:
    /// \short Build mesh, based on a TetMeshFacetedClosedSurface that specifies
    /// the outer boundary of the domain and any number of internal
    /// closed curves, specified by TetMeshFacetedSurfaces.
    /// Also specify target volume for uniform element size.
    RefineableTetgenMesh(
      TetMeshFacetedClosedSurface* const& outer_boundary_pt,
      Vector<TetMeshFacetedSurface*>& internal_closed_surface_pt,
      const double& element_volume,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
      const bool& use_attributes = false,
      const bool& split_corner_elements = false)
      : TetgenMesh<ELEMENT>(outer_boundary_pt,
                            internal_closed_surface_pt,
                            element_volume,
                            time_stepper_pt,
                            use_attributes,
                            split_corner_elements),
        Corner_elements_must_be_split(split_corner_elements)
    {
      // Initialise the data associated with adaptation
      initialise_adaptation_data();
    }


  protected:
    /// \short Specialised constructor used during adaptation only.
    /// Element sizes are specified by vector tetgen_io is passed in
    /// from previous mesh (is then modified to build new mesh)
    /// Ditto with use_attributes, which comes from the previous mesh
    RefineableTetgenMesh(
      const Vector<double>& target_volume,
      tetgenio* const& tetgen_io_pt,
      TetMeshFacetedClosedSurface* const& outer_boundary_pt,
      Vector<TetMeshFacetedSurface*>& internal_surface_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
      const bool& use_attributes = false)
    {
      // NOTE THERE IS A CERTAIN AMOUNT OF DUPLICATION BETWEEN THE
      // CODE IN HERE AND THE ONE IN THE CONSTRUCTOR OF THE TetgenMesh
      // BUT THE FACT THAT WE HAVE TO MODIFY THE TETGENIO STRUCTURE
      // MEANS WE CAN'T QUITE RECYCLE THIS.

      // Mesh can only be built with 3D Telements.
      MeshChecker::assert_geometric_element<TElementGeometricBase, ELEMENT>(3);

      // Initialise the data associated with adaptation
      initialise_adaptation_data();

      // Store Timestepper used to build elements
      this->Time_stepper_pt = time_stepper_pt;

      // Triangulation has been created -- remember to wipe it!
      this->Tetgenio_exists = true;
      this->Tetgenio_pt = new tetgenio;

      // Add the volume constraints to the tetgenio data object.
      // Note that since the tetgenio structure is referred to by pointer
      // we're also modifying the one associated with the (still existing)
      // original mesh. Bit naughty but shouldn't cause any problems
      // since that mesh is already built and about to go out of scope
      // anyway.

      // Create a local copy
      tetgenio* tetgen_input_pt = new tetgenio;
      ;
      this->deep_copy_of_tetgenio(tetgen_io_pt, tetgen_input_pt);

      // Add volume constraints
      tetgen_input_pt->tetrahedronvolumelist =
        new double[tetgen_input_pt->numberoftetrahedra];
      for (int e = 0; e < tetgen_input_pt->numberoftetrahedra; ++e)
      {
        tetgen_input_pt->tetrahedronvolumelist[e] = target_volume[e];
      }

      // Input string
      std::stringstream input_string_stream;
      input_string_stream << "Vqra";

      // Convert to a *char
      char tetswitches[100];
      sprintf(tetswitches, "%s", input_string_stream.str().c_str());

      // Build triangulateio refined object
      tetrahedralize(tetswitches, tetgen_input_pt, this->Tetgenio_pt);
      // Build scaffold
      this->Tmp_mesh_pt = new TetgenScaffoldMesh(*this->Tetgenio_pt);

      // Convert mesh from scaffold to actual mesh
      this->build_from_scaffold(time_stepper_pt, use_attributes);

      // Kill the scaffold
      delete this->Tmp_mesh_pt;
      this->Tmp_mesh_pt = 0;

      // delete the input
      delete tetgen_input_pt;

      // Store the boundary
      this->Outer_boundary_pt = outer_boundary_pt;
      // Setup the reverse lookup scheme
      this->setup_reverse_lookup_schemes_for_faceted_surface(
        this->Outer_boundary_pt);
      // Store the internal boundary
      this->Internal_surface_pt = internal_surface_pt;
      // Setup the reverse lookup schemes
      {
        unsigned n = this->Internal_surface_pt.size();
        for (unsigned i = 0; i < n; i++)
        {
          this->setup_reverse_lookup_schemes_for_faceted_surface(
            this->Internal_surface_pt[i]);
        }
      }

      // Setup boundary coordinates for boundaries
      unsigned nb = nboundary();
      for (unsigned b = 0; b < nb; b++)
      {
        this->template setup_boundary_coordinates<ELEMENT>(b);
      }

      // Now snap onto geometric objects associated with triangular facets
      // (if any!)
      this->snap_nodes_onto_geometric_objects();
    }


  public:
    /// Empty Destructor
    virtual ~RefineableTetgenMesh() {}

    /// Refine mesh uniformly and doc process
    void refine_uniformly(DocInfo& doc_info)
    {
      // hierher do it
      throw OomphLibError("refine_uniformly() not implemented yet",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    /// \short Unrefine mesh uniformly: Return 0 for success,
    /// 1 for failure (if unrefinement has reached the coarsest permitted
    /// level)
    unsigned unrefine_uniformly()
    {
      // hierher do it
      throw OomphLibError("unrefine_uniformly() not implemented yet",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
      // dummy return
      return 0;
    }

    /// Adapt mesh, based on elemental error provided
    void adapt(const Vector<double>& elem_error);


    /// Is projection of old solution onto new mesh disabled?
    bool projection_is_disabled()
    {
      return Projection_is_disabled;
    }

    /// Disable projection of old solution onto new mesh
    void disable_projection()
    {
      Projection_is_disabled = true;
    }

    /// Disable projection of old solution onto new mesh
    void enable_projection()
    {
      Projection_is_disabled = false;
    }


  protected:
    /// Helper function to initialise data associated with adaptation
    void initialise_adaptation_data()
    {
      // Set max and min targets for adaptation
      this->Max_element_size = 1.0;
      this->Min_element_size = 0.001;
      this->Max_permitted_edge_ratio = 2.0;

      /// By default we project solution onto new mesh during adaptation
      Projection_is_disabled = false;
    }

    // Update the surface
    void update_faceted_surface_using_face_mesh(
      TetMeshFacetedSurface*& faceted_surface_pt);

    // Update the inner hole
    void surface_remesh_for_inner_hole_boundaries();

    /// \short Snap the boundary nodes onto any curvilinear boundaries
    void snap_nodes_onto_boundary(RefineableTetgenMesh<ELEMENT>*& new_mesh_pt,
                                  const unsigned& b);

    /// Disable projection of solution onto new mesh during adaptation
    bool Projection_is_disabled;

    /// \short Corner elements which have all of their nodes on the outer
    /// boundary are to be split into elements which have some non-boundary
    /// nodes
    bool Corner_elements_must_be_split;
  };

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  //=========================================================================
  // Unstructured refineable Tetgen Mesh upgraded to solid mesh
  //=========================================================================
  template<class ELEMENT>
  class RefineableSolidTetgenMesh
    : public virtual RefineableTetgenMesh<ELEMENT>,
      public virtual SolidMesh
  {
  public:
    /// \short Build mesh, based on closed curve that specifies
    /// the outer boundary of the domain and any number of internal
    /// closed curves. Specify target area for uniform element size.
    RefineableSolidTetgenMesh(
      TetMeshFacetedClosedSurface* const& outer_boundary_pt,
      Vector<TetMeshFacetedSurface*>& internal_closed_surface_pt,
      const double& element_volume,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
      const bool& use_attributes = false,
      const bool& split_corner_elements = false)
      : TetgenMesh<ELEMENT>(outer_boundary_pt,
                            internal_closed_surface_pt,
                            element_volume,
                            time_stepper_pt,
                            use_attributes,
                            split_corner_elements),
        RefineableTetgenMesh<ELEMENT>(outer_boundary_pt,
                                      internal_closed_surface_pt,
                                      element_volume,
                                      time_stepper_pt,
                                      use_attributes,
                                      split_corner_elements)

    {
      // Assign the Lagrangian coordinates
      set_lagrangian_nodal_coordinates();
    }


    /// \short Build mesh from specified triangulation and
    /// associated target areas for elements in it.
    RefineableSolidTetgenMesh(
      const Vector<double>& target_volume,
      tetgenio* const& tetgen_io_pt,
      TetMeshFacetedClosedSurface* const& outer_boundary_pt,
      Vector<TetMeshFacetedSurface*>& internal_surface_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
      const bool& use_attributes = false)
      : RefineableTetgenMesh<ELEMENT>(target_volume,
                                      tetgen_io_pt,
                                      outer_boundary_pt,
                                      internal_surface_pt,
                                      time_stepper_pt,
                                      use_attributes)

    {
      // Assign the Lagrangian coordinates
      set_lagrangian_nodal_coordinates();
    }

    /// Empty Destructor
    virtual ~RefineableSolidTetgenMesh() {}
  };


} // namespace oomph

#endif
