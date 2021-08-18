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
// Driver for 2D moving block
#ifndef OOMPH_QUAD_FROM_TRIANGLE_MESH_TEMPLATE_HEADER
#define OOMPH_QUAD_FROM_TRIANGLE_MESH_TEMPLATE_HEADER


#ifdef OOMPH_HAS_MPI
// MPI header
#include "mpi.h"
#endif

// Standards
#include <float.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <iomanip>

#ifdef OOMPH_HAS_FPUCONTROLH
#include <fpu_control.h>
#endif

// The mesh
#include "../generic/problem.h"
#include "../generic/quad_mesh.h"
#include "triangle_mesh.template.h"
#include "../generic/triangle_scaffold_mesh.h"
#include "../generic/unstructured_two_d_mesh_geometry_base.h"
#include "../generic/refineable_quad_mesh.h"
#include "../generic/Qelements.h"


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


namespace oomph
{
  //============start_of_quad_triangle_class==============================
  /// Quad mesh built on top of triangle scaffold mesh coming
  /// from the triangle mesh generator Triangle.
  /// http://www.cs.cmu.edu/~quake/triangle.html
  //======================================================================
  template<class ELEMENT>
  class QuadFromTriangleMesh : public virtual UnstructuredTwoDMeshGeometryBase,
                               public virtual QuadMeshBase
  {
  public:
    /// \short Empty constructor
    QuadFromTriangleMesh()
    {
#ifdef OOMPH_HAS_TRIANGLE_LIB
      // By default allow the automatic creation of vertices along the
      // boundaries by Triangle
      this->Allow_automatic_creation_of_vertices_on_boundaries = true;
#endif

      // Mesh can only be built with 2D Qelements.
      MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);
    }

    /// \short Constructor with the input files
    QuadFromTriangleMesh(
      const std::string& node_file_name,
      const std::string& element_file_name,
      const std::string& poly_file_name,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
      const bool& use_attributes = false,
      const bool& allow_automatic_creation_of_vertices_on_boundaries = true)
    {
      // Mesh can only be built with 2D Qelements.
      MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

      // Initialize the value for allowing creation of points on boundaries
      this->Allow_automatic_creation_of_vertices_on_boundaries =
        allow_automatic_creation_of_vertices_on_boundaries;

      // Store Timestepper used to build elements
      this->Time_stepper_pt = time_stepper_pt;

      // Store the attributes
      this->Use_attributes = use_attributes;

      // Build scaffold
      TriangleScaffoldMesh* tmp_mesh_pt = new TriangleScaffoldMesh(
        node_file_name, element_file_name, poly_file_name);

      // Convert mesh from scaffold to actual mesh
      this->build_from_scaffold(tmp_mesh_pt, time_stepper_pt, use_attributes);

      // Kill the scaffold
      delete tmp_mesh_pt;
      tmp_mesh_pt = 0;

      // Setup boundary coordinates for boundaries
      unsigned nbound = nboundary();
      for (unsigned ibound = 0; ibound < nbound; ibound++)
      {
        this->template setup_boundary_coordinates<ELEMENT>(ibound);
      }
    }

#ifdef OOMPH_HAS_TRIANGLE_LIB

    /// \short Build mesh, based on the specifications on
    /// TriangleMeshParameters. All the actual work is done
    /// in UnstructuredTwoDMeshGeometryBase
    QuadFromTriangleMesh(
      TriangleMeshParameters& triangle_mesh_parameters,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
    {
      // Mesh can only be built with 2D Qelements.
      MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

      // Initialize the value for allowing creation of points on boundaries
      this->Allow_automatic_creation_of_vertices_on_boundaries =
        triangle_mesh_parameters
          .is_automatic_creation_of_vertices_on_boundaries_allowed();

      // Store Timestepper used to build elements
      this->Time_stepper_pt = time_stepper_pt;

      // ********************************************************************
      // First part - Get polylines representations
      // ********************************************************************

      // Create the polyline representation of all the boundaries and
      // then create the mesh by calling to "generic_constructor()"

      // Initialise highest boundary id
      unsigned max_boundary_id = 0;

      // *****************************************************************
      // Part 1.1 - Outer boundary
      // *****************************************************************
      // Get the representation of the outer boundaries from the
      // TriangleMeshParameters object
      Vector<TriangleMeshClosedCurve*> outer_boundary_pt =
        triangle_mesh_parameters.outer_boundary_pt();

#ifdef PARANOID

      // Verify that the outer_boundary_object_pt has been set
      if (outer_boundary_pt.size() == 0)
      {
        std::stringstream error_message;
        error_message
          << "There are no outer boundaries defined.\n"
          << "Verify that you have specified the outer boundaries in the\n"
          << "Triangle_mesh_parameter object\n\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      } // if (outer_boundary_pt!=0)

#endif

      // Find the number of outer closed curves
      unsigned n_outer_boundaries = outer_boundary_pt.size();

      // Create the storage for the polygons that define the outer boundaries
      Vector<TriangleMeshPolygon*> outer_boundary_polygon_pt(
        n_outer_boundaries);

      // Loop over the number of outer boundaries
      for (unsigned i = 0; i < n_outer_boundaries; ++i)
      {
        // Get the polygon representation and compute the max boundary_id on
        // each outer polygon. Does nothing (i.e. just returns a pointer to
        // the outer boundary that was input) if the outer boundary is
        // already a polygon
        outer_boundary_polygon_pt[i] = this->closed_curve_to_polygon_helper(
          outer_boundary_pt[i], max_boundary_id);
      }

      // *****************************************************************
      // Part 1.2 - Internal closed boundaries (possible holes)
      // *****************************************************************
      // Get the representation of the internal closed boundaries from the
      // TriangleMeshParameters object
      Vector<TriangleMeshClosedCurve*> internal_closed_curve_pt =
        triangle_mesh_parameters.internal_closed_curve_pt();

      // Find the number of internal closed curves
      unsigned n_internal_closed_curves = internal_closed_curve_pt.size();

      // Create the storage for the polygons that define the internal closed
      // boundaries (again nothing happens (as above) if an internal closed
      // curve is already a polygon)
      Vector<TriangleMeshPolygon*> internal_polygon_pt(
        n_internal_closed_curves);

      // Loop over the number of internal closed curves
      for (unsigned i = 0; i < n_internal_closed_curves; ++i)
      {
        // Get the polygon representation and compute the max boundary_id on
        // each internal polygon
        internal_polygon_pt[i] = this->closed_curve_to_polygon_helper(
          internal_closed_curve_pt[i], max_boundary_id);
      }

      // *****************************************************************
      // Part 1.3 - Internal open boundaries
      // *****************************************************************
      // Get the representation of open boundaries from the
      // TriangleMeshParameteres object
      Vector<TriangleMeshOpenCurve*> internal_open_curve_pt =
        triangle_mesh_parameters.internal_open_curves_pt();

      // Find the number of internal open curves
      unsigned n_internal_open_curves = internal_open_curve_pt.size();

      // Create the storage for the polylines that define the open boundaries
      Vector<TriangleMeshOpenCurve*> internal_open_curve_poly_pt(
        n_internal_open_curves);

      // Loop over the number of internal open curves
      for (unsigned i = 0; i < n_internal_open_curves; i++)
      {
        // Get the open polyline representation and compute the max boundary_id
        // on each open polyline (again, nothing happens if there are curve
        // sections on the current internal open curve)
        internal_open_curve_poly_pt[i] =
          this->create_open_curve_with_polyline_helper(
            internal_open_curve_pt[i], max_boundary_id);
      }

      // ********************************************************************
      // Second part - Get associated geom objects and coordinate limits
      // ********************************************************************

      // ***************************************************************
      // Part 2.1 Outer boundary
      // ***************************************************************
      for (unsigned i = 0; i < n_outer_boundaries; i++)
      {
        this->set_geom_objects_and_coordinate_limits_for_close_curve(
          outer_boundary_pt[i]);
      }

      // ***************************************************************
      // Part 2.2 - Internal closed boundaries (possible holes)
      // ***************************************************************
      for (unsigned i = 0; i < n_internal_closed_curves; i++)
      {
        this->set_geom_objects_and_coordinate_limits_for_close_curve(
          internal_closed_curve_pt[i]);
      }

      // ********************************************************************
      // Part 2.3 - Internal open boundaries
      // ********************************************************************
      for (unsigned i = 0; i < n_internal_open_curves; i++)
      {
        this->set_geom_objects_and_coordinate_limits_for_open_curve(
          internal_open_curve_pt[i]);
      }

      // ********************************************************************
      // Third part - Creates the TriangulateIO object by calling the
      //              "generic_constructor()" function
      // ********************************************************************
      // Get all the other parameters from the TriangleMeshParameters object
      // The maximum element area
      const double element_area = triangle_mesh_parameters.element_area();

      // The holes coordinates
      Vector<Vector<double>> extra_holes_coordinates =
        triangle_mesh_parameters.extra_holes_coordinates();

      // The regions coordinates
      std::map<unsigned, Vector<double>> regions =
        triangle_mesh_parameters.regions_coordinates();

      // If we use regions then we use attributes
      const bool use_attributes = triangle_mesh_parameters.is_use_attributes();

      const bool refine_boundary =
        triangle_mesh_parameters.is_boundary_refinement_allowed();

      const bool refine_internal_boundary =
        triangle_mesh_parameters.is_internal_boundary_refinement_allowed();

      if (!refine_internal_boundary && refine_boundary)
      {
        std::ostringstream error_stream;
        error_stream
          << "You have specified that Triangle may refine the outer boundary, "
             "but\n"
          << "not internal boundaries. Triangle does not support this "
             "combination.\n"
          << "If you do not want Triangle to refine internal boundaries, it "
             "can't\n"
          << "refine outer boundaries either!\n"
          << "Please either disable all boundary refinement\n"
          << "(call TriangleMeshParameters::disable_boundary_refinement()\n"
          << "or enable internal boundary refinement (the default)\n";

        throw OomphLibError(error_stream.str().c_str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      this->generic_constructor(
        outer_boundary_polygon_pt,
        internal_polygon_pt,
        internal_open_curve_poly_pt,
        element_area,
        extra_holes_coordinates,
        regions,
        triangle_mesh_parameters.target_area_for_region(),
        time_stepper_pt,
        use_attributes,
        refine_boundary,
        refine_internal_boundary);

#ifdef OOMPH_HAS_MPI

      // Before calling setup boundary coordinates check if the mesh is
      // marked as distrbuted
      if (triangle_mesh_parameters.is_mesh_distributed())
      {
        // Set the mesh as distributed by passing the communicator
        this->set_communicator_pt(triangle_mesh_parameters.communicator_pt());
      }

#endif

      // Setup boundary coordinates for boundaries
      unsigned nb = nboundary();
      for (unsigned b = 0; b < nb; b++)
      {
        this->template setup_boundary_coordinates<ELEMENT>(b);
      }

      // Snap it!
      this->snap_nodes_onto_geometric_objects();
    }

    /// \short A general-purpose construction function that builds the
    /// mesh once the different specific constructors have assembled the
    /// appropriate information.
    void generic_constructor(
      Vector<TriangleMeshPolygon*>& outer_boundary_pt,
      Vector<TriangleMeshPolygon*>& internal_polygon_pt,
      Vector<TriangleMeshOpenCurve*>& open_polylines_pt,
      const double& element_area,
      Vector<Vector<double>>& extra_holes_coordinates,
      std::map<unsigned, Vector<double>>& regions_coordinates,
      std::map<unsigned, double>& regions_areas,
      TimeStepper* time_stepper_pt,
      const bool& use_attributes,
      const bool& refine_boundary,
      const bool& refine_internal_boundary)
    {
#ifdef PARANOID

      if (element_area < 10e-14)
      {
        std::ostringstream warning_message;
        warning_message
          << "The current elements area was stated to (" << element_area
          << ").\nThe current precision to generate the input to triangle "
          << "is fixed to 14 digits\n\n";
        OomphLibWarning(warning_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
      }

#endif

      // Store the attribute flag
      this->Use_attributes = use_attributes;

      // Store the timestepper used to build elements
      this->Time_stepper_pt = time_stepper_pt;

      // Store outer polygon
      this->Outer_boundary_pt = outer_boundary_pt;

      // Store internal polygons by copy constructor
      this->Internal_polygon_pt = internal_polygon_pt;

      // Store internal polylines by copy constructor
      this->Internal_open_curve_pt = open_polylines_pt;

      // Store the extra holes coordinates
      this->Extra_holes_coordinates = extra_holes_coordinates;

      // Store the extra regions coordinates
      this->Regions_coordinates = regions_coordinates;

      // Create the data structures required to call the triangulate function
      TriangulateIO triangulate_io;
      TriangulateIO triangulate_out;

      // Initialize TriangulateIO structure
      TriangleHelper::initialise_triangulateio(triangulate_io);

      // Convert TriangleMeshPolyLine and TriangleMeshClosedCurvePolyLine
      // to a triangulateio object
      UnstructuredTwoDMeshGeometryBase::build_triangulateio(
        outer_boundary_pt,
        internal_polygon_pt,
        open_polylines_pt,
        extra_holes_coordinates,
        regions_coordinates,
        regions_areas,
        triangulate_io);

      // Initialize TriangulateIO structure
      TriangleHelper::initialise_triangulateio(triangulate_out);

      // Input string for triangle
      std::stringstream input_string_stream;
      input_string_stream.precision(14);
      input_string_stream.setf(std::ios_base::fixed, std::ios_base::floatfield);

      // MH: Used to be:
      // input_string_stream<<"-pA -a" << element_area << " -q30" << std::fixed;
      // The repeated -a allows the specification of areas for different
      // regions (if any)
      input_string_stream << "-pA -a -a" << element_area << " -q30"
                          << std::fixed;

      // Verify if creation of new points on boundaries is allowed
      if (!this->is_automatic_creation_of_vertices_on_boundaries_allowed())
      {
        input_string_stream << " -YY";
      }

      // Suppress insertion of additional points on outer boundary
      if (refine_boundary == false)
      {
        input_string_stream << "-Y";

        // Add the extra flag to suppress additional points on interior segments
        if (refine_internal_boundary == false)
        {
          input_string_stream << "Y";
        }
      }

      // Convert the Input string in *char required by the triangulate function
      char triswitches[100];
      sprintf(triswitches, "%s", input_string_stream.str().c_str());

      // Build the mesh using triangulate function
      triangulate(triswitches, &triangulate_io, &triangulate_out, 0);

#ifdef OOMPH_HAS_FPUCONTROLH
      // Reset flags that are tweaked by triangle; can cause nasty crashes
      fpu_control_t cw = (_FPU_DEFAULT & ~_FPU_EXTENDED) | _FPU_DOUBLE;
      _FPU_SETCW(cw);
#endif

      // Build scaffold
      TriangleScaffoldMesh* tmp_mesh_pt =
        new TriangleScaffoldMesh(triangulate_out);

      // If we have filled holes then we must use the attributes
      if (!regions_coordinates.empty())
      {
        // Convert mesh from scaffold to actual mesh
        build_from_scaffold(tmp_mesh_pt, time_stepper_pt, true);

        // Record the attribute flag
        this->Use_attributes = true;
      }
      // Otherwise use what was asked
      else
      {
        // Convert mesh from scaffold to actual mesh
        build_from_scaffold(tmp_mesh_pt, time_stepper_pt, use_attributes);
      }

      // Kill the scaffold
      delete tmp_mesh_pt;
      tmp_mesh_pt = 0;

      // Cleanup but leave hole and regions alone since it's still used
      bool clear_hole_data = false;
      TriangleHelper::clear_triangulateio(triangulate_io, clear_hole_data);
      TriangleHelper::clear_triangulateio(triangulate_out, clear_hole_data);
    }

#endif // OOMPH_HAS_TRIANGLE_LIB

    /// Broken copy constructor
    QuadFromTriangleMesh(const QuadFromTriangleMesh& dummy)
    {
      BrokenCopy::broken_copy("QuadFromTriangleMesh");
    }

    /// Broken assignment operator
    void operator=(const QuadFromTriangleMesh&)
    {
      BrokenCopy::broken_assign("QuadFromTriangleMesh");
    }


    /// Empty destructor
    ~QuadFromTriangleMesh()
    {
#ifdef OOMPH_HAS_TRIANGLE_LIB

      std::set<TriangleMeshCurveSection*>::iterator it_polyline;
      for (it_polyline = Free_curve_section_pt.begin();
           it_polyline != Free_curve_section_pt.end();
           it_polyline++)
      {
        delete (*it_polyline);
      }

      std::set<TriangleMeshPolygon*>::iterator it_polygon;
      for (it_polygon = Free_polygon_pt.begin();
           it_polygon != Free_polygon_pt.end();
           it_polygon++)
      {
        delete (*it_polygon);
      }

      std::set<TriangleMeshOpenCurve*>::iterator it_open_polyline;
      for (it_open_polyline = Free_open_curve_pt.begin();
           it_open_polyline != Free_open_curve_pt.end();
           it_open_polyline++)
      {
        delete (*it_open_polyline);
      }

#endif
    }

    /// Build the quad mesh from the given scaffold mesh
    void build_from_scaffold(TriangleScaffoldMesh* tmp_mesh_pt,
                             TimeStepper* time_stepper_pt,
                             const bool& use_attributes);

    /// Timestepper used to build elements
    TimeStepper* Time_stepper_pt;

    /// Boolean flag to indicate whether to use attributes or not (required
    /// for multidomain meshes)
    bool Use_attributes;
  };


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //=========================================================================
  /// Unstructured refineable QuadFromTriangleMesh
  //=========================================================================
  template<class ELEMENT>
  class RefineableQuadFromTriangleMesh
    : public virtual QuadFromTriangleMesh<ELEMENT>,
      public virtual RefineableQuadMesh<ELEMENT>
  {
  public:
#ifdef OOMPH_HAS_TRIANGLE_LIB

    /// \short Build mesh, based on the specifications on
    /// TriangleMeshParameters
    RefineableQuadFromTriangleMesh(
      TriangleMeshParameters& triangle_mesh_parameters,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : QuadFromTriangleMesh<ELEMENT>(triangle_mesh_parameters, time_stepper_pt)
    {
      this->setup_quadtree_forest();
    }

#endif

    /// Refine mesh uniformly
    virtual void refine_uniformly()
    {
      DocInfo doc_info;
      doc_info.directory() = "";
      doc_info.disable_doc();
      refine_uniformly(doc_info);
    }

    /// Refine mesh uniformly and doc process
    void refine_uniformly(DocInfo& doc_info)
    {
      // Find the number of elements in the mesh
      unsigned nelem = this->nelement();

      // Set the element error to something big
      Vector<double> elem_error(nelem, DBL_MAX);

      // Refine everything
      adapt(elem_error);
    }

    /// Overload the adapt function (to ensure nodes are snapped to the
    /// boundary)
    void adapt(const Vector<double>& elem_error);

    /// \short Build mesh, based on the polyfiles
    RefineableQuadFromTriangleMesh(
      const std::string& node_file_name,
      const std::string& element_file_name,
      const std::string& poly_file_name,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : QuadFromTriangleMesh<ELEMENT>(
          node_file_name, element_file_name, poly_file_name, time_stepper_pt)
    {
      this->setup_quadtree_forest();
    }

    /// Empty Destructor
    virtual ~RefineableQuadFromTriangleMesh() {}
  };


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //=========================================================================
  /// Unstructured QuadFromTriangleMesh upgraded to solid mesh
  //=========================================================================
  template<class ELEMENT>
  class SolidQuadFromTriangleMesh
    : public virtual QuadFromTriangleMesh<ELEMENT>,
      public virtual SolidMesh
  {
  public:
    SolidQuadFromTriangleMesh(
      const std::string& node_file_name,
      const std::string& element_file_name,
      const std::string& poly_file_name,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
      const bool& use_attributes = false)
      : QuadFromTriangleMesh<ELEMENT>(node_file_name,
                                      element_file_name,
                                      poly_file_name,
                                      time_stepper_pt,
                                      use_attributes)
    {
      // Assign the Lagrangian coordinates
      set_lagrangian_nodal_coordinates();
    }

#ifdef OOMPH_HAS_TRIANGLE_LIB

    /// \short Build mesh, based on closed curve that specifies
    /// the outer boundary of the domain and any number of internal
    /// clsed curves. Specify target area for uniform element size.
    SolidQuadFromTriangleMesh(
      TriangleMeshParameters& triangle_mesh_parameters,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : QuadFromTriangleMesh<ELEMENT>(triangle_mesh_parameters, time_stepper_pt)
    {
      // Assign the Lagrangian coordinates
      set_lagrangian_nodal_coordinates();
    }

#endif

    /// Empty Destructor
    virtual ~SolidQuadFromTriangleMesh() {}
  };


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //=========================================================================
  /// Unstructured refineable QuadFromTriangleMesh upgraded to solid mesh
  //=========================================================================
  template<class ELEMENT>
  class RefineableSolidQuadFromTriangleMesh
    : public virtual RefineableQuadFromTriangleMesh<ELEMENT>,
      public virtual SolidMesh
  {
  public:
    /// \short Build mesh from specified triangulation and associated
    /// target areas for elements in it.
    RefineableSolidQuadFromTriangleMesh(
      const std::string& node_file_name,
      const std::string& element_file_name,
      const std::string& poly_file_name,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
      const bool& use_attributes = false)
      : RefineableQuadFromTriangleMesh<ELEMENT>(node_file_name,
                                                element_file_name,
                                                poly_file_name,
                                                time_stepper_pt,
                                                use_attributes)
    {
      // Assign the Lagrangian coordinates
      set_lagrangian_nodal_coordinates();
    }

#ifdef OOMPH_HAS_TRIANGLE_LIB

    /// \short Build mesh, based on the specifications on
    /// TriangleMeshParameter
    RefineableSolidQuadFromTriangleMesh(
      TriangleMeshParameters& triangle_mesh_parameters,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : QuadFromTriangleMesh<ELEMENT>(triangle_mesh_parameters,
                                      time_stepper_pt),
        RefineableQuadFromTriangleMesh<ELEMENT>(triangle_mesh_parameters,
                                                time_stepper_pt)
    {
      // Assign the Lagrangian coordinates
      set_lagrangian_nodal_coordinates();
    }

#endif

    /// Empty Destructor
    virtual ~RefineableSolidQuadFromTriangleMesh() {}
  };


} // namespace oomph


#endif // OOMPH_QUAD_FROM_TRIANGLE_MESH_HEADER
