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
// Contains the definition of a TriangulateIO object. This is used to
// define the complex geometry of a two-dimensional mesh which is why
// it resides here. The definition of things like TriangleMeshPolygons
// and other classes which define geometrical aspects of a 2D mesh can
// also be found here. The class UnstructuredTwoDMeshGeometryBase is
// defined here. It forms the base class for QuadFromTriangleMesh and
// TriangleMeshBase. This class makes use of previously written code
// to create TriangleScaffoldMeshes and avoid a large amount of code
// duplication.
#ifndef OOMPH_UNSTRUCTURED_TWO_D_MESH_GEOMETRY_BASE_HEADER
#define OOMPH_UNSTRUCTURED_TWO_D_MESH_GEOMETRY_BASE_HEADER

// The scaffold mesh
#include "mesh.h"

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


namespace oomph
{
#ifdef OOMPH_HAS_TRIANGLE_LIB

  //=====================================================================
  /// The Triangle data structure, modified from the triangle.h header
  /// supplied with triangle 1.6. by J. R. Schewchuk. We need to define
  /// this here separately because we can't include a c header directly
  /// into C++ code!
  //=====================================================================
  struct TriangulateIO
  {
    /// Pointer to list of points x coordinate followed by y coordinate
    double* pointlist;

    /// Pointer to list of point attributes
    double* pointattributelist;

    /// Pointer to list of point markers
    int* pointmarkerlist;
    int numberofpoints;
    int numberofpointattributes;

    int* trianglelist;
    double* triangleattributelist;
    double* trianglearealist;
    int* neighborlist;
    int numberoftriangles;
    int numberofcorners;
    int numberoftriangleattributes;

    int* segmentlist;
    int* segmentmarkerlist;
    int numberofsegments;

    double* holelist;
    int numberofholes;

    double* regionlist;
    int numberofregions;

    int* edgelist;
    int* edgemarkerlist; // <---- contains boundary ID (offset by one)
    double* normlist;
    int numberofedges;
  };


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //==================================================================
  /// Helper namespace for triangle meshes
  //==================================================================
  namespace TriangleHelper
  {
    /// Clear TriangulateIO structure
    extern void clear_triangulateio(TriangulateIO& triangulate_io,
                                    const bool& clear_hole_data = true);

    /// Initialise TriangulateIO structure
    extern void initialise_triangulateio(TriangulateIO& triangle_io);

    /// \short Make (partial) deep copy of TriangulateIO object. We only copy
    /// those items we need within oomph-lib's adaptation procedures.
    /// Warnings are issued if triangulate_io contains data that is not
    /// not copied, unless quiet=true;
    extern TriangulateIO deep_copy_of_triangulateio_representation(
      TriangulateIO& triangle_io, const bool& quiet);

    /// \short Write the triangulateio data to disk as a poly file,
    /// mainly used for debugging
    extern void write_triangulateio_to_polyfile(TriangulateIO& triangle_io,
                                                std::ostream& poly_file);

    /// Create a triangulateio data file from ele node and poly
    /// files. This is used if the mesh is generated by using Triangle
    /// externally. The triangulateio structure is required to dump the mesh
    /// topology for restarts.
    extern void create_triangulateio_from_polyfiles(
      const std::string& node_file_name,
      const std::string& element_file_name,
      const std::string& poly_file_name,
      TriangulateIO& triangle_io,
      bool& use_attributes);

    /// \short Write all the triangulateio data to disk in a dump file
    /// that can then be used to restart simulations
    extern void dump_triangulateio(TriangulateIO& triangle_io,
                                   std::ostream& dump_file);

    /// \short Read the triangulateio data from a dump file on
    /// disk, which can then be used to restart simulations
    extern void read_triangulateio(std::istream& restart_file,
                                   TriangulateIO& triangle_io);

  } // namespace TriangleHelper

#endif


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////


  class TriangleMeshPolyLine;
  class TriangleMeshCurviLine;

  //=====================================================================
  /// Base class for defining a triangle mesh boundary, this class has the
  /// methods that allow to connect the initial and final ends to other
  /// triangle mesh boundaries
  //=====================================================================
  class TriangleMeshCurveSection
  {
  public:
    /// Empty constructor. Initialises the curve section as non connected
    TriangleMeshCurveSection()
      : Initial_vertex_connected(false),
        Final_vertex_connected(false),
        Initial_vertex_connected_suspended(false),
        Final_vertex_connected_suspended(false),
        Initial_vertex_connected_to_curviline(false),
        Final_vertex_connected_to_curviline(false),
        Refinement_tolerance(0.08),
        Unrefinement_tolerance(0.04),
        Maximum_length(-1.0)
    {
    }

    /// Empty destructor
    virtual ~TriangleMeshCurveSection() {}

    /// \short Number of segments that this part of the
    /// boundary is to be represented by. This corresponds
    /// to the number of straight-line segments in triangle
    /// representation.
    virtual unsigned nsegment() const = 0;

    /// Boundary id
    virtual unsigned boundary_id() const = 0;

    /// Boundary chunk (Used when a boundary is represented by more
    /// than one polyline
    virtual unsigned boundary_chunk() const = 0;

    /// Number of vertices
    virtual unsigned nvertex() const = 0;

    /// Output the curve_section
    virtual void output(std::ostream& outfile,
                        const unsigned& n_sample = 50) = 0;

    /// \short Enable refinement of curve section to create a better
    /// representation of curvilinear boundaries (e.g. in free-surface
    /// problems). See tutorial for
    /// interpretation of the optional argument which specifies the
    /// refinement tolerance. It defaults to 0.08 and the smaller the
    /// number the finer the surface representation.
    void enable_refinement_tolerance(const double& tolerance = 0.08)
    {
      Refinement_tolerance = tolerance;
    }

    /// \short Set tolerance for refinement of curve sections to create a better
    /// representation of curvilinear boundaries (e.g. in free-surface
    /// problems). See tutorial for
    /// interpretation of the refinement tolerance. (The smaller the
    /// number the finer the surface representation). If set to
    /// a negative value, we're switching off refinement --
    /// equivalent to calling disable_polyline_refinement()
    void set_refinement_tolerance(const double& tolerance)
    {
      Refinement_tolerance = tolerance;
    }

    /// \short Get tolerance for refinement of curve sections to create a better
    /// representation of curvilinear boundaries (e.g. in free-surface
    /// problems). See tutorial for
    /// interpretation. If it's negative refinement is disabled.
    double refinement_tolerance()
    {
      return Refinement_tolerance;
    }

    /// \short Disable refinement of curve section
    void disable_refinement_tolerance()
    {
      Refinement_tolerance = -1.0;
    }

    /// \short Enable unrefinement of curve sections to avoid unnecessarily
    /// large numbers of elements on of curvilinear boundaries (e.g. in
    /// free-surface problems). See tutorial for interpretation of the optional
    /// argument which specifies the unrefinement tolerance. It defaults to 0.04
    /// and the larger the number the more agressive we are when removing
    /// unnecessary vertices on gently curved polylines.
    void enable_unrefinement_tolerance(const double& tolerance = 0.04)
    {
      Unrefinement_tolerance = tolerance;
    }

    /// \short Set tolerance for unrefinement of curve sections
    /// to avoid unnecessarily large
    /// numbers of elements on of curvilinear boundaries (e.g. in free-surface
    /// problems). See tutorial for
    /// interpretation of the optional argument which specifies the
    /// unrefinement tolerance. It defaults to 0.04 and the larger the number
    /// the more agressive we are when removing unnecessary vertices on
    /// gently curved polylines. If set to
    /// a negative value, we're switching off unrefinement --
    /// equivalent to calling disable_curve_section_unrefinement()
    void set_unrefinement_tolerance(const double& tolerance)
    {
      Unrefinement_tolerance = tolerance;
    }

    /// \short Get tolerance for unrefinement of curve section to create a
    /// better representation of curvilinear boundaries (e.g. in free-surface
    /// problems). See tutorial for
    /// interpretation. If it's negative unrefinement is disabled.
    double unrefinement_tolerance()
    {
      return Unrefinement_tolerance;
    }

    /// \short Disable unrefinement of curve sections
    void disable_unrefinement_tolerance()
    {
      Unrefinement_tolerance = -1.0;
    }

    /// \short Allows to specify the maximum distance between two vertices
    /// that define the associated polyline of the curve section, it only
    /// takes effect on the unrefinement and refinement steps
    void set_maximum_length(const double& maximum_length)
    {
      Maximum_length = maximum_length;
    }

    /// \short Disables the use of the maximum length criteria on the
    /// unrefinement or refinement steps
    void disable_use_maximum_length()
    {
      Maximum_length = -1.0;
    }

    /// \short Gets access to the maximum length variable
    double maximum_length()
    {
      return Maximum_length;
    }

    /// Get first vertex coordinates
    virtual void initial_vertex_coordinate(Vector<double>& vertex) = 0;

    /// Get last vertex coordinates
    virtual void final_vertex_coordinate(Vector<double>& vertex) = 0;

    /// \short Connects the initial vertex of the curve section to a desired
    /// target polyline by specifying the vertex number. There is a checking
    /// which verifies that the initial vertex is close enough to the
    /// destination vertex on the target polyline by no more than the specified
    /// tolerance
    void connect_initial_vertex_to_polyline(
      TriangleMeshPolyLine* polyline_pt,
      const unsigned& vertex_number,
      const double& tolerance_for_connection = 1.0e-14);

    /// \short Connects the final vertex of the curve section to a desired
    /// target polyline by specifying the vertex number. There is a checking
    /// which verifies that the final vertex is close enough to the
    /// destination vertex on the target polyline by no more than the specified
    /// tolerance
    void connect_final_vertex_to_polyline(
      TriangleMeshPolyLine* polyline_pt,
      const unsigned& vertex_number,
      const double& tolerance_for_connection = 1.0e-14);

    /// \short Connects the initial vertex of the curve section to a desired
    /// target curviline by specifying the s value (intrinsic value on the
    /// geometric object of the curviline) where to connect on the target
    /// curviline. There is a checking which verifies that the initial vertex
    /// and the coordinates on the given s value are close enough by no more
    /// than the given tolerance
    void connect_initial_vertex_to_curviline(
      TriangleMeshCurviLine* curviline_pt,
      const double& s_value,
      const double& tolerance_for_connection = 1.0e-14);

    /// \short Connects the final vertex of the curve section to a desired
    /// target curviline by specifying the s value (intrinsic value on the
    /// geometric object of the curviline) where to connect on the target
    /// curviline. There is a checking which verifies that the final vertex
    /// and the coordinates on the given s value are close enough by no more
    /// than the given tolerance
    void connect_final_vertex_to_curviline(
      TriangleMeshCurviLine* curviline_pt,
      const double& s_value,
      const double& tolerance_for_connection = 1.0e-14);

    /// Test whether initial vertex is connected or not
    bool is_initial_vertex_connected() const
    {
      return Initial_vertex_connected;
    }

    /// Sets the initial vertex as connected
    void set_initial_vertex_connected()
    {
      Initial_vertex_connected = true;
    }

    /// Sets the initial vertex as non connected
    void unset_initial_vertex_connected()
    {
      Initial_vertex_connected = false;
    }

    /// Set the initial vertex connection as suspended, it will be
    /// resumed when the method to resume the connections is called
    /// This method is only used in a distributed context, when the
    /// boundary to connect is no longer part of the domain in the
    /// processor
    void suspend_initial_vertex_connected()
    {
      if (Initial_vertex_connected)
      {
        Initial_vertex_connected = false;
        Initial_vertex_connected_suspended = true;
      }
    }

    /// Resumes the initial vertex connection, it may be that after load
    /// balancing the boundary to which the connection was suspended be part
    /// of the domain
    void resume_initial_vertex_connected()
    {
      if (Initial_vertex_connected_suspended)
      {
        Initial_vertex_connected = true;
        Initial_vertex_connected_suspended = false;
      }
    }

    /// Test whether final vertex is connected or not
    bool is_final_vertex_connected() const
    {
      return Final_vertex_connected;
    }

    /// Sets the final vertex as connected
    void set_final_vertex_connected()
    {
      Final_vertex_connected = true;
    }

    /// Sets the final vertex as non connected
    void unset_final_vertex_connected()
    {
      Final_vertex_connected = false;
    }

    /// Set the final vertex connection as suspended, it will be
    /// resumed when the method to resume the connections is called
    /// This method is only used in a distributed context, when the
    /// boundary to connect is no longer part of the domain in the
    /// processor
    void suspend_final_vertex_connected()
    {
      if (Final_vertex_connected)
      {
        Final_vertex_connected = false;
        Final_vertex_connected_suspended = true;
      }
    }

    /// Resumes the final vertex connection, it may be that after load
    /// balancing the boundary to which the connection was suspended be part
    /// of the domain
    void resume_final_vertex_connected()
    {
      if (Final_vertex_connected_suspended)
      {
        Final_vertex_connected = true;
        Final_vertex_connected_suspended = false;
      }
    }

    /// Gets the id to which the initial end is connected
    unsigned initial_vertex_connected_bnd_id() const
    {
      return Initial_vertex_connected_bnd_id;
    }

    /// Sets the id to which the initial end is connected
    unsigned& initial_vertex_connected_bnd_id()
    {
      return Initial_vertex_connected_bnd_id;
    }

    /// Gets the vertex number to which the initial end is connected
    unsigned initial_vertex_connected_n_vertex() const
    {
      return Initial_vertex_connected_n_vertex;
    }

    /// Sets the vertex number to which the initial end is connected
    unsigned& initial_vertex_connected_n_vertex()
    {
      return Initial_vertex_connected_n_vertex;
    }

    /// Gets the boundary chunk to which the initial end is connected
    unsigned initial_vertex_connected_n_chunk() const
    {
      return Initial_vertex_connected_n_chunk;
    }

    /// Sets the boundary chunk to which the initial end is connected
    unsigned& initial_vertex_connected_n_chunk()
    {
      return Initial_vertex_connected_n_chunk;
    }

    /// Gets the id to which the final end is connected
    unsigned final_vertex_connected_bnd_id() const
    {
      return Final_vertex_connected_bnd_id;
    }

    /// Sets the id to which the final end is connected
    unsigned& final_vertex_connected_bnd_id()
    {
      return Final_vertex_connected_bnd_id;
    }

    /// Sets the vertex number to which the final end is connected
    unsigned final_vertex_connected_n_vertex() const
    {
      return Final_vertex_connected_n_vertex;
    }

    /// Gets the vertex number to which the final end is connected
    unsigned& final_vertex_connected_n_vertex()
    {
      return Final_vertex_connected_n_vertex;
    }

    /// Gets the boundary chunk to which the final end is connected
    unsigned final_vertex_connected_n_chunk() const
    {
      return Final_vertex_connected_n_chunk;
    }

    /// Sets the boundary chunk to which the final end is connected
    unsigned& final_vertex_connected_n_chunk()
    {
      return Final_vertex_connected_n_chunk;
    }

    /// Test whether the initial vertex is connected to a curviline
    bool is_initial_vertex_connected_to_curviline() const
    {
      return Initial_vertex_connected_to_curviline;
    }

    /// Sets the initial vertex as connected to a curviline
    void set_initial_vertex_connected_to_curviline()
    {
      Initial_vertex_connected_to_curviline = true;
    }

    /// Sets the initial vertex as non connected to a curviline
    void unset_initial_vertex_connected_to_curviline()
    {
      Initial_vertex_connected_to_curviline = false;
    }

    /// Test whether the final vertex is connected to a curviline
    bool is_final_vertex_connected_to_curviline() const
    {
      return Final_vertex_connected_to_curviline;
    }

    /// Sets the final vertex as connected to a curviline
    void set_final_vertex_connected_to_curviline()
    {
      Final_vertex_connected_to_curviline = true;
    }

    /// Sets the final vertex as non connected to a curviline
    void unset_final_vertex_connected_to_curviline()
    {
      Final_vertex_connected_to_curviline = false;
    }

    /// \short Gets the s value to which the initial end is connected
    double initial_s_connection_value() const
    {
      return Initial_s_connection_value;
    }

    /// \short Sets the s value to which the initial end is connected
    double& initial_s_connection_value()
    {
      return Initial_s_connection_value;
    }

    /// \short Gets the s value to which the final end is connected
    double final_s_connection_value() const
    {
      return Final_s_connection_value;
    }

    /// \short Sets the s value to which the final end is connected
    double& final_s_connection_value()
    {
      return Final_s_connection_value;
    }

    /// \short Gets the tolerance value for connections among
    /// curvilines
    double tolerance_for_s_connection() const
    {
      return Tolerance_for_s_connection;
    }

    /// \short Sets the tolerance value for connections among
    /// curvilines
    double& tolerance_for_s_connection()
    {
      return Tolerance_for_s_connection;
    }

  protected:
    /// \short Used for stating if the initial end is connected
    /// to another boundary
    bool Initial_vertex_connected;

    /// \short Used for stating if the final end is connected
    /// to another boundary
    bool Final_vertex_connected;

    /// \short Indicates if the connection is suspended because the
    /// boundary to connect is no longer part of the domain (only used in
    /// a distributed context)
    bool Initial_vertex_connected_suspended;

    /// \short Indicates if the connection is suspended because the
    /// boundary to connect is no longer part of the domain (only used in
    /// a distributed context)
    bool Final_vertex_connected_suspended;

    /// Stores the id to which the initial end is connected
    unsigned Initial_vertex_connected_bnd_id;

    /// \short Stores the vertex number used for connection with
    /// the initial end
    unsigned Initial_vertex_connected_n_vertex;

    /// \short Stores the chunk number of the boundary to which is
    /// connected th initial end
    unsigned Initial_vertex_connected_n_chunk;

    /// Stores the id to which the initial end is connected
    unsigned Final_vertex_connected_bnd_id;

    /// \short Stores the vertex number used for connection with
    /// the final end
    unsigned Final_vertex_connected_n_vertex;

    /// \short Stores the chunk number of the boundary to which is
    /// connected th initial end
    unsigned Final_vertex_connected_n_chunk;

    /// States if the initial vertex is connected to a curviline
    bool Initial_vertex_connected_to_curviline;

    /// States if the final vertex is connected to a curviline
    bool Final_vertex_connected_to_curviline;

    /// \short Stores the s value used for connecting the
    /// initial end with a curviline
    double Initial_s_connection_value;

    /// \short Stores the s value used for connecting the
    /// final end with a curviline
    double Final_s_connection_value;

    /// Tolerance used for connecting the ends to a curviline
    double Tolerance_for_s_connection;

  private:
    /// Tolerance for refinement of curve sections (neg if refinement is
    /// disabled)
    double Refinement_tolerance;

    /// Tolerance for unrefinement of curve sections (neg if refinement
    /// is disabled)
    double Unrefinement_tolerance;

    /// Maximum allowed distance between two vertices on the polyline
    /// representation of the curve section
    double Maximum_length;
  };


  //=====================================================================
  /// Class definining a curvilinear triangle mesh boundary in terms
  /// of a GeomObject. Curvlinear equivalent of PolyLine.
  //=====================================================================
  class TriangleMeshCurviLine : public virtual TriangleMeshCurveSection
  {
  public:
    /// \short Constructor: Specify GeomObject, the start and end coordinates
    /// of the relevant boundary in terms of the GeomObject's intrinsic
    /// coordinate, the number of (initially straight-line) segments that
    /// this GeomObject is to be split up into, and the boundary ID.
    /// The final optional boolean argument specifies if vertices in
    /// polygonhal represenation are spaced
    /// (approximately) evenly in arclength along the GeomObject [true,
    /// default] or in equal increments in zeta.
    /// This is the curvlinear equivalent of PolyLine.
    TriangleMeshCurviLine(GeomObject* geom_object_pt,
                          const double& zeta_start,
                          const double& zeta_end,
                          const unsigned& nsegment,
                          const unsigned& boundary_id,
                          const bool& space_vertices_evenly_in_arclength = true,
                          const unsigned& boundary_chunk = 0)
      : TriangleMeshCurveSection(),
        Geom_object_pt(geom_object_pt),
        Zeta_start(zeta_start),
        Zeta_end(zeta_end),
        Nsegment(nsegment),
        Boundary_id(boundary_id),
        Space_vertices_evenly_in_arclength(space_vertices_evenly_in_arclength),
        Boundary_chunk(boundary_chunk)
    {
    }


    /// \short Empty Destuctor
    virtual ~TriangleMeshCurviLine() {}

    /// Pointer to GeomObject that represents this part of the boundary
    GeomObject* geom_object_pt()
    {
      return Geom_object_pt;
    }

    /// Start coordinate in terms of the GeomObject's intrinsic coordinate
    double zeta_start()
    {
      return Zeta_start;
    }

    /// End coordinate in terms of the GeomObject's intrinsic coordinate
    double zeta_end()
    {
      return Zeta_end;
    }

    /// \short Number of (initially straight-line) segments that this part of
    /// the boundary is to be represented by
    unsigned nsegment() const
    {
      return Nsegment;
    }

    /// \short Number of (initially straight-line) segments that this part of
    /// the boundary is to be represented by. This version allows the change of
    /// the number of segments
    unsigned& nsegment()
    {
      return Nsegment;
    }

    /// Boundary ID
    unsigned boundary_id() const
    {
      return Boundary_id;
    }

    /// Boundary chunk (Used when a boundary is represented by more
    /// than one polyline
    unsigned boundary_chunk() const
    {
      return Boundary_chunk;
    }

    /// Output curvilinear boundary at n_sample (default: 50) points
    void output(std::ostream& outfile, const unsigned& n_sample = 50)
    {
      outfile << "ZONE T=\"Boundary " << Boundary_id << "\"\n";
      Vector<double> zeta(1);
      Vector<double> r(2);
      for (unsigned i = 0; i < n_sample; i++)
      {
        zeta[0] = Zeta_start +
                  (Zeta_end - Zeta_start) * double(i) / double(n_sample - 1);
        Geom_object_pt->position(zeta, r);
        outfile << r[0] << " " << r[1] << std::endl;
      }
    }

    /// \short Boolean to indicate if vertices in polygonal representation
    /// of the Curvline are spaced (approximately) evenly in arclength
    /// along the GeomObject [true] or in equal increments in zeta [false]
    bool space_vertices_evenly_in_arclength() const
    {
      return Space_vertices_evenly_in_arclength;
    }

    /// Number of vertices
    unsigned nvertex() const
    {
      return 2;
    }

    /// Get first vertex coordinates
    void initial_vertex_coordinate(Vector<double>& vertex)
    {
      Vector<double> z(1);
      z[0] = Zeta_start;
      Geom_object_pt->position(z, vertex);
    }

    /// Get last vertex coordinates
    void final_vertex_coordinate(Vector<double>& vertex)
    {
      Vector<double> z(1);
      z[0] = Zeta_end;
      Geom_object_pt->position(z, vertex);
    }

    /// \short Does the vector for storing connections has elements?
    bool are_there_connection_points()
    {
      return !Connection_points_pt.empty();
    }

    /// \short Returns the connection points vector
    Vector<double>* connection_points_pt()
    {
      return &Connection_points_pt;
    }

    /// Add the connection point (z_value) to the connection
    /// points that receive the curviline
    void add_connection_point(const double& z_value,
                              const double& tol = 1.0e-12)
    {
      // If we are trying to connect to the initial or final
      // point then it is not necessary to add this point
      // to the list since it will explicitly be created when
      // transforming the curviline to polyline
      if (std::fabs(z_value - Zeta_start) < tol ||
          std::fabs(z_value - Zeta_end) < tol)
      {
        return;
      }

      // We need to deal with repeated connection points,
      // if the connection point is already in the list then
      // we will not add it!!!
      // Search for repeated elements
      unsigned n_size = Connection_points_pt.size();
      for (unsigned i = 0; i < n_size; i++)
      {
        if (fabs(z_value - Connection_points_pt[i]) < tol)
        {
          return;
        }
      }

      // Only add the connection point if it is not the
      // initial or final zeta value and it is not already on the
      // list
      Connection_points_pt.push_back(z_value);
    }

  private:
    /// Pointer to GeomObject that represents this part of the boundary
    GeomObject* Geom_object_pt;

    /// Start coordinate in terms of the GeomObject's intrinsic coordinate
    double Zeta_start;

    /// End coordinate in terms of the GeomObject's intrinsic coordinate
    double Zeta_end;

    /// Number of (initially straight-line) segments that this part of the
    /// boundary is to be represented by
    unsigned Nsegment;

    /// Boundary ID
    unsigned Boundary_id;

    /// \short Boolean to indicate if vertices in polygonal representation
    /// of the Curviline are spaced (approximately) evenly in arclength
    /// along the GeomObject [true] or in equal increments in zeta [false]
    bool Space_vertices_evenly_in_arclength;

    /// Boundary chunk (Used when a boundary is represented by more
    /// than one polyline
    unsigned Boundary_chunk;

    /// \short Stores the information for connections received on the
    /// curviline. Used when converting to polyline
    Vector<double> Connection_points_pt;
  };


  //=====================================================================
  /// Class defining a polyline for use in Triangle Mesh generation
  //=====================================================================
  class TriangleMeshPolyLine : public virtual TriangleMeshCurveSection
  {
  public:
    /// \short Constructor: Takes vectors of vertex coordinates in order
    /// Also allows the optional specification of a boundary ID -- useful
    /// in a mesh generation context. If not specified it defaults to zero.
    TriangleMeshPolyLine(const Vector<Vector<double>>& vertex_coordinate,
                         const unsigned& boundary_id,
                         const unsigned& boundary_chunk = 0)
      : TriangleMeshCurveSection(),
        Vertex_coordinate(vertex_coordinate),
        Boundary_id(boundary_id),
        Boundary_chunk(boundary_chunk)
    {
#ifdef PARANOID
      unsigned nvert = Vertex_coordinate.size();
      for (unsigned i = 0; i < nvert; i++)
      {
        if (Vertex_coordinate[i].size() != 2)
        {
          std::ostringstream error_stream;
          error_stream << "TriangleMeshPolyLine should only be used in 2D!\n"
                       << "Your Vector of coordinates, contains data for "
                       << Vertex_coordinate[i].size()
                       << "-dimensional coordinates." << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
#endif
    }

    /// Empty destructor
    virtual ~TriangleMeshPolyLine() {}

    /// Number of vertices
    unsigned nvertex() const
    {
      return Vertex_coordinate.size();
    }

    /// Number of segments
    unsigned nsegment() const
    {
      return Vertex_coordinate.size() - 1;
    }

    /// Boundary id
    unsigned boundary_id() const
    {
      return Boundary_id;
    }

    /// Boundary chunk (Used when a boundary is represented by more
    /// than one polyline
    unsigned boundary_chunk() const
    {
      return Boundary_chunk;
    }

    /// Coordinate vector of i-th vertex (const version)
    Vector<double> vertex_coordinate(const unsigned& i) const
    {
      return Vertex_coordinate[i];
    }

    /// Coordinate vector of i-th vertex
    Vector<double>& vertex_coordinate(const unsigned& i)
    {
      return Vertex_coordinate[i];
    }

    /// Get first vertex coordinates
    void initial_vertex_coordinate(Vector<double>& vertex)
    {
      vertex = Vertex_coordinate[0];
    }

    /// Get last vertex coordinates
    void final_vertex_coordinate(Vector<double>& vertex)
    {
      vertex = Vertex_coordinate[nvertex() - 1];
    }

    /// Output the polyline -- n_sample is ignored
    void output(std::ostream& outfile, const unsigned& n_sample = 50)
    {
      outfile << "ZONE T=\"TriangleMeshPolyLine with boundary ID" << Boundary_id
              << "\"" << std::endl;
      unsigned nvert = Vertex_coordinate.size();
      for (unsigned i = 0; i < nvert; i++)
      {
        outfile << Vertex_coordinate[i][0] << " " << Vertex_coordinate[i][1]
                << std::endl;
      }
    }

    /// Reverse the polyline, this includes the connection information
    /// and the vertices order
    void reverse()
    {
      // Do the reversing of the connection information

      // Is there a connection to the initial vertex
      const bool initial_connection = is_initial_vertex_connected();

      // Is there a connection to the initial vertex
      const bool final_connection = is_final_vertex_connected();

      // If there are any connection at the ends that info. needs to be
      // reversed
      if (initial_connection || final_connection)
      {
        // Backup the connection information

        // -------------------------------------------------------------------
        // Backup the initial vertex connection information
        // The boundary id
        const unsigned backup_initial_vertex_connected_bnd_id =
          initial_vertex_connected_bnd_id();
        // The chunk number
        const unsigned backup_initial_vertex_connected_chunk =
          initial_vertex_connected_n_chunk();
        // The vertex number
        const unsigned backup_initial_vertex_connected_n_vertex =
          initial_vertex_connected_n_vertex();
        // Is it connected to a curviline
        const bool backup_initial_vertex_connected_to_curviline =
          is_initial_vertex_connected_to_curviline();
        // The s value for the curviline connection
        const double backup_initial_s_connection = initial_s_connection_value();
        // The tolerance
        const double backup_initial_s_tolerance = tolerance_for_s_connection();

        // -------------------------------------------------------------------
        // Backup the final vertex connection information
        // The boundary id
        const unsigned backup_final_vertex_connected_bnd_id =
          final_vertex_connected_bnd_id();
        // The chunk number
        const unsigned backup_final_vertex_connected_chunk =
          final_vertex_connected_n_chunk();
        // The vertex number
        const unsigned backup_final_vertex_connected_n_vertex =
          final_vertex_connected_n_vertex();
        // Is it connected to a curviline
        const bool backup_final_vertex_connected_to_curviline =
          is_final_vertex_connected_to_curviline();
        // The s value for the curviline connection
        const double backup_final_s_connection = final_s_connection_value();
        // The tolerance
        const double backup_final_s_tolerance = tolerance_for_s_connection();
        // -------------------------------------------------------------------

        // Disconnect the polyline

        // Disconnect the initial vertex
        unset_initial_vertex_connected();
        unset_initial_vertex_connected_to_curviline();

        // Disconnect the final vertex
        unset_final_vertex_connected();
        unset_final_vertex_connected_to_curviline();

        // Now reconnected but in inverted order
        if (initial_connection)
        {
          // Set the final vertex as connected
          set_final_vertex_connected();
          // Copy the boundary id
          final_vertex_connected_bnd_id() =
            backup_initial_vertex_connected_bnd_id;
          // Copy the chunk number
          final_vertex_connected_n_chunk() =
            backup_initial_vertex_connected_chunk;
          // Copy the vertex number
          final_vertex_connected_n_vertex() =
            backup_initial_vertex_connected_n_vertex;
          // Is it connected to a curviline
          if (backup_initial_vertex_connected_to_curviline)
          {
            // Set the final vertex as connected to curviline
            set_final_vertex_connected_to_curviline();
            // Copy the s value to connected
            final_s_connection_value() = backup_initial_s_connection;
            // Copy the tolerance
            tolerance_for_s_connection() = backup_initial_s_tolerance;
          } // if (backup_initial_vertex_connected_to_curviline)

        } // if (initial_connection)

        if (final_connection)
        {
          // Set the initial vertex as connected
          set_initial_vertex_connected();
          // Copy the boundary id
          initial_vertex_connected_bnd_id() =
            backup_final_vertex_connected_bnd_id;
          // Copy the chunk number
          initial_vertex_connected_n_chunk() =
            backup_final_vertex_connected_chunk;
          // Copy the vertex number
          initial_vertex_connected_n_vertex() =
            backup_final_vertex_connected_n_vertex;
          // Is it connected to a curviline
          if (backup_final_vertex_connected_to_curviline)
          {
            // Set the initial vertex as connected to curviline
            set_initial_vertex_connected_to_curviline();
            // Copy the s value to connected
            initial_s_connection_value() = backup_final_s_connection;
            // Copy the tolerance
            tolerance_for_s_connection() = backup_final_s_tolerance;
          } // if (backup_final_vertex_connected_to_curviline)

        } // if (final_connection)

      } // if (initial_connection || final_connection)

      // Do the reversing of the vertices
      std::reverse(Vertex_coordinate.begin(), Vertex_coordinate.end());
    }

  private:
    /// Vector of Vector of vertex coordinates
    Vector<Vector<double>> Vertex_coordinate;

    /// Boundary ID
    unsigned Boundary_id;

    /// Boundary chunk (Used when a boundary is represented by more
    /// than one polyline
    unsigned Boundary_chunk;
  };


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //===================================================================
  /// \short Namespace that allows the specification of a tolerance
  /// between vertices at the ends of polylines that are supposed
  /// to be at the same position.
  //===================================================================
  namespace ToleranceForVertexMismatchInPolygons
  {
    /// \short Acceptable discrepancy for mismatch in vertex coordinates.
    /// In paranoid mode, the code will die if the beginning/end of
    /// two adjacent polylines differ by more than that. If the
    /// discrepancy is smaller (but nonzero) one of the vertex coordinates
    /// get adjusted to match perfectly; without paranoia the vertex
    /// coordinates are taken as they come...
    extern double Tolerable_error;
  } // namespace ToleranceForVertexMismatchInPolygons

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //=====================================================================
  // \short Class defining triangle mesh curves. Abstract class for
  /// closed curves and open curves. All TriangleMeshCurves are composed
  /// of a Vector of TriangleMeshCurveSections
  //=====================================================================
  class TriangleMeshCurve
  {
  public:
    /// Empty constructor
    TriangleMeshCurve(const Vector<TriangleMeshCurveSection*>& curve_section_pt)
      : Curve_section_pt(curve_section_pt),
        Polyline_refinement_tolerance(0.08),
        Polyline_unrefinement_tolerance(0.04)
    {
    }

    /// Empty destructor
    virtual ~TriangleMeshCurve() {}

    /// Number of vertices
    virtual unsigned nvertices() const = 0;

    /// Total number of segments
    virtual unsigned nsegments() const = 0;

    /// Return max boundary id of associated curves
    unsigned max_boundary_id()
    {
      unsigned max = 0;
      unsigned n_curve_section = ncurve_section();
      for (unsigned i = 0; i < n_curve_section; i++)
      {
        unsigned boundary_id = Curve_section_pt[i]->boundary_id();
        if (boundary_id > max)
        {
          max = boundary_id;
        }
      }
      return max;
    }

    /// Number of constituent curves
    virtual unsigned ncurve_section() const
    {
      return Curve_section_pt.size();
    }

    /// \short Enable refinement of polylines to create a better
    /// representation of curvilinear boundaries (e.g. in free-surface
    /// problems). See tutorial for
    /// interpretation of the optional argument which specifies the
    /// refinement tolerance. It defaults to 0.08 and the smaller the
    /// number the finer the surface representation.
    void enable_polyline_refinement(const double& tolerance = 0.08)
    {
      Polyline_refinement_tolerance = tolerance;
      // Establish the refinement tolerance for all the
      // curve sections on the TriangleMeshCurve
      unsigned n_curve_sections = Curve_section_pt.size();
      for (unsigned i = 0; i < n_curve_sections; i++)
      {
        Curve_section_pt[i]->set_refinement_tolerance(
          Polyline_refinement_tolerance);
      }
    }

    /// \short Set tolerance for refinement of polylines to create a better
    /// representation of curvilinear boundaries (e.g. in free-surface
    /// problems). See tutorial for
    /// interpretation of the refinement tolerance. (The smaller the
    /// number the finer the surface representation). If set to
    /// a negative value, we're switching off refinement --
    /// equivalent to calling disable_polyline_refinement()
    void set_polyline_refinement_tolerance(const double& tolerance)
    {
      Polyline_refinement_tolerance = tolerance;
      // Establish the refinement tolerance for all the
      // curve sections on the TriangleMeshCurve
      unsigned n_curve_sections = Curve_section_pt.size();
      for (unsigned i = 0; i < n_curve_sections; i++)
      {
        Curve_section_pt[i]->set_refinement_tolerance(
          Polyline_refinement_tolerance);
      }
    }

    /// \short Get tolerance for refinement of polylines to create a better
    /// representation of curvilinear boundaries (e.g. in free-surface
    /// problems). See tutorial for
    /// interpretation. If it's negative refinement is disabled.
    double polyline_refinement_tolerance()
    {
      return Polyline_refinement_tolerance;
    }

    /// \short Disable refinement of polylines
    void disable_polyline_refinement()
    {
      Polyline_refinement_tolerance = -1.0;
      // Disable the refinement tolerance for all the
      // curve sections on the TriangleMeshCurve
      unsigned n_curve_sections = Curve_section_pt.size();
      for (unsigned i = 0; i < n_curve_sections; i++)
      {
        Curve_section_pt[i]->disable_refinement_tolerance();
      }
    }

    /// \short Enable unrefinement of polylines to avoid unnecessarily large
    /// numbers of elements on of curvilinear boundaries (e.g. in free-surface
    /// problems). See tutorial for
    /// interpretation of the optional argument which specifies the
    /// unrefinement tolerance. It defaults to 0.04 and the larger the number
    /// the more agressive we are when removing unnecessary vertices on
    /// gently curved polylines.
    void enable_polyline_unrefinement(const double& tolerance = 0.04)
    {
      Polyline_unrefinement_tolerance = tolerance;
      // Establish the unrefinement tolerance for all the
      // curve sections on the TriangleMeshCurve
      unsigned n_curve_sections = Curve_section_pt.size();
      for (unsigned i = 0; i < n_curve_sections; i++)
      {
        Curve_section_pt[i]->set_unrefinement_tolerance(
          Polyline_unrefinement_tolerance);
      }
    }

    /// \short Set tolerance for unrefinement of polylines
    /// to avoid unnecessarily large
    /// numbers of elements on of curvilinear boundaries (e.g. in free-surface
    /// problems). See tutorial for
    /// interpretation of the optional argument which specifies the
    /// unrefinement tolerance. It defaults to 0.04 and the larger the number
    /// the more agressive we are when removing unnecessary vertices on
    /// gently curved polylines. If set to
    /// a negative value, we're switching off unrefinement --
    /// equivalent to calling disable_polyline_unrefinement()
    void set_polyline_unrefinement_tolerance(const double& tolerance)
    {
      Polyline_unrefinement_tolerance = tolerance;
      // Establish the unrefinement tolerance for all the
      // curve sections on the TriangleMeshCurve
      unsigned n_curve_sections = Curve_section_pt.size();
      for (unsigned i = 0; i < n_curve_sections; i++)
      {
        Curve_section_pt[i]->set_unrefinement_tolerance(
          Polyline_unrefinement_tolerance);
      }
    }

    /// \short Get tolerance for unrefinement of polylines to create a better
    /// representation of curvilinear boundaries (e.g. in free-surface
    /// problems). See tutorial for
    /// interpretation. If it's negative unrefinement is disabled.
    double polyline_unrefinement_tolerance()
    {
      return Polyline_unrefinement_tolerance;
    }

    /// \short Disable unrefinement of polylines
    void disable_polyline_unrefinement()
    {
      Polyline_unrefinement_tolerance = -1.0;
      // Disable the unrefinement tolerance for all the
      // curve sections on the TriangleMeshCurve
      unsigned n_curve_sections = Curve_section_pt.size();
      for (unsigned i = 0; i < n_curve_sections; i++)
      {
        Curve_section_pt[i]->disable_unrefinement_tolerance();
      }
    }

    /// Output each sub-boundary at n_sample (default: 50) points
    virtual void output(std::ostream& outfile,
                        const unsigned& n_sample = 50) = 0;

    /// Pointer to i-th constituent curve section
    virtual TriangleMeshCurveSection* curve_section_pt(const unsigned& i) const
    {
      return Curve_section_pt[i];
    }

    /// Pointer to i-th constituent curve section
    virtual TriangleMeshCurveSection*& curve_section_pt(const unsigned& i)
    {
      return Curve_section_pt[i];
    }

  protected:
    /// Vector of curve sections
    Vector<TriangleMeshCurveSection*> Curve_section_pt;

  private:
    /// Tolerance for refinement of polylines (neg if refinement is disabled)
    double Polyline_refinement_tolerance;

    /// Tolerance for unrefinement of polylines (neg if refinement is disabled)
    double Polyline_unrefinement_tolerance;
  };

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  //=====================================================================
  /// Base class defining a closed curve for the Triangle mesh generation
  //=====================================================================
  class TriangleMeshClosedCurve : public virtual TriangleMeshCurve
  {
  public:
    /// Constructor prototype
    TriangleMeshClosedCurve(
      const Vector<TriangleMeshCurveSection*>& curve_section_pt,
      const Vector<double>& internal_point_pt = Vector<double>(0),
      const bool& is_internal_point_fixed = false);

    /// Empty destructor
    virtual ~TriangleMeshClosedCurve() {}

    /// Number of vertices
    unsigned nvertices() const
    {
      unsigned n_curve_section = ncurve_section();
      unsigned n_vertices = 0;
      for (unsigned j = 0; j < n_curve_section; j++)
      {
        // Storing the number of the vertices
        n_vertices += Curve_section_pt[j]->nvertex() - 1;
      }
      // If there's just one boundary. All the vertices should be counted
      if (n_curve_section == 1)
      {
        n_vertices += 1;
      }
      return n_vertices;
    }

    /// Total number of segments
    unsigned nsegments() const
    {
      unsigned n_curve_section = ncurve_section();
      unsigned nseg = 0;
      for (unsigned j = 0; j < n_curve_section; j++)
      {
        nseg += Curve_section_pt[j]->nsegment();
      }
      // If there's just one boundary poly line we have another segment
      // connecting the last vertex to the first one
      if (n_curve_section == 1)
      {
        nseg += 1;
      }
      return nseg;
    }

    /// Output each sub-boundary at n_sample (default: 50) points
    void output(std::ostream& outfile, const unsigned& n_sample = 50)
    {
      unsigned nb = Curve_section_pt.size();
      for (unsigned i = 0; i < nb; i++)
      {
        Curve_section_pt[i]->output(outfile, n_sample);
      }

      if (!Internal_point_pt.empty())
      {
        outfile << "ZONE T=\"Internal point\"\n";
        outfile << Internal_point_pt[0] << " " << Internal_point_pt[1] << "\n";
      }
    }

    /// Coordinates of the internal point
    Vector<double> internal_point() const
    {
      return Internal_point_pt;
    }

    /// Coordinates of the internal point
    Vector<double>& internal_point()
    {
      return Internal_point_pt;
    }

    /// Fix the internal point (i.e. do not allow our automatic machinery
    /// to update it)
    void fix_internal_point()
    {
      Is_internal_point_fixed = true;
    }

    /// Unfix the internal point (i.e. allow our automatic machinery
    /// to update it)
    void unfix_internal_point()
    {
      Is_internal_point_fixed = false;
    }

    /// Test whether the internal point is fixed
    bool is_internal_point_fixed() const
    {
      return Is_internal_point_fixed;
    }

  protected:
    /// Vector of vertex coordinates
    Vector<double> Internal_point_pt;

    /// Indicate whether the internal point should be updated automatically
    bool Is_internal_point_fixed;
  };

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //=====================================================================
  /// Class defining a closed polygon for the Triangle mesh generation
  //=====================================================================
  class TriangleMeshPolygon : public virtual TriangleMeshClosedCurve
  {
  public:
    /// \short Constructor: Specify vector of pointers to
    /// TriangleMeshCurveSection that define the boundary of the segments of the
    /// polygon. Each TriangleMeshCurveSection has its own boundary ID and can
    /// contain multiple (straight-line) segments. For consistency across the
    /// various uses of this class, we insist that the closed boundary
    /// is represented by at least two separate TriangleMeshCurveSection
    /// whose joint vertices must be specified in both.
    /// (This is to allow the setup of unique boundary coordinate(s)
    /// around the polygon.) This may seem slightly annoying
    /// in cases where a polygon really only represents a single
    /// boundary, but...
    /// Note: The specified vector of pointers must consist of only
    /// TriangleMeshPolyLine elements. There is a checking on the PARANOID
    /// mode for this constraint
    TriangleMeshPolygon(
      const Vector<TriangleMeshCurveSection*>& boundary_polyline_pt,
      const Vector<double>& internal_point_pt = Vector<double>(0),
      const bool& is_internal_point_fixed = false);

    /// Empty virtual destructor
    virtual ~TriangleMeshPolygon() {}

    /// Number of constituent curves
    unsigned ncurve_section() const
    {
      return npolyline();
    }

    /// Number of constituent polylines
    unsigned npolyline() const
    {
      return Curve_section_pt.size();
    }

    /// Pointer to i-th constituent polyline
    TriangleMeshPolyLine* polyline_pt(const unsigned& i) const
    {
      TriangleMeshPolyLine* tmp_polyline =
        dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[i]);
#ifdef PARANOID
      if (tmp_polyline == NULL)
      {
        std::ostringstream error_stream;
        error_stream
          << "The (" << i << ") TriangleMeshCurveSection is not a "
          << "TriangleMeshPolyLine\nThe TriangleMeshPolygon object"
          << "is constituent of only TriangleMeshPolyLine objects.\n"
          << "The problem could be generated when changing the constituent "
          << "objects of the TriangleMeshPolygon.\nCheck where you got "
          << "access to this objects and review that you are not introducing "
          << "any other objects than TriangleMeshPolyLines" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return tmp_polyline;
    }

    /// Pointer to i-th constituent polyline
    TriangleMeshPolyLine* polyline_pt(const unsigned& i)
    {
      TriangleMeshPolyLine* tmp_polyline =
        dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[i]);
#ifdef PARANOID
      if (tmp_polyline == NULL)
      {
        std::ostringstream error_stream;
        error_stream
          << "The (" << i << ") TriangleMeshCurveSection is not a "
          << "TriangleMeshPolyLine\nThe TriangleMeshPolygon object"
          << "is constituent of only TriangleMeshPolyLine objects.\n"
          << "The problem could be generated when changing the constituent "
          << "objects of the TriangleMeshPolygon.\nCheck where you got "
          << "access to this objects and review that you are not introducing "
          << "any other objects than TriangleMeshPolyLines" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return tmp_polyline;
    }

    /// Return vector of boundary ids of associated polylines
    Vector<unsigned> polygon_boundary_id()
    {
      // Get the number of polylines
      unsigned nline = npolyline();
      Vector<unsigned> boundary_id(nline);

      // Loop over the polyline to get the id
      for (unsigned iline = 0; iline < nline; iline++)
      {
        boundary_id[iline] = Curve_section_pt[iline]->boundary_id();
      }
      return boundary_id;
    }

    /// \short Is re-distribution of polyline segments in the curve
    /// between different boundaries during adaptation enabled?
    bool is_redistribution_of_segments_between_polylines_enabled()
    {
      return Enable_redistribution_of_segments_between_polylines;
    }

    /// \short Enable re-distribution of polyline segments in the curve
    /// between different boundaries during adaptation
    void enable_redistribution_of_segments_between_polylines()
    {
      Enable_redistribution_of_segments_between_polylines = true;
    }

    /// \short Disable re-distribution of polyline segments in the curve
    /// between different boundaries during adaptation
    void disable_redistribution_of_segments_between_polylines()
    {
      Enable_redistribution_of_segments_between_polylines = false;
    }

    /// \short Test whether curve can update reference
    bool can_update_reference_configuration() const
    {
      return Can_update_configuration;
    }

    /// \short Virtual function that should be overloaded to update the polygons
    /// reference configuration
    virtual void reset_reference_configuration()
    {
      std::ostringstream error_stream;
      error_stream
        << "Broken Default Called\n"
        << "This function should be overloaded if Can_update_configuration = "
           "true,"
        << "\nwhich indicates that the curve in it polylines parts can update "
           "its "
        << "own position (i.e. it\n"
        << "may be a rigid body. Otherwise the update will be via a FaceMesh \n"
        << "representation of the boundary which is appropriate for \n"
        << "general deforming surfaces\n";

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    /// \short Test whether the polygon is fixed or not
    bool is_fixed() const
    {
      return Polygon_fixed;
    }

    /// \short Set the polygon to be fixed
    void set_fixed()
    {
      Polygon_fixed = true;
    }

    /// \short Set the polygon to be allowed to move (default)
    void set_unfixed()
    {
      Polygon_fixed = false;
    }

  protected:
    /// \short Is re-distribution of polyline segments between different
    /// boundaries during adaptation enabled? (Default: false)
    bool Enable_redistribution_of_segments_between_polylines;

    ///\short Boolean flag to indicate whether the polygon can update its
    /// own reference configuration after it has moved i.e. if it is
    /// upgraded to a rigid body rather than being a free surface (default
    /// false)
    bool Can_update_configuration;

  private:
    ///\short Boolean flag to indicate whether the polygon can move
    /// (default false because if it doesn't move this will just lead to
    ///  wasted work)
    bool Polygon_fixed;
  };

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  //=====================================================================
  /// Base class defining an open curve for the Triangle mesh generation
  /// Basically used to define internal boundaries on the mesh
  //=====================================================================
  class TriangleMeshOpenCurve : public virtual TriangleMeshCurve
  {
  public:
    /// Constructor
    TriangleMeshOpenCurve(
      const Vector<TriangleMeshCurveSection*>& curve_section_pt);

    /// Empty destructor
    virtual ~TriangleMeshOpenCurve() {}

    /// Number of vertices
    unsigned nvertices() const
    {
      unsigned n_vertices = 0;
      unsigned n_curve_section = ncurve_section();
      for (unsigned i = 0; i < n_curve_section; i++)
        n_vertices += Curve_section_pt[i]->nvertex();
      // If there's just one boundary. All the vertices should be counted
      if (n_curve_section == 1)
      {
        n_vertices += 1;
      }
      return n_vertices;
    }

    /// Total number of segments
    unsigned nsegments() const
    {
      unsigned n_curve_section = ncurve_section();
      unsigned nseg = 0;
      for (unsigned j = 0; j < n_curve_section; j++)
      {
        nseg += Curve_section_pt[j]->nsegment();
      }
      return nseg;
    }

    /// Output each sub-boundary at n_sample (default: 50) points
    void output(std::ostream& outfile, const unsigned& n_sample = 50)
    {
      unsigned nb = Curve_section_pt.size();
      for (unsigned i = 0; i < nb; i++)
      {
        Curve_section_pt[i]->output(outfile, n_sample);
      }
    }

    /// Pointer to i-th constituent polyline
    TriangleMeshPolyLine* polyline_pt(const unsigned& i) const
    {
      TriangleMeshPolyLine* tmp_polyline =
        dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[i]);
#ifdef PARANOID
      if (tmp_polyline == NULL)
      {
        std::ostringstream error_stream;
        error_stream << "The (" << i
                     << ")-th TriangleMeshCurveSection is not a "
                     << "TriangleMeshPolyLine.\nPlease make sure that when you"
                     << "first created this object the (" << i << ")-th\n"
                     << "TriangleCurveSection is a TriangleMeshPolyLine."
                     << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return tmp_polyline;
    }

    /// Pointer to i-th constituent polyline
    TriangleMeshPolyLine* polyline_pt(const unsigned& i)
    {
      TriangleMeshPolyLine* tmp_polyline =
        dynamic_cast<TriangleMeshPolyLine*>(Curve_section_pt[i]);
#ifdef PARANOID
      if (tmp_polyline == NULL)
      {
        std::ostringstream error_stream;
        error_stream << "The (" << i
                     << ")-th TriangleMeshCurveSection is not a "
                     << "TriangleMeshPolyLine.\nPlease make sure that when you"
                     << "first created this object the (" << i << ")-th\n"
                     << "TriangleCurveSection is a TriangleMeshPolyLine."
                     << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return tmp_polyline;
    }
  };

  //==============start_of_geometry_helper_functions_class================
  /// Contains functions which define the geometry of the mesh, i.e.
  /// regions, boundaries, etc.
  //======================================================================
  class UnstructuredTwoDMeshGeometryBase : public virtual Mesh
  {
  public:
    /// Public static flag to suppress warning; defaults to false
    static bool Suppress_warning_about_regions_and_boundaries;

    /// \short Empty constructor
    UnstructuredTwoDMeshGeometryBase() {}

    /// Broken copy constructor
    UnstructuredTwoDMeshGeometryBase(
      const UnstructuredTwoDMeshGeometryBase& dummy) = delete;

    /// Broken assignment operator
    void operator=(const UnstructuredTwoDMeshGeometryBase&) = delete;

    /// Empty destructor
    ~UnstructuredTwoDMeshGeometryBase() {}

    /// Return the number of regions specified by attributes
    unsigned nregion()
    {
      return Region_element_pt.size();
    }

    /// Return the number of elements in the i-th region
    unsigned nregion_element(const unsigned& i)
    {
      // Create an iterator to iterate over Region_element_pt
      std::map<unsigned, Vector<FiniteElement*>>::iterator it;

      // Find the entry of Region_element_pt associated with the i-th region
      it = Region_element_pt.find(i);

      // If there is an entry associated with the i-th region
      if (it != Region_element_pt.end())
      {
        return (it->second).size();
      }
      else
      {
        return 0;
      }
    } // End of nregion_element

    /// Return the e-th element in the i-th region
    FiniteElement* region_element_pt(const unsigned& i, const unsigned& e)
    {
      // Create an iterator to iterate over Region_element_pt
      std::map<unsigned, Vector<FiniteElement*>>::iterator it;

      // Find the entry of Region_element_pt associated with the i-th region
      it = Region_element_pt.find(i);

      // If there is an entry associated with the i-th region
      if (it != Region_element_pt.end())
      {
        // Return a pointer to the e-th element in the i-th region
        return (it->second)[e];
      }
      else
      {
        // Create a stringstream object
        std::stringstream error_message;

        // Create the error message
        error_message << "There are no regions associated with "
                      << "region ID " << i << ".";

        // Throw an error
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    } // End of region_element_pt

    /// Return the number of attributes used in the mesh
    unsigned nregion_attribute()
    {
      return Region_attribute.size();
    }

    /// Return the attribute associated with region i
    double region_attribute(const unsigned& i)
    {
      return Region_attribute[i];
    }

    /// \short Return the geometric object associated with the b-th boundary or
    /// null if the boundary has associated geometric object.
    GeomObject* boundary_geom_object_pt(const unsigned& b)
    {
      std::map<unsigned, GeomObject*>::iterator it;
      it = Boundary_geom_object_pt.find(b);
      if (it == Boundary_geom_object_pt.end())
      {
        return 0;
      }
      else
      {
        return (*it).second;
      }
    }

    /// \short Return direct access to the geometric object storage
    std::map<unsigned, GeomObject*>& boundary_geom_object_pt()
    {
      return Boundary_geom_object_pt;
    }

    /// \short Return access to the vector of boundary coordinates associated
    /// with each geometric object
    std::map<unsigned, Vector<double>>& boundary_coordinate_limits()
    {
      return Boundary_coordinate_limits;
    }

    /// \short Return access to the coordinate limits associated with
    /// the geometric object associated with boundary b
    Vector<double>& boundary_coordinate_limits(const unsigned& b)
    {
      std::map<unsigned, Vector<double>>::iterator it;
      it = Boundary_coordinate_limits.find(b);
      if (it == Boundary_coordinate_limits.end())
      {
        throw OomphLibError(
          "No coordinate limits associated with this boundary\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        return (*it).second;
      }
    }

    /// Return the number of elements adjacent to boundary b in region r
    inline unsigned nboundary_element_in_region(const unsigned& b,
                                                const unsigned& r) const
    {
      // Need to use a constant iterator here to keep the function "const".
      // Return an iterator to the appropriate entry, if we find it
      std::map<unsigned, Vector<FiniteElement*>>::const_iterator it =
        Boundary_region_element_pt[b].find(r);
      if (it != Boundary_region_element_pt[b].end())
      {
        return (it->second).size();
      }
      // Otherwise there are no elements adjacent to boundary b in the region r
      else
      {
        return 0;
      }
    }

    /// Return pointer to the e-th element adjacent to boundary b in region r
    FiniteElement* boundary_element_in_region_pt(const unsigned& b,
                                                 const unsigned& r,
                                                 const unsigned& e) const
    {
      // Use a constant iterator here to keep function "const" overall
      std::map<unsigned, Vector<FiniteElement*>>::const_iterator it =
        Boundary_region_element_pt[b].find(r);
      if (it != Boundary_region_element_pt[b].end())
      {
        return (it->second)[e];
      }
      else
      {
        return 0;
      }
    }

    /// Return face index of the e-th element adjacent to boundary b in region r
    int face_index_at_boundary_in_region(const unsigned& b,
                                         const unsigned& r,
                                         const unsigned& e) const
    {
      // Use a constant iterator here to keep function "const" overall
      std::map<unsigned, Vector<int>>::const_iterator it =
        Face_index_region_at_boundary[b].find(r);
      if (it != Face_index_region_at_boundary[b].end())
      {
        return (it->second)[e];
      }
      else
      {
        std::ostringstream error_message;
        error_message << "Face indices not set up for boundary " << b
                      << " in region " << r << "\n";
        error_message << "This probably means that the boundary is not "
                         "adjacent to region\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }

    /// \short Return pointer to the current polyline that describes
    /// the b-th mesh boundary
    TriangleMeshPolyLine* boundary_polyline_pt(const unsigned& b)
    {
      std::map<unsigned, TriangleMeshCurveSection*>::iterator it =
        Boundary_curve_section_pt.find(b);
      // Search for the polyline associated with the given boundary
      if (it != Boundary_curve_section_pt.end())
      {
        return dynamic_cast<TriangleMeshPolyLine*>(
          Boundary_curve_section_pt[b]);
      }
      // If the boundary was not established then return 0, or a null pointer
      return 0;
    }

    /// \short Gets a pointer to a set with all the nodes related with a
    /// boundary
    std::map<unsigned, std::set<Node*>>& nodes_on_boundary_pt()
    {
      return Nodes_on_boundary_pt;
    }

    /// \short Gets the vertex number on the destination polyline (used
    /// to create the connections among shared boundaries)
    const bool get_connected_vertex_number_on_destination_polyline(
      TriangleMeshPolyLine* dst_polyline_pt,
      Vector<double>& vertex_coordinates,
      unsigned& vertex_number);

    /// \short Sort the polylines coming from joining them. Check whether
    /// it is necessary to reverse them or not. Used when joining two curve
    /// polylines in order to create a polygon
    void check_contiguousness_on_polylines_helper(
      Vector<TriangleMeshPolyLine*>& polylines_pt, unsigned& index);

    /// \short Sort the polylines coming from joining them. Check whether
    /// it is necessary to reverse them or not. Used when joining polylines
    /// and they still do not create a polygon
    void check_contiguousness_on_polylines_helper(
      Vector<TriangleMeshPolyLine*>& polylines_pt,
      unsigned& index_halo_start,
      unsigned& index_halo_end);

    /// Helper function that checks if a given point is inside a polygon
    /// (a set of sorted vertices that connected create a polygon)
    bool is_point_inside_polygon_helper(Vector<Vector<double>> polygon_vertices,
                                        Vector<double> point);

    /// Enables the creation of points (by Triangle) on the outer and
    /// internal boundaries
    void enable_automatic_creation_of_vertices_on_boundaries()
    {
      Allow_automatic_creation_of_vertices_on_boundaries = true;
    }

    /// Disables the creation of points (by Triangle) on the outer and
    /// internal boundaries
    void disable_automatic_creation_of_vertices_on_boundaries()
    {
      Allow_automatic_creation_of_vertices_on_boundaries = false;
    }

    /// Returns the status of the variable
    /// Allow_automatic_creation_of_vertices_on_boundaries
    bool is_automatic_creation_of_vertices_on_boundaries_allowed()
    {
      return Allow_automatic_creation_of_vertices_on_boundaries;
    }

#ifdef OOMPH_HAS_MPI

    /// Flush the boundary segment node storage
    void flush_boundary_segment_node(const unsigned& b)
    {
      Boundary_segment_node_pt[b].clear();
    }

    /// Set the number of segments associated with a boundary
    void set_nboundary_segment_node(const unsigned& b, const unsigned& s)
    {
      Boundary_segment_node_pt[b].resize(s);
    }

    /// Return the number of segments associated with a boundary
    unsigned nboundary_segment(const unsigned& b)
    {
      return Boundary_segment_node_pt[b].size();
    }

    /// Return the number of segments associated with a boundary
    unsigned long nboundary_segment_node(const unsigned& b)
    {
      unsigned ntotal_nodes = 0;
      unsigned nsegments = Boundary_segment_node_pt[b].size();
      for (unsigned is = 0; is < nsegments; is++)
      {
        ntotal_nodes += nboundary_segment_node(b, is);
      }
      return ntotal_nodes;
    }

    /// Return the number of nodes associated with a given segment of a
    /// boundary
    unsigned long nboundary_segment_node(const unsigned& b, const unsigned& s)
    {
      return Boundary_segment_node_pt[b][s].size();
    }

    /// Add the node node_pt to the b-th boundary and the s-th segment of
    /// the mesh
    void add_boundary_segment_node(const unsigned& b,
                                   const unsigned& s,
                                   Node* const& node_pt)
    {
      // Get the size of the Boundary_node_pt vector
      unsigned nbound_seg_node = nboundary_segment_node(b, s);
      bool node_already_on_this_boundary_segment = false;

      // Loop over the vector
      for (unsigned n = 0; n < nbound_seg_node; n++)
      {
        // Is the current node here already?
        if (node_pt == Boundary_segment_node_pt[b][s][n])
        {
          node_already_on_this_boundary_segment = true;
        }
      }

      // Add the base node pointer to the vector if it's not there already
      if (!node_already_on_this_boundary_segment)
      {
        Boundary_segment_node_pt[b][s].push_back(node_pt);
      }
    }

    /// \short Flag used at the setup_boundary_coordinate function to know
    /// if initial zeta values for segments have been assigned
    std::map<unsigned, bool> Assigned_segments_initial_zeta_values;

    /// \short Return direct access to the initial coordinates of a boundary
    std::map<unsigned, Vector<double>>& boundary_initial_coordinate()
    {
      return Boundary_initial_coordinate;
    }

    /// \short Return direct access to the final coordinates of a boundary
    std::map<unsigned, Vector<double>>& boundary_final_coordinate()
    {
      return Boundary_final_coordinate;
    }

    /// \short Return direct access to the initial zeta coordinate of a
    /// boundary
    std::map<unsigned, Vector<double>>& boundary_initial_zeta_coordinate()
    {
      return Boundary_initial_zeta_coordinate;
    }

    /// \short Return direct access to the final zeta coordinates of a
    /// boundary
    std::map<unsigned, Vector<double>>& boundary_final_zeta_coordinate()
    {
      return Boundary_final_zeta_coordinate;
    }

    /// \short Return the info. to know if it is necessary to reverse the
    /// segment based on a previous mesh
    std::map<unsigned, Vector<unsigned>>& boundary_segment_inverted()
    {
      return Boundary_segment_inverted;
    }

    /// \short Return direct access to the initial coordinates for the
    /// segments that are part of a boundary
    std::map<unsigned, Vector<Vector<double>>>& boundary_segment_initial_coordinate()
    {
      return Boundary_segment_initial_coordinate;
    }

    /// \short Return direct access to the final coordinates for the
    /// segments that are part of a boundary
    std::map<unsigned, Vector<Vector<double>>>& boundary_segment_final_coordinate()
    {
      return Boundary_segment_final_coordinate;
    }

    /// \short Return direct access to the initial arclength for the
    /// segments that are part of a boundary
    std::map<unsigned, Vector<double>>& boundary_segment_initial_arclength()
    {
      return Boundary_segment_initial_arclength;
    }

    /// \short Return direct access to the final arclength for the
    /// segments that are part of a boundary
    std::map<unsigned, Vector<double>>& boundary_segment_final_arclength()
    {
      return Boundary_segment_final_arclength;
    }

    /// \short Return direct access to the initial zeta for the
    /// segments that are part of a boundary
    std::map<unsigned, Vector<double>>& boundary_segment_initial_zeta()
    {
      return Boundary_segment_initial_zeta;
    }

    /// \short Return direct access to the final zeta for the
    /// segments that are part of a boundary
    std::map<unsigned, Vector<double>>& boundary_segment_final_zeta()
    {
      return Boundary_segment_final_zeta;
    }

    /// \short Return the initial zeta for the segments that are
    /// part of a boundary
    Vector<double>& boundary_segment_initial_zeta(const unsigned& b)
    {
      std::map<unsigned, Vector<double>>::iterator it =
        Boundary_segment_initial_zeta.find(b);

#ifdef PARANOID

      if (it == Boundary_segment_initial_zeta.end())
      {
        std::stringstream error_message;
        error_message << "The boundary (" << b
                      << ") has no segments associated with it!!\n\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif // PARANOID

      return (*it).second;
    }

    /// \short Return the final zeta for the segments that are
    /// part of a boundary
    Vector<double>& boundary_segment_final_zeta(const unsigned& b)
    {
      std::map<unsigned, Vector<double>>::iterator it =
        Boundary_segment_final_zeta.find(b);

#ifdef PARANOID

      if (it == Boundary_segment_final_zeta.end())
      {
        std::stringstream error_message;
        error_message << "The boundary (" << b
                      << ") has no segments associated with it!!\n\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif // PARANOID

      return (*it).second;
    }

    /// \short Return the initial coordinates for the boundary
    Vector<double>& boundary_initial_coordinate(const unsigned& b)
    {
      std::map<unsigned, Vector<double>>::iterator it =
        Boundary_initial_coordinate.find(b);

#ifdef PARANOID

      if (it == Boundary_initial_coordinate.end())
      {
        std::stringstream error_message;
        error_message << "The boundary (" << b
                      << ") has not established initial coordinates\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif
      return (*it).second;
    }

    /// \short Return the final coordinates for the boundary
    Vector<double>& boundary_final_coordinate(const unsigned& b)
    {
      std::map<unsigned, Vector<double>>::iterator it =
        Boundary_final_coordinate.find(b);

#ifdef PARANOID

      if (it == Boundary_final_coordinate.end())
      {
        std::stringstream error_message;
        error_message << "The boundary (" << b
                      << ") has not established final coordinates\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif

      return (*it).second;
    }

    /// \short Return the info. to know if it is necessary to reverse the
    /// segment based on a previous mesh
    const Vector<unsigned> boundary_segment_inverted(const unsigned& b) const
    {
      std::map<unsigned, Vector<unsigned>>::const_iterator it =
        Boundary_segment_inverted.find(b);

#ifdef PARANOID

      if (it == Boundary_segment_inverted.end())
      {
        std::stringstream error_message;
        error_message << "The boundary (" << b
                      << ") has not established inv. segments info\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif

      return (*it).second;
    }

    /// \short Return the info. to know if it is necessary to reverse the
    /// segment based on a previous mesh
    Vector<unsigned>& boundary_segment_inverted(const unsigned& b)
    {
      std::map<unsigned, Vector<unsigned>>::iterator it =
        Boundary_segment_inverted.find(b);

#ifdef PARANOID

      if (it == Boundary_segment_inverted.end())
      {
        std::stringstream error_message;
        error_message << "The boundary (" << b
                      << ") has not established inv. segments info\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif

      return (*it).second;
    }

    /// \short Return the initial zeta coordinate for the boundary
    Vector<double>& boundary_initial_zeta_coordinate(const unsigned& b)
    {
      std::map<unsigned, Vector<double>>::iterator it =
        Boundary_initial_zeta_coordinate.find(b);

#ifdef PARANOID

      if (it == Boundary_initial_zeta_coordinate.end())
      {
        std::stringstream error_message;
        error_message << "The boundary (" << b
                      << ") has not established initial zeta "
                      << "coordinate\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif

      return (*it).second;
    }

    /// \short Return the final zeta coordinate for the boundary
    Vector<double>& boundary_final_zeta_coordinate(const unsigned& b)
    {
      std::map<unsigned, Vector<double>>::iterator it =
        Boundary_final_zeta_coordinate.find(b);

#ifdef PARANOID

      if (it == Boundary_final_zeta_coordinate.end())
      {
        std::stringstream error_message;
        error_message << "The boundary (" << b
                      << ") has not established final zeta coordinate\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif

      return (*it).second;
    }

    /// \short Return the initial arclength for the segments that are
    /// part of a boundary
    Vector<double>& boundary_segment_initial_arclength(const unsigned& b)
    {
      std::map<unsigned, Vector<double>>::iterator it =
        Boundary_segment_initial_arclength.find(b);

#ifdef PARANOID

      if (it == Boundary_segment_initial_arclength.end())
      {
        std::stringstream error_message;
        error_message << "The boundary (" << b
                      << ") has no segments associated with it!!\n\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif

      return (*it).second;
    }

    /// \short Return the final arclength for the segments that are
    /// part of a boundary
    Vector<double>& boundary_segment_final_arclength(const unsigned& b)
    {
      std::map<unsigned, Vector<double>>::iterator it =
        Boundary_segment_final_arclength.find(b);

#ifdef PARANOID

      if (it == Boundary_segment_final_arclength.end())
      {
        std::stringstream error_message;
        error_message << "The boundary (" << b
                      << ") has no segments associated with it!!\n\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif

      return (*it).second;
    }

    /// \short Return the initial coordinates for the segments that are
    /// part of a boundary
    Vector<Vector<double>>& boundary_segment_initial_coordinate(
      const unsigned& b)
    {
      std::map<unsigned, Vector<Vector<double>>>::iterator it =
        Boundary_segment_initial_coordinate.find(b);

#ifdef PARANOID

      if (it == Boundary_segment_initial_coordinate.end())
      {
        std::stringstream error_message;
        error_message << "The boundary (" << b
                      << ") has no segments associated with it!!\n\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif

      return (*it).second;
    }

    /// \short Return the final coordinates for the segments that are
    /// part of a boundary
    Vector<Vector<double>>& boundary_segment_final_coordinate(const unsigned& b)
    {
      std::map<unsigned, Vector<Vector<double>>>::iterator it =
        Boundary_segment_final_coordinate.find(b);

#ifdef PARANOID

      if (it == Boundary_segment_final_coordinate.end())
      {
        std::stringstream error_message;
        error_message << "The boundary (" << b
                      << ") has no segments associated with it!!\n\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif

      return (*it).second;
    }

#endif // OOMPH_HAS_MPI

    /// \short Setup boundary coordinate on boundary b.
    /// Boundary coordinate increases continously along
    /// polygonal boundary. It's zero at the lowest left
    /// node on the boundary.
    template<class ELEMENT>
    void setup_boundary_coordinates(const unsigned& b)
    {
      // Dummy file
      std::ofstream some_file;
      setup_boundary_coordinates<ELEMENT>(b, some_file);
    }

    /// \short Setup boundary coordinate on boundary b. Doc Faces
    /// in outfile.
    /// Boundary coordinate increases continously along
    /// polygonal boundary. It's zero at the lowest left
    /// node on the boundary.
    template<class ELEMENT>
    void setup_boundary_coordinates(const unsigned& b, std::ofstream& outfile);


  protected:
#ifdef OOMPH_HAS_TRIANGLE_LIB

    /// \short Create TriangulateIO object from outer boundaries,
    /// internal boundaries, and open curves. Add the holes and regions
    /// information as well
    void build_triangulateio(
      Vector<TriangleMeshPolygon*>& outer_polygons_pt,
      Vector<TriangleMeshPolygon*>& internal_polygons_pt,
      Vector<TriangleMeshOpenCurve*>& open_curves_pt,
      Vector<Vector<double>>& extra_holes_coordinates,
      std::map<unsigned, Vector<double>>& regions_coordinates,
      std::map<unsigned, double>& regions_areas,
      TriangulateIO& triangulate_io);

    /// \short Data structure filled when the connection matrix is created, for
    /// each polyline, there are two vertex_connection_info structures,
    /// one for each end
    struct vertex_connection_info
    {
      bool is_connected;
      unsigned boundary_id_to_connect;
      unsigned boundary_chunk_to_connect;
      unsigned vertex_number_to_connect;
    }; // vertex_connection_info

    /// \short Data structure to store the base vertex info, initial or final
    /// vertex in the polylines have an associated base vertex
    struct base_vertex_info
    {
      bool has_base_vertex_assigned;
      bool is_base_vertex;
      unsigned boundary_id;
      unsigned boundary_chunk;
      unsigned vertex_number;
    }; // base_vertex_info

    /// \short Helps to add information to the connection matrix of the
    /// given polyline
    void add_connection_matrix_info_helper(
      TriangleMeshPolyLine* polyline_pt,
      std::map<unsigned, std::map<unsigned, Vector<vertex_connection_info>>>&
        connection_matrix,
      TriangleMeshPolyLine* next_polyline_pt = 0);

    /// \short Initialise the base vertex structure, set every vertex to
    /// no visited and not being a base vertex
    void initialise_base_vertex(
      TriangleMeshPolyLine* polyline_pt,
      std::map<unsigned, std::map<unsigned, Vector<base_vertex_info>>>&
        base_vertices);

    /// \short Helps to identify the base vertex of the given polyline
    void add_base_vertex_info_helper(
      TriangleMeshPolyLine* polyline_pt,
      std::map<unsigned, std::map<unsigned, Vector<base_vertex_info>>>&
        base_vertices,
      std::map<unsigned, std::map<unsigned, Vector<vertex_connection_info>>>&
        connection_matrix,
      std::map<unsigned, std::map<unsigned, unsigned>>&
        boundary_chunk_nvertices);

#endif

#ifdef OOMPH_HAS_MPI

    /// \short Used to store the nodes associated to a boundary and to an
    /// specific segment (this only applies in distributed meshes where the
    /// boundary is splitted in segments)
    std::map<unsigned, Vector<Vector<Node*>>> Boundary_segment_node_pt;

    /// \short Stores the initial zeta coordinate for the segments that
    /// appear when a boundary is splited among processors
    std::map<unsigned, Vector<double>> Boundary_segment_initial_zeta;

    /// \short Stores the final zeta coordinate for the segments that
    /// appear when a boundary is splited among processors
    std::map<unsigned, Vector<double>> Boundary_segment_final_zeta;

    /// \short Stores the initial coordinates for the boundary
    std::map<unsigned, Vector<double>> Boundary_initial_coordinate;

    /// \short Stores the final coordinates for the boundary
    std::map<unsigned, Vector<double>> Boundary_final_coordinate;

    /// \short Stores the info. to know if it is necessary to reverse the
    /// segment based on a previous mesh
    std::map<unsigned, Vector<unsigned>> Boundary_segment_inverted;

    /// \short Stores the initial zeta coordinate for the boundary
    std::map<unsigned, Vector<double>> Boundary_initial_zeta_coordinate;

    /// \short Stores the final zeta coordinate for the boundary
    std::map<unsigned, Vector<double>> Boundary_final_zeta_coordinate;

    /// \short Stores the initial arclength for the segments that appear when
    /// a boundary is splited among processors
    std::map<unsigned, Vector<double>> Boundary_segment_initial_arclength;

    /// \short Stores the final arclength for the segments that appear when
    /// a boundary is splited among processors
    std::map<unsigned, Vector<double>> Boundary_segment_final_arclength;

    /// \short Stores the initial coordinates for the segments that appear
    /// when a boundary is splited among processors
    std::map<unsigned, Vector<Vector<double>>>
      Boundary_segment_initial_coordinate;

    /// \short Stores the final coordinates for the segments that appear
    /// when a boundary is splited among processors
    std::map<unsigned, Vector<Vector<double>>>
      Boundary_segment_final_coordinate;

#endif

    /// \short Flag to indicate whether the automatic creation of vertices
    /// along the boundaries by Triangle is allowed
    bool Allow_automatic_creation_of_vertices_on_boundaries;

    /// \short Snap the boundary nodes onto any curvilinear boundaries
    /// defined by geometric objects
    void snap_nodes_onto_geometric_objects();

    /// \short Vector of elements in each region differentiated by attribute
    /// (the key of the map is the attribute)
    std::map<unsigned, Vector<FiniteElement*>> Region_element_pt;

    /// Vector of attributes associated with the elements in each region
    Vector<double> Region_attribute;

    /// \short Storage for the geometric objects associated with any boundaries
    std::map<unsigned, GeomObject*> Boundary_geom_object_pt;

    /// Storage for the limits of the boundary coordinates defined by the use
    /// of geometric objects. Only used for curvilinear boundaries.
    std::map<unsigned, Vector<double>> Boundary_coordinate_limits;

    /// Polygon that defines outer boundaries
    Vector<TriangleMeshPolygon*> Outer_boundary_pt;

    /// Vector of polygons that define internal polygons
    Vector<TriangleMeshPolygon*> Internal_polygon_pt;

    /// Vector of open polylines that define internal curves
    Vector<TriangleMeshOpenCurve*> Internal_open_curve_pt;

    /// Storage for extra coordinates for holes
    Vector<Vector<double>> Extra_holes_coordinates;

    /// Storage for extra coordinates for regions. The key on the map
    /// is the region id
    std::map<unsigned, Vector<double>> Regions_coordinates;

    /// A map that stores the associated curve section of the specified boundary
    /// id
    std::map<unsigned, TriangleMeshCurveSection*> Boundary_curve_section_pt;

    /// Storage for elements adjacent to a boundary in a particular region
    Vector<std::map<unsigned, Vector<FiniteElement*>>>
      Boundary_region_element_pt;

    /// Storage for the face index adjacent to a boundary in a particular region
    Vector<std::map<unsigned, Vector<int>>> Face_index_region_at_boundary;

    /// \short Storage for pairs of doubles representing:
    /// .first: the arclength along the polygonal representation of
    ///         the curviline
    /// .second: the corresponding intrinsic coordinate on the associated
    ///          geometric object
    /// at which the vertices on the specified boundary are located.
    /// Only used for boundaries represented by geom objects.
    std::map<unsigned, Vector<std::pair<double, double>>>
      Polygonal_vertex_arclength_info;

    /// \short Stores a pointer to a set with all the nodes
    /// related with a boundary
    std::map<unsigned, std::set<Node*>> Nodes_on_boundary_pt;

    /// \short A set that contains the curve sections created by this object
    /// therefore it is necessary to free their associated memory
    std::set<TriangleMeshCurveSection*> Free_curve_section_pt;

    /// \short A set that contains the polygons created by this object
    /// therefore it is necessary to free their associated memory
    std::set<TriangleMeshPolygon*> Free_polygon_pt;

    /// \short A set that contains the open curves created by this
    /// object therefore it is necessary to free their associated memory
    std::set<TriangleMeshOpenCurve*> Free_open_curve_pt;

    /// \short Helper function to copy the connection information from
    /// the input curve(polyline or curviline) to the output polyline
    void copy_connection_information(TriangleMeshCurveSection* input_curve_pt,
                                     TriangleMeshCurveSection* output_curve_pt);

    /// \short Helper function to copy the connection information from
    /// the input curve(polyline or curviline) to the output sub-polyline
    void copy_connection_information_to_sub_polylines(
      TriangleMeshCurveSection* input_curve_pt,
      TriangleMeshCurveSection* output_curve_pt);


#ifdef PARANOID

    // Used to verify if any of the polygons (closedcurves) that define
    // the mesh are of type ImmersedRigidBodyTriangleMeshPolygon, if
    // that is the case it may lead to problems in case of using load
    // balance
    bool Immersed_rigid_body_triangle_mesh_polygon_used;

#endif

#ifdef OOMPH_HAS_TRIANGLE_LIB

    /// \short Helper function to create polyline vertex coordinates for
    /// curvilinear boundary specified by boundary_pt, using either
    /// equal increments in zeta or in (approximate) arclength
    /// along the curviline. vertex_coord[i_vertex][i_dim] stores
    /// i_dim-th coordinate of i_vertex-th vertex.
    /// polygonal_vertex_arclength_info[i_vertex] contains the pair of doubles
    /// made of the arclength of the i_vertex-th vertex along the polygonal
    /// representation (.first), and the corresponding coordinate on the
    /// GeomObject (.second)
    void create_vertex_coordinates_for_polyline_no_connections(
      TriangleMeshCurviLine* boundary_pt,
      Vector<Vector<double>>& vertex_coord,
      Vector<std::pair<double, double>>& polygonal_vertex_arclength_info)
    {
      // Intrinsic coordinate along GeomObjects
      Vector<double> zeta(1);

      // Position vector to point on GeomObject
      Vector<double> posn(2);

      // Start coordinate
      double zeta_initial = boundary_pt->zeta_start();

      // How many segments do we want on this polyline?
      unsigned n_seg = boundary_pt->nsegment();
      vertex_coord.resize(n_seg + 1);
      polygonal_vertex_arclength_info.resize(n_seg + 1);
      polygonal_vertex_arclength_info[0].first = 0.0;
      polygonal_vertex_arclength_info[0].second = zeta_initial;

      // Vertices placed in equal zeta increments
      if (!(boundary_pt->space_vertices_evenly_in_arclength()))
      {
        // Read the values of the limiting coordinates, assuming equal
        // spacing of the nodes
        double zeta_increment =
          (boundary_pt->zeta_end() - boundary_pt->zeta_start()) /
          (double(n_seg));

        // Loop over the n_seg+1 points bounding the segments
        for (unsigned s = 0; s < n_seg + 1; s++)
        {
          // Get the coordinates
          zeta[0] = zeta_initial + zeta_increment * double(s);
          boundary_pt->geom_object_pt()->position(zeta, posn);
          vertex_coord[s] = posn;

          // Bump up the polygonal arclength
          if (s > 0)
          {
            polygonal_vertex_arclength_info[s].first =
              polygonal_vertex_arclength_info[s - 1].first +
              sqrt(pow(vertex_coord[s][0] - vertex_coord[s - 1][0], 2) +
                   pow(vertex_coord[s][1] - vertex_coord[s - 1][1], 2));
            polygonal_vertex_arclength_info[s].second = zeta[0];
          }
        }
      }
      // Vertices placed in equal increments in (approximate) arclength
      else
      {
        // Number of sampling points to compute arclength and
        // arclength increments
        unsigned nsample_per_segment = 100;
        unsigned nsample = nsample_per_segment * n_seg;

        // Work out start and increment
        double zeta_increment =
          (boundary_pt->zeta_end() - boundary_pt->zeta_start()) /
          (double(nsample));

        // Get coordinate of first point
        Vector<double> start_point(2);
        zeta[0] = zeta_initial;

        boundary_pt->geom_object_pt()->position(zeta, start_point);

        // Storage for coordinates of end point
        Vector<double> end_point(2);

        // Compute total arclength
        double total_arclength = 0.0;
        for (unsigned i = 1; i < nsample; i++)
        {
          // Next point
          zeta[0] += zeta_increment;

          // Get coordinate of end point
          boundary_pt->geom_object_pt()->position(zeta, end_point);

          // Increment arclength
          total_arclength += sqrt(pow(end_point[0] - start_point[0], 2) +
                                  pow(end_point[1] - start_point[1], 2));

          // Shift back
          start_point = end_point;
        }

        // Desired arclength increment
        double target_s_increment = total_arclength / (double(n_seg));

        // Get coordinate of first point again
        zeta[0] = zeta_initial;
        boundary_pt->geom_object_pt()->position(zeta, start_point);

        // Assign as coordinate
        vertex_coord[0] = start_point;

        // Start sampling point
        unsigned i_lo = 1;

        // Loop over the n_seg-1 internal points bounding the segments
        for (unsigned s = 1; s < n_seg; s++)
        {
          // Visit potentially all sample points until we've found
          // the one at which we exceed the target arclength increment
          double arclength_increment = 0.0;
          for (unsigned i = i_lo; i < nsample; i++)
          {
            // Next point
            zeta[0] += zeta_increment;

            // Get coordinate of end point
            boundary_pt->geom_object_pt()->position(zeta, end_point);

            // Increment arclength increment
            arclength_increment += sqrt(pow(end_point[0] - start_point[0], 2) +
                                        pow(end_point[1] - start_point[1], 2));

            // Shift back
            start_point = end_point;

            // Are we there yet?
            if (arclength_increment > target_s_increment)
            {
              // Remember how far we've got
              i_lo = i;

              // And bail out
              break;
            }
          }

          // Store the coordinates
          vertex_coord[s] = end_point;

          // Bump up the polygonal arclength
          if (s > 0)
          {
            polygonal_vertex_arclength_info[s].first =
              polygonal_vertex_arclength_info[s - 1].first +
              sqrt(pow(vertex_coord[s][0] - vertex_coord[s - 1][0], 2) +
                   pow(vertex_coord[s][1] - vertex_coord[s - 1][1], 2));
            polygonal_vertex_arclength_info[s].second = zeta[0];
          }
        }

        // Final point
        unsigned s = n_seg;
        zeta[0] = boundary_pt->zeta_end();
        boundary_pt->geom_object_pt()->position(zeta, end_point);
        vertex_coord[s] = end_point;
        polygonal_vertex_arclength_info[s].first =
          polygonal_vertex_arclength_info[s - 1].first +
          sqrt(pow(vertex_coord[s][0] - vertex_coord[s - 1][0], 2) +
               pow(vertex_coord[s][1] - vertex_coord[s - 1][1], 2));
        polygonal_vertex_arclength_info[s].second = zeta[0];
      }
    }

    /// \short Helper function to create polyline vertex coordinates for
    /// curvilinear boundary specified by boundary_pt, using either
    /// equal increments in zeta or in (approximate) arclength
    /// along the curviline. vertex_coord[i_vertex][i_dim] stores
    /// i_dim-th coordinate of i_vertex-th vertex.
    /// polygonal_vertex_arclength_info[i_vertex] contains the pair of doubles
    /// made of the arclength of the i_vertex-th vertex along the polygonal
    /// representation (.first), and the corresponding coordinate on the
    /// GeomObject (.second)
    void create_vertex_coordinates_for_polyline_connections(
      TriangleMeshCurviLine* boundary_pt,
      Vector<Vector<double>>& vertex_coord,
      Vector<std::pair<double, double>>& polygonal_vertex_arclength_info)
    {
      // Start coordinate
      double zeta_initial = boundary_pt->zeta_start();
      // Final coordinate
      double zeta_final = boundary_pt->zeta_end();

      Vector<double>* connection_points_pt =
        boundary_pt->connection_points_pt();

      unsigned n_connections = connection_points_pt->size();

      // We need to sort the connection points
      if (n_connections > 1)
      {
        std::sort(connection_points_pt->begin(), connection_points_pt->end());
      }

#ifdef PARANOID
      // Are the connection points out of range of the polyline
      bool out_of_range_connection_points = false;
      std::ostringstream error_message;
      // Check if the curviline should be created on a reversed way
      bool reversed = false;
      if (zeta_final < zeta_initial)
      {
        reversed = true;
      }
      if (!reversed)
      {
        if (zeta_initial > (*connection_points_pt)[0])
        {
          error_message
            << "One of the specified connection points is out of the\n"
            << "curviline limits. We found that the point ("
            << (*connection_points_pt)[0] << ") is\n"
            << "less than the"
            << "initial s value which is (" << zeta_initial << ").\n"
            << "Initial value: (" << zeta_initial << ")\n"
            << "Final value: (" << zeta_final << ")\n"
            << std::endl;
          out_of_range_connection_points = true;
        }

        if (zeta_final < (*connection_points_pt)[n_connections - 1])
        {
          error_message
            << "One of the specified connection points is out of the\n"
            << "curviline limits. We found that the point ("
            << (*connection_points_pt)[n_connections - 1] << ") is\n"
            << "greater than the final s value which is (" << zeta_final
            << ").\n"
            << "Initial value: (" << zeta_initial << ")\n"
            << "Final value: (" << zeta_final << ")\n"
            << std::endl;
          out_of_range_connection_points = true;
        }
      }
      else
      {
        if (zeta_initial < (*connection_points_pt)[0])
        {
          error_message
            << "One of the specified connection points is out of the\n"
            << "curviline limits. We found that the point ("
            << (*connection_points_pt)[0] << ") is\n"
            << "greater than the"
            << "initial s value which is (" << zeta_initial << ").\n"
            << "Initial value: (" << zeta_initial << ")\n"
            << "Final value: (" << zeta_final << ")\n"
            << std::endl;
          out_of_range_connection_points = true;
        }

        if (zeta_final > (*connection_points_pt)[n_connections - 1])
        {
          error_message
            << "One of the specified connection points is out of the\n"
            << "curviline limits. We found that the point ("
            << (*connection_points_pt)[n_connections - 1] << ") is\n"
            << "less than the final s value which is (" << zeta_final << ").\n"
            << "Initial value: (" << zeta_initial << ")\n"
            << "Final value: (" << zeta_final << ")\n"
            << std::endl;
          out_of_range_connection_points = true;
        }
      }

      if (out_of_range_connection_points)
      {
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif // PARANOID

      // Intrinsic coordinate along GeomObjects
      Vector<double> zeta(1);

      // Position vector to point on GeomObject
      Vector<double> posn(2);

      // How many segments do we want on this polyline?
      unsigned n_seg = boundary_pt->nsegment();

      // How many connection vertices have we already created
      unsigned i_connection = 0;
      Vector<double> zeta_connection(1);

      // If we have more connection points than the generated
      // by the number of segments then we have to change the
      // number of segments and create all the vertices
      // according to the connection points list
      if (n_connections >= n_seg - 1)
      {
        std::ostringstream warning_message;
        std::string output_string = "UnstructuredTwoDMeshGeometryBase::";
        output_string += "create_vertex_coordinates_for_polyline_connections()";

        warning_message
          << "The number of segments specified for the curviline with\n"
          << "boundary id (" << boundary_pt->boundary_id() << ") is less "
          << "(or equal) than the ones that will be\ngenerated by using "
          << "the specified number of connection points.\n"
          << "You specified (" << n_seg << ") segments but ("
          << n_connections + 1 << ") segments\nwill be generated." << std::endl;
        OomphLibWarning(
          warning_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);

        // We have to explicitly change the number of segments
        boundary_pt->nsegment() = n_connections + 1;
        n_seg = boundary_pt->nsegment();
        vertex_coord.resize(n_seg + 1);

        // Initial coordinate and initial values
        zeta[0] = zeta_initial;
        boundary_pt->geom_object_pt()->position(zeta, posn);
        vertex_coord[0] = posn;

        polygonal_vertex_arclength_info.resize(n_seg + 1);
        polygonal_vertex_arclength_info[0].first = 0.0;
        polygonal_vertex_arclength_info[0].second = zeta_initial;

        // Loop over the n_connections points bounding the segments
        for (i_connection = 0; i_connection < n_connections; i_connection++)
        {
          // Get the coordinates
          zeta[0] = (*connection_points_pt)[i_connection];
          boundary_pt->geom_object_pt()->position(zeta, posn);
          vertex_coord[i_connection + 1] = posn;

          // Bump up the polygonal arclength
          polygonal_vertex_arclength_info[i_connection + 1].first =
            polygonal_vertex_arclength_info[i_connection].first +
            sqrt(pow(vertex_coord[i_connection + 1][0] -
                       vertex_coord[i_connection][0],
                     2) +
                 pow(vertex_coord[i_connection + 1][1] -
                       vertex_coord[i_connection][1],
                     2));
          polygonal_vertex_arclength_info[i_connection + 1].second = zeta[0];
        }

        // Final coordinate and final values
        zeta[0] = zeta_final;
        boundary_pt->geom_object_pt()->position(zeta, posn);
        vertex_coord[n_seg] = posn;

        polygonal_vertex_arclength_info[n_seg].first =
          polygonal_vertex_arclength_info[n_seg - 1].first +
          sqrt(pow(vertex_coord[n_seg][0] - vertex_coord[n_seg - 1][0], 2) +
               pow(vertex_coord[n_seg][1] - vertex_coord[n_seg - 1][1], 2));
        polygonal_vertex_arclength_info[n_seg].second = zeta_final;
      }
      else
      {
        // Total number of vertices
        unsigned n_t_vertices = n_seg + 1;

        // Number of vertices left for creation
        unsigned l_vertices = n_t_vertices;

        // Total number of already created vertices
        unsigned n_assigned_vertices = 0;

        // Stores the distance between current vertices in the list
        // Edge vertices + Connection points - 1
        Vector<double> delta_z(2 + n_connections - 1);

        std::list<double> zeta_values_pt;
        zeta_values_pt.push_back(zeta_initial);
        for (unsigned s = 0; s < n_connections; s++)
        {
          zeta_values_pt.push_back((*connection_points_pt)[s]);
        }
        zeta_values_pt.push_back(zeta_final);

        l_vertices -= 2; // Edge vertices
        l_vertices -= n_connections; // Connection points
        n_assigned_vertices += 2; // Edge vertices
        n_assigned_vertices += n_connections; // Connection points

        // Vertices placed in equal zeta increments
        if (!(boundary_pt->space_vertices_evenly_in_arclength()))
        {
          double local_zeta_initial;
          double local_zeta_final;
          double local_zeta_increment;
          double local_zeta_insert;

          // How many vertices for each section
          unsigned local_n_vertices;

          std::list<double>::iterator l_it = zeta_values_pt.begin();
          std::list<double>::iterator r_it = zeta_values_pt.begin();
          r_it++;

          for (unsigned h = 0; r_it != zeta_values_pt.end();
               l_it++, r_it++, h++)
          {
            delta_z[h] = *r_it - *l_it;
          }

          l_it = r_it = zeta_values_pt.begin();
          r_it++;

          for (unsigned h = 0; r_it != zeta_values_pt.end(); h++)
          {
            local_n_vertices =
              static_cast<unsigned>(((double)n_t_vertices * delta_z[h]) /
                                    std::fabs(zeta_final - zeta_initial));

            local_zeta_initial = *l_it;
            local_zeta_final = *r_it;
            local_zeta_increment = (local_zeta_final - local_zeta_initial) /
                                   (double)(local_n_vertices + 1);

            for (unsigned s = 0; s < local_n_vertices; s++)
            {
              local_zeta_insert =
                local_zeta_initial + local_zeta_increment * double(s + 1);
              zeta_values_pt.insert(r_it, local_zeta_insert);
              n_assigned_vertices++;
            }
            // Moving to the next segment
            l_it = r_it;
            r_it++;
          }

          // Finishing it ...!!!
#ifdef PARANOID
          // Counting the vertices number and the total of
          // assigned vertices values
          unsigned s = zeta_values_pt.size();

          if (s != n_assigned_vertices)
          {
            error_message
              << "The total number of assigned vertices is different from\n"
              << "the number of elements in the z_values list. The number"
              << "of\nelements in the z_values list is (" << s << ") but "
              << "the number\n"
              << "of assigned vertices is (" << n_assigned_vertices << ")."
              << std::endl
              << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif // PARANOID

          vertex_coord.resize(n_assigned_vertices);
          polygonal_vertex_arclength_info.resize(n_assigned_vertices);
          polygonal_vertex_arclength_info[0].first = 0.0;
          polygonal_vertex_arclength_info[0].second = zeta_initial;

          // Creating the vertices with the corresponding z_values
          l_it = zeta_values_pt.begin();
          for (unsigned s = 0; l_it != zeta_values_pt.end(); s++, l_it++)
          {
            // Get the coordinates
            zeta[0] = *l_it;
            boundary_pt->geom_object_pt()->position(zeta, posn);
            vertex_coord[s] = posn;

            // Bump up the polygonal arclength
            if (s > 0)
            {
              polygonal_vertex_arclength_info[s].first =
                polygonal_vertex_arclength_info[s - 1].first +
                sqrt(pow(vertex_coord[s][0] - vertex_coord[s - 1][0], 2) +
                     pow(vertex_coord[s][1] - vertex_coord[s - 1][1], 2));
            }
          }
        }
        // Vertices placed in equal increments in (approximate) arclength
        else
        {
          // Compute the total arclength
          // Number of sampling points to compute arclength and
          // arclength increments
          unsigned nsample_per_segment = 100;
          unsigned nsample = nsample_per_segment * n_seg;

          // Work out start and increment
          double zeta_increment =
            (zeta_final - zeta_initial) / (double(nsample));

          // Get coordinate of first point
          Vector<double> start_point(2);
          zeta[0] = zeta_initial;
          boundary_pt->geom_object_pt()->position(zeta, start_point);

          // Storage for coordinates of end point
          Vector<double> end_point(2);

          // Compute total arclength
          double total_arclength = 0.0;
          for (unsigned i = 1; i < nsample; i++)
          {
            // Next point
            zeta[0] += zeta_increment;

            // Get coordinate of end point
            boundary_pt->geom_object_pt()->position(zeta, end_point);

            // Increment arclength
            total_arclength += sqrt(pow(end_point[0] - start_point[0], 2) +
                                    pow(end_point[1] - start_point[1], 2));

            // Shift back
            start_point = end_point;
          }

          double local_zeta_initial;
          double local_zeta_final;
          double local_zeta_increment;

          // How many vertices per section
          unsigned local_n_vertices;

          std::list<double>::iterator l_it = zeta_values_pt.begin();
          std::list<double>::iterator r_it = zeta_values_pt.begin();
          r_it++;

          for (unsigned h = 0; r_it != zeta_values_pt.end(); h++)
          {
            // There is no need to move the r_it iterator since it is
            // moved at the final of this loop
            local_zeta_initial = *l_it;
            local_zeta_final = *r_it;
            local_zeta_increment =
              (local_zeta_final - local_zeta_initial) / (double)(nsample);

            // Compute local arclength
            // Get coordinate of first point
            zeta[0] = local_zeta_initial;
            boundary_pt->geom_object_pt()->position(zeta, start_point);

            delta_z[h] = 0.0;

            for (unsigned i = 1; i < nsample; i++)
            {
              // Next point
              zeta[0] += local_zeta_increment;

              // Get coordinate of end point
              boundary_pt->geom_object_pt()->position(zeta, end_point);

              // Increment arclength
              delta_z[h] += sqrt(pow(end_point[0] - start_point[0], 2) +
                                 pow(end_point[1] - start_point[1], 2));

              // Shift back
              start_point = end_point;
            }

            local_n_vertices = static_cast<unsigned>(
              ((double)n_t_vertices * delta_z[h]) / (total_arclength));

            // Desired arclength increment
            double local_target_s_increment =
              delta_z[h] / double(local_n_vertices + 1);

            // Get coordinate of first point again
            zeta[0] = local_zeta_initial;
            boundary_pt->geom_object_pt()->position(zeta, start_point);

            // Start sampling point
            unsigned i_lo = 1;

            // Loop over the n_seg-1 internal points bounding the segments
            for (unsigned s = 0; s < local_n_vertices; s++)
            {
              // Visit potentially all sample points until we've found
              // the one at which we exceed the target arclength increment
              double local_arclength_increment = 0.0;
              for (unsigned i = i_lo; i < nsample; i++)
              // for (unsigned i=i_lo;i<nsample_per_segment;i++)
              {
                // Next point
                zeta[0] += local_zeta_increment;

                // Get coordinate of end point
                boundary_pt->geom_object_pt()->position(zeta, end_point);

                // Increment arclength increment
                local_arclength_increment +=
                  sqrt(pow(end_point[0] - start_point[0], 2) +
                       pow(end_point[1] - start_point[1], 2));

                // Shift back
                start_point = end_point;

                // Are we there yet?
                if (local_arclength_increment > local_target_s_increment)
                {
                  // Remember how far we've got
                  i_lo = i;

                  // And bail out
                  break;
                }
              }

              zeta_values_pt.insert(r_it, zeta[0]);
              n_assigned_vertices++;
            }
            // Moving to the next segments
            l_it = r_it;
            r_it++;
          }

          // Finishing it ... !!!
#ifdef PARANOID
          // Counting the vertices number and the total of
          // assigned vertices values
          unsigned h = zeta_values_pt.size();

          if (h != n_assigned_vertices)
          {
            error_message
              << "The total number of assigned vertices is different from\n"
              << "the number of elements in the z_values list. The number of\n"
              << "elements in the z_values list is (" << h
              << ") but the number\n"
              << "of assigned vertices is (" << n_assigned_vertices << ")."
              << std::endl
              << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif // PARANOID

          vertex_coord.resize(n_assigned_vertices);
          polygonal_vertex_arclength_info.resize(n_assigned_vertices);
          polygonal_vertex_arclength_info[0].first = 0.0;
          polygonal_vertex_arclength_info[0].second = zeta_initial;

          // Creating the vertices with the corresponding z_values
          l_it = zeta_values_pt.begin();
          for (unsigned s = 0; l_it != zeta_values_pt.end(); s++, l_it++)
          {
            // Get the coordinates
            zeta[0] = *l_it;
            boundary_pt->geom_object_pt()->position(zeta, posn);
            vertex_coord[s] = posn;

            // Bump up the polygonal arclength
            if (s > 0)
            {
              polygonal_vertex_arclength_info[s].first =
                polygonal_vertex_arclength_info[s - 1].first +
                sqrt(pow(vertex_coord[s][0] - vertex_coord[s - 1][0], 2) +
                     pow(vertex_coord[s][1] - vertex_coord[s - 1][1], 2));
              polygonal_vertex_arclength_info[s].second = zeta[0];
            }
          }
        } // Arclength uniformly spaced
      } // Less number of insertion points than vertices
    }

    /// \short Helper function that returns a polygon representation for
    /// the given closed curve, it also computes the maximum boundary id of
    /// the constituent curves.
    /// If the TriangleMeshClosedCurve is already a TriangleMeshPolygon
    /// we simply return a pointer to it. Otherwise a new TrilangleMeshPolygon
    /// is created -- this is deleted automatically when the TriangleMesh
    /// destructor is called, so no external book-keeping is required.
    TriangleMeshPolygon* closed_curve_to_polygon_helper(
      TriangleMeshClosedCurve* closed_curve_pt, unsigned& max_bnd_id_local)
    {
      // How many separate boundaries do we have
      unsigned nb = closed_curve_pt->ncurve_section();

#ifdef PARANOID
      if (nb < 2)
      {
        std::ostringstream error_message;
        error_message << "TriangleMeshClosedCurve that defines outer boundary\n"
                      << "must be made up of at least two "
                      << "TriangleMeshCurveSections\n"
                      << "to allow the automatic set up boundary coordinates.\n"
                      << "Yours only has (" << nb << ")" << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Provide storage for accompanying polylines
      Vector<TriangleMeshCurveSection*> my_boundary_polyline_pt(nb);

      // Store refinement tolerance
      Vector<double> refinement_tolerance(nb);

      // Store unrefinement tolerance
      Vector<double> unrefinement_tolerance(nb);

      // Store max. length
      Vector<double> max_length(nb);

      // Loop over boundaries that make up this boundary
      for (unsigned b = 0; b < nb; b++)
      {
        // Get pointer to the curve segment boundary that makes up
        // this part of the boundary
        TriangleMeshCurviLine* curviline_pt =
          dynamic_cast<TriangleMeshCurviLine*>(
            closed_curve_pt->curve_section_pt(b));

        TriangleMeshPolyLine* polyline_pt = dynamic_cast<TriangleMeshPolyLine*>(
          closed_curve_pt->curve_section_pt(b));

        if (curviline_pt != 0)
        {
          // Boundary id
          unsigned bnd_id = curviline_pt->boundary_id();

          // Build associated polyline
          my_boundary_polyline_pt[b] =
            curviline_to_polyline(curviline_pt, bnd_id);

          // Copy the unrefinement tolerance
          unrefinement_tolerance[b] = curviline_pt->unrefinement_tolerance();

          // Copy the refinement tolerance
          refinement_tolerance[b] = curviline_pt->refinement_tolerance();

          // Copy the maximum length
          max_length[b] = curviline_pt->maximum_length();

          // Updates bnd_id<--->curve section map
          Boundary_curve_section_pt[bnd_id] = my_boundary_polyline_pt[b];

          // Keep track of curve sections that need to be deleted!!!
          Free_curve_section_pt.insert(my_boundary_polyline_pt[b]);

          // Keep track...
          if (bnd_id > max_bnd_id_local)
          {
            max_bnd_id_local = bnd_id;
          }
        }
        else if (polyline_pt != 0)
        {
          // Boundary id
          unsigned bnd_id = polyline_pt->boundary_id();

          // Pass the pointer of the already existing polyline
          my_boundary_polyline_pt[b] = polyline_pt;

          // Copy the unrefinement tolerance
          unrefinement_tolerance[b] = polyline_pt->unrefinement_tolerance();

          // Copy the refinement tolerance
          refinement_tolerance[b] = polyline_pt->refinement_tolerance();

          // Copy the maximum length
          max_length[b] = polyline_pt->maximum_length();

          // Updates bnd_id<--->curve section map
          Boundary_curve_section_pt[bnd_id] = my_boundary_polyline_pt[b];

          // Keep track...
          if (bnd_id > max_bnd_id_local)
          {
            max_bnd_id_local = bnd_id;
          }
        }
        else
        {
          std::ostringstream error_stream;
          error_stream << "The 'curve_segment' is not a curviline neither a\n "
                       << "polyline: What is it?\n"
                       << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      } // end of loop over boundaries

      // Create a new polygon by using the new created polylines
      TriangleMeshPolygon* output_polygon_pt =
        new TriangleMeshPolygon(my_boundary_polyline_pt,
                                closed_curve_pt->internal_point(),
                                closed_curve_pt->is_internal_point_fixed());

      // Keep track of new created polygons that need to be deleted!!!
      Free_polygon_pt.insert(output_polygon_pt);

      // Pass on refinement information
      output_polygon_pt->set_polyline_refinement_tolerance(
        closed_curve_pt->polyline_refinement_tolerance());
      output_polygon_pt->set_polyline_unrefinement_tolerance(
        closed_curve_pt->polyline_unrefinement_tolerance());

      // Loop over boundaries that make up this boundary and copy
      // refinement, unrefinement and max length information
      for (unsigned b = 0; b < nb; b++)
      {
        // Set the unrefinement and refinement information
        my_boundary_polyline_pt[b]->set_unrefinement_tolerance(
          unrefinement_tolerance[b]);

        my_boundary_polyline_pt[b]->set_refinement_tolerance(
          refinement_tolerance[b]);

        // Copy the maximum length constraint
        my_boundary_polyline_pt[b]->set_maximum_length(max_length[b]);
      }
      return output_polygon_pt;
    }

    /// \short Helper function that creates and returns an open curve with
    /// the polyline representation of its constituent curve sections. The
    /// new created open curve is deleted when the TriangleMesh destructor
    /// is called
    TriangleMeshOpenCurve* create_open_curve_with_polyline_helper(
      TriangleMeshOpenCurve* open_curve_pt, unsigned& max_bnd_id_local)
    {
      unsigned nb = open_curve_pt->ncurve_section();

      // Provide storage for accompanying polylines
      Vector<TriangleMeshCurveSection*> my_boundary_polyline_pt(nb);

      // Store refinement tolerance
      Vector<double> refinement_tolerance(nb);

      // Store unrefinement tolerance
      Vector<double> unrefinement_tolerance(nb);

      // Store max. length
      Vector<double> max_length(nb);

      // Loop over the number of curve sections on the open curve
      for (unsigned i = 0; i < nb; i++)
      {
        // Get pointer to the curve segment boundary that makes up
        // this part of the boundary
        TriangleMeshCurviLine* curviline_pt =
          dynamic_cast<TriangleMeshCurviLine*>(
            open_curve_pt->curve_section_pt(i));
        TriangleMeshPolyLine* polyline_pt = dynamic_cast<TriangleMeshPolyLine*>(
          open_curve_pt->curve_section_pt(i));

        if (curviline_pt != 0)
        {
          // Boundary id
          unsigned bnd_id = curviline_pt->boundary_id();

          // Build associated polyline
          my_boundary_polyline_pt[i] =
            curviline_to_polyline(curviline_pt, bnd_id);

          // Copy the unrefinement tolerance
          unrefinement_tolerance[i] = curviline_pt->unrefinement_tolerance();

          // Copy the refinement tolerance
          refinement_tolerance[i] = curviline_pt->refinement_tolerance();

          // Copy the maximum length
          max_length[i] = curviline_pt->maximum_length();

          // Pass the connection information to the polyline representation
          copy_connection_information(curviline_pt, my_boundary_polyline_pt[i]);

          // Updates bnd_id<--->curve section map
          Boundary_curve_section_pt[bnd_id] = my_boundary_polyline_pt[i];

          // Keep track of curve sections that need to be deleted!!!
          Free_curve_section_pt.insert(my_boundary_polyline_pt[i]);

          // Keep track...
          if (bnd_id > max_bnd_id_local)
          {
            max_bnd_id_local = bnd_id;
          }
        }
        else if (polyline_pt != 0)
        {
          // Boundary id
          unsigned bnd_id = polyline_pt->boundary_id();

          // Storage pointer
          my_boundary_polyline_pt[i] = polyline_pt;

          // Copy the unrefinement tolerance
          unrefinement_tolerance[i] = polyline_pt->unrefinement_tolerance();

          // Copy the refinement tolerance
          refinement_tolerance[i] = polyline_pt->refinement_tolerance();

          // Copy the maximum length
          max_length[i] = polyline_pt->maximum_length();

          // Pass the connection information to the polyline representation
          copy_connection_information(polyline_pt, my_boundary_polyline_pt[i]);

          // Updates bnd_id<--->curve section map
          Boundary_curve_section_pt[bnd_id] = my_boundary_polyline_pt[i];

          // Keep track...
          if (bnd_id > max_bnd_id_local)
          {
            max_bnd_id_local = bnd_id;
          }
        }
        else
        {
          std::ostringstream error_stream;
          error_stream
            << "The 'curve_segment' (open) is not a curviline neither a\n "
            << "polyline: What is it?\n"
            << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      } // end of loop over boundaries

      // Create open curve with polylines boundaries
      TriangleMeshOpenCurve* output_open_polyline_pt =
        new TriangleMeshOpenCurve(my_boundary_polyline_pt);

      // Keep track of open polylines that need to be deleted!!!
      Free_open_curve_pt.insert(output_open_polyline_pt);

      // Pass on refinement information
      output_open_polyline_pt->set_polyline_refinement_tolerance(
        open_curve_pt->polyline_refinement_tolerance());
      output_open_polyline_pt->set_polyline_unrefinement_tolerance(
        open_curve_pt->polyline_unrefinement_tolerance());

      // Loop over boundaries that make up this boundary and copy
      // refinement, unrefinement and max length information
      for (unsigned b = 0; b < nb; b++)
      {
        // Set the unrefinement and refinement information
        my_boundary_polyline_pt[b]->set_unrefinement_tolerance(
          unrefinement_tolerance[b]);

        my_boundary_polyline_pt[b]->set_refinement_tolerance(
          refinement_tolerance[b]);

        // Copy the maximum length constraint
        my_boundary_polyline_pt[b]->set_maximum_length(max_length[b]);
      }
      return output_open_polyline_pt;
    }

    /// \short Stores the geometric objects associated to the
    /// curve sections that compound the closed curve. It also
    /// stores the limits defined by these geometric objects
    void set_geom_objects_and_coordinate_limits_for_close_curve(
      TriangleMeshClosedCurve* input_closed_curve_pt)
    {
      unsigned nb = input_closed_curve_pt->ncurve_section();

#ifdef PARANOID

      if (nb < 2)
      {
        std::ostringstream error_message;
        error_message << "TriangleMeshCurve that defines closed boundary\n"
                      << "must be made up of at least two "
                      << "TriangleMeshCurveSection\n"
                      << "to allow the automatic set up boundary coordinates.\n"
                      << "Yours only has " << nb << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif

      // TODO: Used for the ImmersedRigidBodyTriangleMeshPolygon objects only
      // ImmersedRigidBodyTriangleMeshPolygon* bound_geom_obj_pt
      //= dynamic_cast<ImmersedRigidBodyTriangleMeshPolygon*>
      // (input_closed_curve_pt);
      GeomObject* bound_geom_obj_pt =
        dynamic_cast<GeomObject*>(input_closed_curve_pt);

      // If cast successful set up the coordinates
      if (bound_geom_obj_pt != 0)
      {
        unsigned n_poly = input_closed_curve_pt->ncurve_section();
        for (unsigned p = 0; p < n_poly; p++)
        {
          // Read out the index of the boundary from the polyline
          unsigned b_index =
            input_closed_curve_pt->curve_section_pt(p)->boundary_id();

          // Set the geometric object
          Boundary_geom_object_pt[b_index] = bound_geom_obj_pt;

          // The coordinates along each polyline boundary are scaled to
          // of unit length so the total coordinate limits are simply
          // (p,p+1)
          Boundary_coordinate_limits[b_index].resize(2);
          Boundary_coordinate_limits[b_index][0] = p;
          Boundary_coordinate_limits[b_index][1] = p + 1.0;
        }

#ifdef PARANOID
        // If we are using parallel mesh adaptation and load balancing,
        // we are going to need to check for the use of this type of
        // Polygon at this stage, so switch on the flag
        Immersed_rigid_body_triangle_mesh_polygon_used = true;
#endif
      }
      else
      {
        // Loop over curve sections that make up this boundary
        for (unsigned b = 0; b < nb; b++)
        {
          TriangleMeshCurviLine* curviline_pt =
            dynamic_cast<TriangleMeshCurviLine*>(
              input_closed_curve_pt->curve_section_pt(b));

          if (curviline_pt != 0)
          {
            // Read the values of the limiting coordinates
            Vector<double> zeta_bound(2);
            zeta_bound[0] = curviline_pt->zeta_start();
            zeta_bound[1] = curviline_pt->zeta_end();

            // Boundary id
            unsigned bnd_id = curviline_pt->boundary_id();

            // Set the boundary geometric object and limits
            Boundary_geom_object_pt[bnd_id] = curviline_pt->geom_object_pt();
            Boundary_coordinate_limits[bnd_id] = zeta_bound;
          }
        } // for
      } // else
    } // function

    /// \short Stores the geometric objects associated to the
    /// curve sections that compound the open curve. It also
    /// stores the limits defined by these geometric objects
    void set_geom_objects_and_coordinate_limits_for_open_curve(
      TriangleMeshOpenCurve* input_open_curve_pt)
    {
      unsigned nb = input_open_curve_pt->ncurve_section();

      // Loop over curve sections that make up this boundary
      for (unsigned b = 0; b < nb; b++)
      {
        TriangleMeshCurviLine* curviline_pt =
          dynamic_cast<TriangleMeshCurviLine*>(
            input_open_curve_pt->curve_section_pt(b));

        if (curviline_pt != 0)
        {
          // ead the values of the limiting coordinates
          Vector<double> zeta_bound(2);
          zeta_bound[0] = curviline_pt->zeta_start();
          zeta_bound[1] = curviline_pt->zeta_end();

          // Boundary id
          unsigned bnd_id = curviline_pt->boundary_id();

          // Set the boundary geometric object and limits
          Boundary_geom_object_pt[bnd_id] = curviline_pt->geom_object_pt();
          Boundary_coordinate_limits[bnd_id] = zeta_bound;
        }
      } // for
    } // function

#endif // OOMPH_HAS_TRIANGLE_LIB

  private:
#ifdef OOMPH_HAS_TRIANGLE_LIB

    /// \short Helper function that creates the associated polyline
    /// representation for curvilines
    TriangleMeshCurveSection* curviline_to_polyline(
      TriangleMeshCurviLine*& curviline_pt, unsigned& bnd_id)
    {
      // Create vertex coordinates for polygonal representation
      Vector<Vector<double>> bound;
      Vector<std::pair<double, double>> polygonal_vertex_arclength;

      if (curviline_pt->are_there_connection_points())
      {
        this->create_vertex_coordinates_for_polyline_connections(
          curviline_pt, bound, polygonal_vertex_arclength);
      }
      else
      {
        this->create_vertex_coordinates_for_polyline_no_connections(
          curviline_pt, bound, polygonal_vertex_arclength);
      }

      // Store the vertex-arclength information
      Polygonal_vertex_arclength_info[bnd_id] = polygonal_vertex_arclength;

      // Build associated polyline
      return new TriangleMeshPolyLine(bound, bnd_id);
    }

    /// \short Get the associated vertex to the given s value by looking to
    /// the list of s values created when changing from curviline to polyline
    unsigned get_associated_vertex_to_svalue(double& target_s_value,
                                             unsigned& bnd_id)
    {
      double s_tolerance = 1.0e-14;
      return get_associated_vertex_to_svalue(
        target_s_value, bnd_id, s_tolerance);
    }

    /// \short Get the associated vertex to the given s value by looking to
    /// the list of s values created when changing from curviline to polyline
    unsigned get_associated_vertex_to_svalue(double& target_s_value,
                                             unsigned& bnd_id,
                                             double& s_tolerance)
    {
      // Create a pointer to the list of s coordinates and arclength values
      // associated with a vertex
      Vector<std::pair<double, double>>* vertex_info =
        &Polygonal_vertex_arclength_info[bnd_id];

      // Total vertex number
      unsigned vector_size = vertex_info->size();

      // Counter for current vertex number
      unsigned n_vertex = 0;

      // Find the associated value to the given s value
      do
      {
        // Store the current zeta value
        double s = (*vertex_info)[n_vertex].second;

        // When find it save the vertex number and return it
        if (std::fabs(s - target_s_value) < s_tolerance)
        {
          break;
        }

        // Increment n_vertex
        n_vertex++;
      } while (n_vertex < vector_size);

#ifdef PARANOID

      if (n_vertex >= vector_size)
      {
        std::ostringstream error_message;
        error_message << "Could not find the associated vertex number in\n"
                      << "boundary " << bnd_id << " with the given s\n"
                      << "connection value (" << target_s_value << ") using\n"
                      << "this tolerance: " << s_tolerance << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

#endif
      return n_vertex;
    }

#endif // OOMPH_HAS_TRIANGLE_LIB
  };


  //======================================================================
  /// Setup boundary coordinate on boundary b. Doc Faces
  /// in outfile. Boundary coordinate increases continously along
  /// polygonal boundary. It's zero at the lexicographically
  /// smallest node on the boundary.
  //======================================================================
  template<class ELEMENT>
  void UnstructuredTwoDMeshGeometryBase::setup_boundary_coordinates(
    const unsigned& b, std::ofstream& outfile)
  {
    // Temporary storage for face elements
    Vector<FiniteElement*> face_el_pt;

    // Temporary storage for number of elements adjacent to the boundary
    unsigned nel = 0;

    // =================================================================
    // BEGIN: Get face elements from boundary elements
    // =================================================================

    // Temporary storage for elements adjacent to the boundary that have
    // an common edge (related with internal boundaries)
    unsigned n_repeated_ele = 0;

    unsigned n_regions = this->nregion();

#ifdef OOMPH_HAS_MPI
    // map to associate the face element to the bulk element, this info.
    // is only necessary for the setup of boundary coordinates in a
    // distributed mesh when we need to extract the halo/haloed info.
    std::map<FiniteElement*, FiniteElement*> face_to_bulk_element_pt;
#endif

    // Temporary storage for already done nodes
    Vector<std::pair<Node*, Node*>> done_nodes_pt;

    // If there is more than one region then only use boundary
    // coordinates from the bulk side (region 0)
    if (n_regions > 1)
    {
      for (unsigned rr = 0; rr < n_regions; rr++)
      {
        unsigned region_id = static_cast<unsigned>(this->Region_attribute[rr]);

#ifdef PARANOID
        double diff =
          fabs(Region_attribute[rr] - static_cast<double>(static_cast<unsigned>(
                                        this->Region_attribute[rr])));
        if (diff > 0.0)
        {
          std::ostringstream error_message;
          error_message << "Region attributes should be unsigneds because we \n"
                        << "only use them to set region ids\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Loop over all elements on boundaries in region rr
        unsigned nel_in_region =
          this->nboundary_element_in_region(b, region_id);
        unsigned nel_repeated_in_region = 0;

#ifdef PARANOID
        if (!Suppress_warning_about_regions_and_boundaries)
        {
          if (nel_in_region == 0)
          {
            std::ostringstream warning_message;
            std::string output_string =
              "UnstructuredTwoDMeshGeometryBase::setup_boundary_coordinates()";
            warning_message
              << "There are no elements associated with boundary (" << b
              << ")\n"
              << "in region (" << region_id << "). This could happen because:\n"
              << "1) You did not specify boundaries with this boundary id.\n"
              << "---- Review carefully the indexing of your boundaries.\n"
              << "2) The boundary (" << b << ") is not associated with region ("
              << region_id << ").\n"
              << "---- The boundary does not touch the region.\n"
              << "You can suppress this warning by setting the static public "
                 "bool\n\n"
              << "   "
                 "UnstructuredTwoDMeshGeometryBase::Suppress_warning_about_"
                 "regions_and_boundaries\n\n"
              << "to true.\n";
            OomphLibWarning(
              warning_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
          }
        }
#endif

        // Only bother to do anything else, if there are elements
        // associated with the boundary and the current region
        if (nel_in_region > 0)
        {
          // Flag that activates when a repeated face element is found,
          // possibly because we are dealing with an internal boundary
          bool repeated = false;

          // Loop over the bulk elements adjacent to boundary b
          for (unsigned e = 0; e < nel_in_region; e++)
          {
            // Get pointer to the bulk element that is adjacent to boundary b
            FiniteElement* bulk_elem_pt =
              this->boundary_element_in_region_pt(b, region_id, e);

#ifdef OOMPH_HAS_MPI
            // In a distributed mesh only work with nonhalo elements
            if (this->is_mesh_distributed() && bulk_elem_pt->is_halo())
            {
              // Increase the number of repeated elements
              n_repeated_ele++;
              // Skip this element and go for the next one
              continue;
            }
#endif

            // Find the index of the face of element e along boundary b
            int face_index =
              this->face_index_at_boundary_in_region(b, region_id, e);

            // Before adding the new element we need to be sure that
            // the edge that this element represent has not been
            // already added
            FiniteElement* tmp_ele_pt =
              new DummyFaceElement<ELEMENT>(bulk_elem_pt, face_index);

            const unsigned n_nodes = tmp_ele_pt->nnode();

            std::pair<Node*, Node*> tmp_pair = std::make_pair(
              tmp_ele_pt->node_pt(0), tmp_ele_pt->node_pt(n_nodes - 1));

            std::pair<Node*, Node*> tmp_pair_inverse = std::make_pair(
              tmp_ele_pt->node_pt(n_nodes - 1), tmp_ele_pt->node_pt(0));

            // Search for repeated nodes
            const unsigned n_done_nodes = done_nodes_pt.size();
            for (unsigned l = 0; l < n_done_nodes; l++)
            {
              if (tmp_pair == done_nodes_pt[l] ||
                  tmp_pair_inverse == done_nodes_pt[l])
              {
                nel_repeated_in_region++;
                repeated = true;
                break;
              }
            }

            // Create new face element?
            if (!repeated)
            {
              // Add the pair of nodes (edge) to the node dones
              done_nodes_pt.push_back(tmp_pair);
#ifdef OOMPH_HAS_MPI
              // If the mesh is distributed then create a map from the
              // temporary face element to the bulk element, further
              // info. will be extracted from the bulk element for setup
              // of boundary coordinates in a distributed mesh
              if (this->is_mesh_distributed())
              {
                face_to_bulk_element_pt[tmp_ele_pt] = bulk_elem_pt;
              }
#endif
              // Add the face element to the storage
              face_el_pt.push_back(tmp_ele_pt);
            }
            else
            {
              // Clean up
              delete tmp_ele_pt;
              tmp_ele_pt = 0;
            }

            // Re-start
            repeated = false;

            // Output faces?
            if (outfile.is_open())
            {
              face_el_pt[face_el_pt.size() - 1]->output(outfile);
            }
          } // for nel

          nel += nel_in_region;

          n_repeated_ele += nel_repeated_in_region;

        } // if (nel_in_region > 0)

      } // for (rr < n_regions)

    } // if (n_regions > 1)
    // Otherwise it's just the normal boundary functions
    else
    {
      // Loop over all elements on boundaries
      nel = this->nboundary_element(b);

#ifdef PARANOID
      if (!Suppress_warning_about_regions_and_boundaries)
      {
        if (nel == 0)
        {
          std::ostringstream warning_message;
          std::string output_string =
            "UnstructuredTwoDMeshGeometryBase::setup_boundary_coordinates()";
          warning_message
            << "There are no elements associated with boundary (" << b << ").\n"
            << "This could happen because you did not specify boundaries with\n"
            << "this boundary id. Review carefully the indexing of your\n"
            << "boundaries.";
          OomphLibWarning(
            warning_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
        }
      }
#endif

      // Only bother to do anything else, if there are elements
      if (nel > 0)
      {
        // Flag that activates when a repeated face element is found,
        // possibly because we are dealing with an internal boundary
        bool repeated = false;

        // Loop over the bulk elements adjacent to boundary b
        for (unsigned e = 0; e < nel; e++)
        {
          // Get pointer to the bulk element that is adjacent to boundary b
          FiniteElement* bulk_elem_pt = this->boundary_element_pt(b, e);

#ifdef OOMPH_HAS_MPI

          // In a distributed mesh only work with nonhalo elements
          if (this->is_mesh_distributed() && bulk_elem_pt->is_halo())
          {
            // Increase the number of repeated elements
            n_repeated_ele++;
            // Skip this element and go for the next one
            continue;
          }

#endif

          // Find the index of the face of element e along boundary b
          int face_index = this->face_index_at_boundary(b, e);

          // Before adding the new element we need to be sure that the
          // edge that this element represent has not been already added
          FiniteElement* tmp_ele_pt =
            new DummyFaceElement<ELEMENT>(bulk_elem_pt, face_index);

          const unsigned n_nodes = tmp_ele_pt->nnode();

          std::pair<Node*, Node*> tmp_pair = std::make_pair(
            tmp_ele_pt->node_pt(0), tmp_ele_pt->node_pt(n_nodes - 1));

          std::pair<Node*, Node*> tmp_pair_inverse = std::make_pair(
            tmp_ele_pt->node_pt(n_nodes - 1), tmp_ele_pt->node_pt(0));

          // Search for repeated nodes
          const unsigned n_done_nodes = done_nodes_pt.size();
          for (unsigned l = 0; l < n_done_nodes; l++)
          {
            if (tmp_pair == done_nodes_pt[l] ||
                tmp_pair_inverse == done_nodes_pt[l])
            {
              n_repeated_ele++;
              repeated = true;
              break;
            }
          }

          // Create new face element
          if (!repeated)
          {
            // Add the pair of nodes (edge) to the node dones
            done_nodes_pt.push_back(tmp_pair);
#ifdef OOMPH_HAS_MPI
            // Create a map from the temporary face element to the bulk
            // element, further info. will be extracted from the bulk
            // element for setup of boundary coordinates in a
            // distributed mesh
            if (this->is_mesh_distributed())
            {
              face_to_bulk_element_pt[tmp_ele_pt] = bulk_elem_pt;
            }
#endif
            face_el_pt.push_back(tmp_ele_pt);
          }
          else
          {
            // Free the repeated bulk element!!
            delete tmp_ele_pt;
            tmp_ele_pt = 0;
          }

          // Re-start
          repeated = false;

          // Output faces?
          if (outfile.is_open())
          {
            face_el_pt[face_el_pt.size() - 1]->output(outfile);
          }

        } // for (e < nel)

      } // if (nel > 0)

    } // else if (n_regions > 1)

    // Do not consider the repeated elements
    nel -= n_repeated_ele;

#ifdef PARANOID
    if (nel != face_el_pt.size())
    {
      std::ostringstream error_message;
      error_message
        << "The independent counting of face elements (" << nel << ") for "
        << "boundary (" << b << ") is different\n"
        << "from the real number of face elements in the container ("
        << face_el_pt.size() << ")\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // =================================================================
    // END: Get face elements from boundary elements
    // =================================================================

    // Only bother to do anything else, if there are elements
    if (nel > 0)
    {
      // A flag vector to mark those face elements that are considered
      // as halo in the current processor
      std::vector<bool> is_halo_face_element(nel, false);

      // Count the total number of non halo face elements
      unsigned nnon_halo_face_elements = 0;

#ifdef OOMPH_HAS_MPI
      // Only mark the face elements as halo if the mesh is marked as
      // distributed
      if (this->is_mesh_distributed())
      {
        for (unsigned ie = 0; ie < nel; ie++)
        {
          FiniteElement* face_ele_pt = face_el_pt[ie];
          // Get the bulk element
          FiniteElement* tmp_bulk_ele_pt = face_to_bulk_element_pt[face_ele_pt];
          // Check if the bulk element is halo
          if (!tmp_bulk_ele_pt->is_halo())
          {
            // Mark the face element as nonhalo
            is_halo_face_element[ie] = false;
            // Increase the counter for nonhalo elements
            nnon_halo_face_elements++;
          }
          else
          {
            // Mark the face element as halo
            is_halo_face_element[ie] = true;
          }
        } // for (ie < nel)
      } // if (this->is_mesh_distributed())
      else
      {
#endif // OOMPH_HAS_MPI

        // If the mesh is not distributed then the number of non halo
        // elements is the same as the number of elements
        nnon_halo_face_elements = nel;

#ifdef OOMPH_HAS_MPI
      } // else if (this->is_mesh_distributed())
#endif

#ifdef PARANOID
      // Get the total number of halo face elements
      const unsigned nhalo_face_element = nel - nnon_halo_face_elements;

      if (nhalo_face_element > 0)
      {
        std::ostringstream error_message;
        error_message
          << "There should not be halo face elements since they were not\n"
          << "considered when computing the face elements.\n"
          << "The number of found halo face elements is: ("
          << nhalo_face_element << ").\n\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // =================================================================
      // BEGIN: Sort face elements
      // =================================================================

      // The vector of lists to store the "segments" that compound the
      // boundaries (segments may appear only in a distributed mesh
      // because the boundary may have been split across multiple
      // processors)
      Vector<std::list<FiniteElement*>> segment_sorted_ele_pt;

      // Number of already sorted face elements (only nonhalo face
      // elements for a distributed mesh)
      unsigned nsorted_face_elements = 0;

      // Keep track of who's done (in a distributed mesh this apply to
      // nonhalo only)
      std::map<FiniteElement*, bool> done_el;

      // Keep track of which element is inverted (in distributed mesh
      // the elements may be inverted with respect to the segment they
      // belong)
      std::map<FiniteElement*, bool> is_inverted;

      // Iterate until all possible segments have been created. In a
      // non distributed mesh there is only one segment which defines
      // the complete boundary
      while (nsorted_face_elements < nnon_halo_face_elements)
      {
        // The sorted list of face elements (in a distributed mesh a
        // collection of continuous face elements define a segment)
        std::list<FiniteElement*> sorted_el_pt;

        FiniteElement* ele_face_pt = 0;

#ifdef PARANOID
        // Select an initial element for the segment
        bool found_initial_face_element = false;
#endif

        // Store the index of the initial face element
        unsigned iface = 0;
#ifdef OOMPH_HAS_MPI
        if (this->is_mesh_distributed())
        {
          for (iface = 0; iface < nel; iface++)
          {
            // Only work with nonhalo face elements
            if (!is_halo_face_element[iface])
            {
              ele_face_pt = face_el_pt[iface];
              // If not done then take it as initial face element
              if (!done_el[ele_face_pt])
              {
#ifdef PARANOID
                // Set the flag to indicate the initial element was
                // found
                found_initial_face_element = true;
#endif
                // Increase the number of sorted face elements
                nsorted_face_elements++;
                // Set the index to the next face element
                iface++;
                // Add the face element in the container
                sorted_el_pt.push_back(ele_face_pt);
                // Mark as done
                done_el[ele_face_pt] = true;
                break;
              } // if (!done_el[ele_face_pt])
            } // if (!is_halo_face_element[iface])
          } // for (iface < nel)
        } // if (this->is_mesh_distributed())
        else
        {
#endif // #ifdef OOMPH_HAS_MPI

          // When the mesh is not distributed just take the first
          // element and put it in the sorted list
          ele_face_pt = face_el_pt[0];
#ifdef PARANOID
          // Set the flag to indicate the initial element was found
          found_initial_face_element = true;
#endif
          // Increase the number of sorted face elements
          nsorted_face_elements++;
          // Set the index to the next face element
          iface = 1;
          // Add the face element in the container
          sorted_el_pt.push_back(ele_face_pt);
          // Mark as done
          done_el[ele_face_pt] = true;

#ifdef OOMPH_HAS_MPI
        } // else if (this->is_mesh_distributed())
#endif

#ifdef PARANOID
        if (!found_initial_face_element)
        {
          std::ostringstream error_message;
          error_message << "Could not find an initial face element for the "
                           "current segment\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Number of nodes of the initial face element
        const unsigned nnod = ele_face_pt->nnode();

        // Left and rightmost nodes (the left and right nodes of the
        // current face element)
        Node* left_node_pt = ele_face_pt->node_pt(0);
        Node* right_node_pt = ele_face_pt->node_pt(nnod - 1);

        // Continue iterating if a new face element has been added to
        // the list
        bool face_element_added = false;

        // While a new face element has been added to the set of sorted
        // face elements continue iterating
        do
        {
          // Start from the next face element since we have already
          // added the previous one as the initial face element (any
          // previous face element had to be added on previous
          // iterations)
          for (unsigned iiface = iface; iiface < nel; iiface++)
          {
            // Re-start flag
            face_element_added = false;

            // Get the candidate element
            ele_face_pt = face_el_pt[iiface];

            // Check that the candidate element has not been done and
            // is not a halo element
            if (!(done_el[ele_face_pt] || is_halo_face_element[iiface]))
            {
              // Get the left and right nodes of the current element
              Node* local_left_node_pt = ele_face_pt->node_pt(0);
              Node* local_right_node_pt = ele_face_pt->node_pt(nnod - 1);

              // New element fits at the left of segment and is not inverted
              if (left_node_pt == local_right_node_pt)
              {
                left_node_pt = local_left_node_pt;
                sorted_el_pt.push_front(ele_face_pt);
                is_inverted[ele_face_pt] = false;
                face_element_added = true;
              }
              // New element fits at the left of segment and is inverted
              else if (left_node_pt == local_left_node_pt)
              {
                left_node_pt = local_right_node_pt;
                sorted_el_pt.push_front(ele_face_pt);
                is_inverted[ele_face_pt] = true;
                face_element_added = true;
              }
              // New element fits on the right of segment and is not inverted
              else if (right_node_pt == local_left_node_pt)
              {
                right_node_pt = local_right_node_pt;
                sorted_el_pt.push_back(ele_face_pt);
                is_inverted[ele_face_pt] = false;
                face_element_added = true;
              }
              // New element fits on the right of segment and is inverted
              else if (right_node_pt == local_right_node_pt)
              {
                right_node_pt = local_left_node_pt;
                sorted_el_pt.push_back(ele_face_pt);
                is_inverted[ele_face_pt] = true;
                face_element_added = true;
              }

              if (face_element_added)
              {
                done_el[ele_face_pt] = true;
                nsorted_face_elements++;
                break;
              }

            } // if (!(done_el[ele_face_pt] || is_halo_face_element[iiface]))
          } // for (iiface<nnon_halo_face_element)
        } while (face_element_added &&
                 (nsorted_face_elements < nnon_halo_face_elements));

        // Store the created segment in the vector of segments
        segment_sorted_ele_pt.push_back(sorted_el_pt);

      } // while(nsorted_face_elements < nnon_halo_face_elements);

#ifdef OOMPH_HAS_MPI
      if (!this->is_mesh_distributed())
      {
#endif
        // Are we done?
        if (nsorted_face_elements != nel || segment_sorted_ele_pt.size() != 1)
        {
          std::ostringstream error_message;
          error_message << "Was only able to setup boundary coordinate on "
                        << "boundary " << b << "\nfor " << nsorted_face_elements
                        << " of " << nel
                        << " face elements. This usually means\n"
                        << "that the boundary is not simply connected.\n\n"
                        << "Re-run the setup_boundary_coordintes() function\n"
                        << "with an output file specified "
                        << "as the second argument.\n"
                        << "This file will contain FaceElements that\n"
                        << "oomph-lib believes to be located on the boundary.\n"
                        << std::endl;
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#ifdef OOMPH_HAS_MPI
      } // if (!this->is_mesh_distributed())
#endif

      // =================================================================
      // END: Sort face elements
      // =================================================================

      // ----------------------------------------------------------------

      // =================================================================
      // BEGIN: Assign global/local (non distributed mesh/distributed
      // mesh) boundary coordinates to nodes
      // =================================================================

      // Compute the (local) boundary coordinates of the nodes in the
      // segments. In a distributed mesh this info. will be used by a
      // root processor to compute the (global) boundary coordinates.

      // Vector of sets that stores the nodes of each segment based on
      // a lexicographically order starting from the bottom left node
      // of each segment
      Vector<std::set<Node*>> segment_all_nodes_pt;

      // The number of segments in this processor
      const unsigned nsegments = segment_sorted_ele_pt.size();

#ifdef PARANOID
      if (nnon_halo_face_elements > 0 && nsegments == 0)
      {
        std::ostringstream error_message;
        error_message
          << "The number of segments is zero, but the number of nonhalo\n"
          << "elements is: (" << nnon_halo_face_elements << ")\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      } // if (nnon_halo_face_elements > 0 && nsegments == 0)

#endif

      // Store the arclength of each segment in the current processor
      Vector<double> segment_arclength(nsegments);

#ifdef OOMPH_HAS_MPI

      // Clear the boundary segment nodes storage
      this->flush_boundary_segment_node(b);

      // Set the number of segments for the current boundary
      this->set_nboundary_segment_node(b, nsegments);

#endif // #ifdef OOMPH_HAS_MPI

      // Go through all the segments and compute their (local) boundary
      // coordinates
      for (unsigned is = 0; is < nsegments; is++)
      {
#ifdef PARANOID

        if (segment_sorted_ele_pt[is].size() == 0)
        {
          std::ostringstream error_message;
          std::string output_string =
            "UnstructuredTwoDMeshGeometryBase::setup_boundary_coordinates()";
          error_message << "The (" << is << ")-th segment has no elements\n";
          throw OomphLibError(
            error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
        } // if (segment_sorted_ele_pt[is].size() == 0)

#endif

        // Get access to the first element on the segment
        FiniteElement* first_ele_pt = segment_sorted_ele_pt[is].front();

        // Number of nodes
        unsigned nnod = first_ele_pt->nnode();

        // Get the first node of the current segment
        Node* first_node_pt = first_ele_pt->node_pt(0);
        if (is_inverted[first_ele_pt])
        {
          first_node_pt = first_ele_pt->node_pt(nnod - 1);
        }

        // Coordinates of left node
        double x_left = first_node_pt->x(0);
        double y_left = first_node_pt->x(1);

        // Initialise boundary coordinate (local boundary coordinate
        // for boundaries with more than one segment)
        Vector<double> zeta(1, 0.0);

        // Set boundary coordinate
        first_node_pt->set_coordinates_on_boundary(b, zeta);

        // Lexicographically bottom left node
        std::set<Node*> local_nodes_pt;
        local_nodes_pt.insert(first_node_pt);

#ifdef OOMPH_HAS_MPI

        // Insert the node in the look-up scheme for nodes in segments
        this->add_boundary_segment_node(b, is, first_node_pt);

#endif // #ifdef OOMPH_HAS_MPI

        // Now loop over nodes in order
        for (std::list<FiniteElement*>::iterator it =
               segment_sorted_ele_pt[is].begin();
             it != segment_sorted_ele_pt[is].end();
             it++)
        {
          // Get element
          FiniteElement* el_pt = *it;

          // Start node and increment
          unsigned k_nod = 1;
          int nod_diff = 1;
          if (is_inverted[el_pt])
          {
            k_nod = nnod - 2;
            nod_diff = -1;
          }

          // Loop over nodes
          for (unsigned j = 1; j < nnod; j++)
          {
            Node* nod_pt = el_pt->node_pt(k_nod);
            k_nod += nod_diff;

            // Coordinates of right node
            double x_right = nod_pt->x(0);
            double y_right = nod_pt->x(1);

            // Increment boundary coordinate
            zeta[0] += sqrt((x_right - x_left) * (x_right - x_left) +
                            (y_right - y_left) * (y_right - y_left));

            // Set boundary coordinate
            nod_pt->set_coordinates_on_boundary(b, zeta);

            // Increment reference coordinate
            x_left = x_right;
            y_left = y_right;

            // Get lexicographically bottom left node but only
            // use vertex nodes as candidates
            local_nodes_pt.insert(nod_pt);

#ifdef OOMPH_HAS_MPI

            // Insert the node in the look-up scheme for nodes in segments
            this->add_boundary_segment_node(b, is, nod_pt);

#endif // #ifdef OOMPH_HAS_MPI

          } // for (j < nnod)
        } // iterator over the elements in the segment

        // Assign the arclength of the current segment
        segment_arclength[is] = zeta[0];

        // Add the nodes for the corresponding segment in the container
        segment_all_nodes_pt.push_back(local_nodes_pt);

      } // for (is < nsegments)

      // =================================================================
      // END: Assign global/local (non distributed mesh/distributed
      // mesh) boundary coordinates to nodes
      // =================================================================

      // Store the initial arclength for each segment of boundary in
      // the current processor, initialise to zero in case we have a
      // non distributed mesh
      Vector<double> initial_segment_arclength(nsegments, 0.0);

      // Store the final arclength for each segment of boundary in the
      // current processor, initalise to zero in case we have a non
      // distributed mesh
      Vector<double> final_segment_arclength(nsegments, 0.0);

      // Store the initial zeta for each segment of boundary in the
      // current processor, initalise to zero in case we have a non
      // distributed mesh
      Vector<double> initial_segment_zeta(nsegments, 0.0);

      // Store the final zeta for each segment of boundary in the
      // current processor, initalise to zero in case we have a non
      // distributed mesh
      Vector<double> final_segment_zeta(nsegments, 0.0);

      // --------------------------------------------------------------
      // DISTRIBUTED MESH: BEGIN - Check orientation of boundaries and
      // assign boundary coordinates accordingly
      // --------------------------------------------------------------

#ifdef OOMPH_HAS_MPI

      // Check that the mesh is distributed and that the initial zeta
      // values for the boundary segments have been already
      // assigned. When the method is called by the first time (when
      // the mesh is just created) the initial zeta values for the
      // boundary segments are not available
      if (this->is_mesh_distributed() &&
          Assigned_segments_initial_zeta_values[b])
      {
        // Get the initial and final coordinates of each segment

        // For each segment in the processor check whether it was
        // inverted in the root processor, if that is the case then it
        // is necessary to re-sort the face elements and invert them
        for (unsigned is = 0; is < nsegments; is++)
        {
          // Check if we need/or not to invert the current zeta values
          // on the segment (only apply for boundaries with GeomObject
          // associated)

          // Does the boundary HAS a GeomObject associated?
          if (this->boundary_geom_object_pt(b) != 0)
          {
            // Get the first and last node of the current segment and
            // their zeta values

            // Get access to the first element on the segment
            FiniteElement* first_ele_pt = segment_sorted_ele_pt[is].front();

            // Number of nodes
            const unsigned nnod = first_ele_pt->nnode();

            // Get the first node of the current segment
            Node* first_node_pt = first_ele_pt->node_pt(0);
            if (is_inverted[first_ele_pt])
            {
              first_node_pt = first_ele_pt->node_pt(nnod - 1);
            }

            // Get access to the last element on the segment
            FiniteElement* last_ele_pt = segment_sorted_ele_pt[is].back();

            // Get the last node of the current segment
            Node* last_node_pt = last_ele_pt->node_pt(nnod - 1);
            if (is_inverted[last_ele_pt])
            {
              last_node_pt = last_ele_pt->node_pt(0);
            }

            // Get the zeta coordinates for the first and last node
            Vector<double> current_segment_initial_zeta(1);
            Vector<double> current_segment_final_zeta(1);

            first_node_pt->get_coordinates_on_boundary(
              b, current_segment_initial_zeta);
            last_node_pt->get_coordinates_on_boundary(
              b, current_segment_final_zeta);

#ifdef PARANOID
            // Check whether the zeta values in the segment are in
            // increasing or decreasing order
            if (current_segment_initial_zeta[0] < current_segment_final_zeta[0])
            {
              // do nothing
            }
            else if (current_segment_initial_zeta[0] >
                     current_segment_final_zeta[0])
            {
              std::stringstream error_message;
              std::string output_string = "UnstructuredTwoDMeshGeometryBase::"
                                          "setup_boundary_coordinates()";
              error_message
                << "The zeta values are in decreasing order, this is weird\n"
                << "since they have just been set-up in increasing order few\n"
                << "lines above\n"
                << "Boundary: (" << b << ")\n"
                << "Segment: (" << is << ")\n"
                << "Initial zeta value: (" << current_segment_initial_zeta[0]
                << ")\n"
                << "Initial coordinate: (" << first_node_pt->x(0) << ", "
                << first_node_pt->x(1) << ")\n"
                << "Final zeta value: (" << current_segment_final_zeta[0]
                << ")\n"
                << "Initial coordinate: (" << last_node_pt->x(0) << ", "
                << last_node_pt->x(1) << ")\n";
              throw OomphLibError(
                error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
            }
            else
            {
              std::stringstream error_message;
              std::string output_string = "UnstructuredTwoDMeshGeometryBase::"
                                          "setup_boundary_coordinates()";
              error_message
                << "It was not possible to determine whether the zeta values on"
                << "the current segment\nof the boundary are in"
                << "increasing or decreasing order\n\n"
                << "Boundary: (" << b << ")\n"
                << "Segment: (" << is << ")\n"
                << "Arclength: (" << segment_arclength[is] << "\n"
                << "Initial zeta value: (" << current_segment_initial_zeta[0]
                << ")\n"
                << "Initial coordinate: (" << first_node_pt->x(0) << ", "
                << first_node_pt->x(1) << ")\n"
                << "Final zeta value: (" << current_segment_final_zeta[0]
                << ")\n"
                << "Initial coordinate: (" << last_node_pt->x(0) << ", "
                << last_node_pt->x(1) << ")\n";
              throw OomphLibError(
                error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
            }
#endif // #ifdef PARANOID

            // Now get the original zeta values and check if they are in
            // increasing or decreasing order
            const double original_segment_initial_zeta =
              boundary_segment_initial_zeta(b)[is];
            const double original_segment_final_zeta =
              boundary_segment_final_zeta(b)[is];

            // Now check if the zeta values go in increase or decrease
            // order
            if (original_segment_final_zeta > original_segment_initial_zeta)
            {
              // The original values go in increasing order, only need
              // to change the values if the original segment is marked
              // as inverted
              if (boundary_segment_inverted(b)[is])
              {
                // The original segment is inverted, then we need to
                // reverse the boundary coordinates. Go through all the
                // nodes and reverse its value
                std::set<Node*> all_nodes_pt = segment_all_nodes_pt[is];
                for (std::set<Node*>::iterator it = all_nodes_pt.begin();
                     it != all_nodes_pt.end();
                     it++)
                {
                  Vector<double> zeta(1);
                  // Get the node
                  Node* nod_pt = (*it);
                  // Get the boundary coordinate associated to the node
                  nod_pt->get_coordinates_on_boundary(b, zeta);
                  // Compute its new value
                  zeta[0] = segment_arclength[is] - zeta[0];
                  // Set the new boundary coordinate value
                  nod_pt->set_coordinates_on_boundary(b, zeta);
                } // Loop over the nodes in the segment
              } // if (boundary_segment_inverted(b)[is])
              else
              {
                // The boundary is not inverted, do nothing!!!
              }

            } // original zeta values in increasing order
            else if (original_segment_final_zeta <
                     original_segment_initial_zeta)
            {
              // The original values go in decreasing order, only need
              // to change the values if the original segment is NOT
              // marked as inverted
              if (boundary_segment_inverted(b)[is])
              {
                // The boundary is inverted, do nothing!!!
              }
              else
              {
                // The original segment is NOT inverted, then we need
                // to reverse the boundary coordinates. Go through all
                // the nodes and reverse its value
                std::set<Node*> all_nodes_pt = segment_all_nodes_pt[is];
                for (std::set<Node*>::iterator it = all_nodes_pt.begin();
                     it != all_nodes_pt.end();
                     it++)
                {
                  Vector<double> zeta(1);
                  // Get the node
                  Node* nod_pt = (*it);
                  // Get the boundary coordinate associated to the node
                  nod_pt->get_coordinates_on_boundary(b, zeta);
                  // Compute its new value
                  zeta[0] = segment_arclength[is] - zeta[0];
                  // Set the new boundary coordinate value
                  nod_pt->set_coordinates_on_boundary(b, zeta);
                } // Loop over the nodes in the segment
              } // else if (boundary_segment_inverted(b)[is])
            } // original zeta values in decreasing order
#ifdef PARANOID
            else
            {
              std::stringstream error_message;
              std::string output_string = "UnstructuredTwoDMeshGeometryBase::"
                                          "setup_boundary_coordinates()";
              error_message
                << "It was not possible to identify if the zeta values on the\n"
                << "current segment in the boundary should go in increasing\n"
                << "or decreasing order.\n\n"
                << "Boundary: (" << b << ")\n"
                << "Segment: (" << is << ")\n"
                << "Initial zeta value: (" << original_segment_initial_zeta
                << ")\n"
                << "Final zeta value: (" << original_segment_final_zeta
                << ")\n";
              throw OomphLibError(
                error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
            }
#endif

          } // if (this->boundary_geom_object_pt(b)!=0)
          else
          {
            // No GeomObject associated, do nothing!!!

          } // else if (this->boundary_geom_object_pt(b)!=0)
        } // for (is < nsegments)
      } // if (this->is_mesh_distributed() &&
      // Assigned_segments_initial_zeta_values[b])

#endif // #ifdef OOMPH_HAS_MPI

      // Get the number of sets for nodes
#ifdef PARANOID
      if (segment_all_nodes_pt.size() != nsegments)
      {
        std::ostringstream error_message;
        std::string output_string =
          "UnstructuredTwoDMeshGeometryBase::setup_boundary_coordinates()";
        error_message << "The number of segments (" << nsegments
                      << ") and the number of "
                      << "sets of nodes (" << segment_all_nodes_pt.size()
                      << ") representing\n"
                      << "the\nsegments is different!!!\n\n";
        throw OomphLibError(
          error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // --------------------------------------------------------------
      // DISTRIBUTED MESH: END - Check orientation of boundaries and
      // assign boundary coordinates accordingly
      // --------------------------------------------------------------

      // =================================================================
      // BEGIN: Get the lenght of the boundaries or segments of the
      // boundary
      // =================================================================

      // The nodes have been assigned arc-length coordinates from one
      // end or the other of the segments. If the mesh is distributed
      // the values are set so that they agree with the increasing or
      // decreasing order of the zeta values for the segments

      // Storage for the coordinates of the first and last nodes
      Vector<double> first_coordinate(2);
      Vector<double> last_coordinate(2);
      // Storage for the zeta coordinate of the first and last nodes
      Vector<double> first_node_zeta_coordinate(1, 0.0);
      Vector<double> last_node_zeta_coordinate(1, 0.0);

      // Store the accumulated arclength, used at the end to check the
      // max arclength of the boundary
      double boundary_arclength = 0.0;

      // If the mesh is marked as distributed and the initial zeta
      // values for the segments have been computed then get the
      // info. regarding the initial and final nodes coordinates, same
      // as the zeta boundary values for those nodes

#ifdef OOMPH_HAS_MPI
      if (this->is_mesh_distributed() &&
          Assigned_segments_initial_zeta_values[b])
      {
        // --------------------------------------------------------------
        // DISTRIBUTED MESH: BEGIN
        // --------------------------------------------------------------

        // Get the initial and final coordinates for the complete boundary
        first_coordinate = boundary_initial_coordinate(b);
        last_coordinate = boundary_final_coordinate(b);
        // Get the initial and final zeta values for the boundary
        // (arclength)
        first_node_zeta_coordinate = boundary_initial_zeta_coordinate(b);
        last_node_zeta_coordinate = boundary_final_zeta_coordinate(b);

        // The total arclength is the maximum between the initial and
        // final zeta coordinate
        boundary_arclength =
          std::max(first_node_zeta_coordinate[0], last_node_zeta_coordinate[0]);

#ifdef PARANOID
        if (boundary_arclength == 0)
        {
          std::ostringstream error_message;
          std::string output_string =
            "UnstructuredTwoDMeshGeometryBase::setup_boundary_coordinates()";
          error_message << "The boundary arclength is zero for boundary (" << b
                        << ")\n";
          throw OomphLibError(
            error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Check if there is a GeomObject associated to the boundary,
        // and assign the corresponding zeta (arclength) values
        if (this->boundary_geom_object_pt(b) != 0)
        {
          initial_segment_zeta = boundary_segment_initial_zeta(b);
          final_segment_zeta = boundary_segment_final_zeta(b);
        }
        else
        {
          initial_segment_arclength = boundary_segment_initial_arclength(b);
          final_segment_arclength = boundary_segment_final_arclength(b);
        }

        // --------------------------------------------------------------
        // DISTRIBUTED MESH: END
        // --------------------------------------------------------------

      } // if (this->is_mesh_distributed() &&
      //     Assigned_segments_initial_zeta_values[b])
      else
      {
#endif // #ifdef OOMPH_HAS_MPI

        // If the mesh is NOT distributed or the initial and final zeta
        // values of the segments have not been assigned then perform
        // as in the serial case
        if (this->boundary_geom_object_pt(b) != 0)
        {
          initial_segment_zeta[0] = this->boundary_coordinate_limits(b)[0];
          final_segment_zeta[0] = this->boundary_coordinate_limits(b)[1];
        }
        else
        {
          // When the mesh is or not distributed but the initial
          // segment's zeta values HAVE NOT been established then
          // initalize the initial segment to zero
          initial_segment_arclength[0] = 0.0;
        }
#ifdef OOMPH_HAS_MPI
      } // The mesh is NOT distributed or the zeta values for the
      // segments have not been established

#endif // #ifdef OOMPH_HAS_MPI

      // =================================================================
      // END: Get the lenght of the boundaries or segments of the
      // boundary
      // =================================================================

      // =================================================================
      // BEGIN: Scale the boundary coordinates based on whether the
      // boundary has an associated GeomObj or not
      // =================================================================

      // Go through all the segments and assign the scaled boundary
      // coordinates
      for (unsigned is = 0; is < nsegments; is++)
      {
#ifdef PARANOID
        if (segment_sorted_ele_pt[is].size() == 0)
        {
          std::ostringstream error_message;
          std::string output_string =
            "UnstructuredTwoDMeshGeometryBase::setup_boundary_coordinates()";
          error_message << "The (" << is << ")-th segment has no elements\n";
          throw OomphLibError(
            error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
        } // if (segment_sorted_ele_pt[is].size() == 0)x
#endif

        // Get the first and last nodes coordinates for the current
        // segment

#ifdef OOMPH_HAS_MPI
        // Check if the initial zeta values of the segments have been
        // established, if they have not then get the first and last
        // coordinates from the current segments, same as the zeta
        // values
        if (!Assigned_segments_initial_zeta_values[b])
        {
#endif // #ifdef OOMPH_HAS_MPI

          // Get access to the first element on the segment
          FiniteElement* first_ele_pt = segment_sorted_ele_pt[is].front();

          // Number of nodes
          const unsigned nnod = first_ele_pt->nnode();

          // Get the first node of the current segment
          Node* first_node_pt = first_ele_pt->node_pt(0);
          if (is_inverted[first_ele_pt])
          {
            first_node_pt = first_ele_pt->node_pt(nnod - 1);
          }

          // Get access to the last element on the segment
          FiniteElement* last_ele_pt = segment_sorted_ele_pt[is].back();

          // Get the last node of the current segment
          Node* last_node_pt = last_ele_pt->node_pt(nnod - 1);
          if (is_inverted[last_ele_pt])
          {
            last_node_pt = last_ele_pt->node_pt(0);
          }

          // Get the coordinates for the first and last node
          for (unsigned i = 0; i < 2; i++)
          {
            first_coordinate[i] = first_node_pt->x(i);
            last_coordinate[i] = last_node_pt->x(i);
          }

          // Get the zeta coordinates for the first and last node
          first_node_pt->get_coordinates_on_boundary(
            b, first_node_zeta_coordinate);
          last_node_pt->get_coordinates_on_boundary(b,
                                                    last_node_zeta_coordinate);

#ifdef OOMPH_HAS_MPI
        } // if (!Assigned_segments_initial_zeta_values[b])
#endif // #ifdef OOMPH_HAS_MPI

        // Get the nodes of the current segment
        std::set<Node*> all_nodes_pt = segment_all_nodes_pt[is];

        // If the boundary has a geometric object representation then
        // scale the coordinates to match those of the geometric object
        GeomObject* const geom_object_pt = this->boundary_geom_object_pt(b);
        if (geom_object_pt != 0)
        {
          Vector<double> bound_coord_limits =
            this->boundary_coordinate_limits(b);
          // Get the position of the ends of the geometric object
          Vector<double> zeta(1);
          Vector<double> first_geom_object_location(2);
          Vector<double> last_geom_object_location(2);

          // Get the zeta value for the initial coordinates of the
          // GeomObject
          zeta[0] = bound_coord_limits[0];
          // zeta[0] = initial_segment_zeta[is];

          // Get the coordinates from the initial zeta value from the
          // GeomObject
          geom_object_pt->position(zeta, first_geom_object_location);

          // Get the zeta value for the final coordinates of the
          // GeomObject
          zeta[0] = bound_coord_limits[1];
          // zeta[0] = final_segment_zeta[is];

          // Get the coordinates from the final zeta value from the
          // GeomObject
          geom_object_pt->position(zeta, last_geom_object_location);

          // Calculate the errors in position between the first and
          // last nodes and the endpoints of the geometric object
          double error = 0.0;
          double tmp_error = 0.0;
          for (unsigned i = 0; i < 2; i++)
          {
            const double dist =
              first_geom_object_location[i] - first_coordinate[i];
            tmp_error += dist * dist;
          }
          error += sqrt(tmp_error);
          tmp_error = 0.0;
          for (unsigned i = 0; i < 2; i++)
          {
            const double dist =
              last_geom_object_location[i] - last_coordinate[i];
            tmp_error += dist * dist;
          }
          error += sqrt(tmp_error);

          // Calculate the errors in position between the first and
          // last nodes and the endpoints of the geometric object if
          // reversed
          double rev_error = 0.0;
          tmp_error = 0.0;
          for (unsigned i = 0; i < 2; i++)
          {
            const double dist =
              first_geom_object_location[i] - last_coordinate[i];
            tmp_error += dist * dist;
          }
          rev_error += sqrt(tmp_error);
          tmp_error = 0.0;
          for (unsigned i = 0; i < 2; i++)
          {
            const double dist =
              last_geom_object_location[i] - first_coordinate[i];
            tmp_error += dist * dist;
          }
          rev_error += sqrt(tmp_error);

          // Number of polyline vertices along this boundary
          const unsigned n_vertex = Polygonal_vertex_arclength_info[b].size();

          // Get polygonal vertex data
          Vector<std::pair<double, double>> polygonal_vertex_arclength =
            Polygonal_vertex_arclength_info[b];

          // If the (normal) error is small than reversed then we have
          // the coordinate direction correct. If not then we must
          // reverse it bool reversed = false;
          if (error < rev_error)
          {
            // Coordinates are aligned (btw: don't delete this block --
            // there's a final else below to catch errors!)

            // reversed = false;
          }
          else if (error > rev_error)
          {
            // reversed = true;

            // Reverse the limits of the boundary coordinates along the
            // geometric object
            double temp = bound_coord_limits[0];
            bound_coord_limits[0] = bound_coord_limits[1];
            bound_coord_limits[1] = temp;

#ifdef OOMPH_HAS_MPI
            // If we are working with a NON distributed mesh then
            // re-assign the initial and final zeta values
            if (!this->is_mesh_distributed())
            {
#endif
              temp = initial_segment_zeta[is];
              initial_segment_zeta[is] = final_segment_zeta[is];
              final_segment_zeta[is] = temp;
#ifdef OOMPH_HAS_MPI
            }
#endif
            // Reverse the vertices information
            for (unsigned v = 0; v < n_vertex; v++)
            {
              polygonal_vertex_arclength[v].first =
                Polygonal_vertex_arclength_info[b][v].first;

              polygonal_vertex_arclength[v].second =
                Polygonal_vertex_arclength_info[b][n_vertex - v - 1].second;
            } // for (v < n_vertex)
          }
          else
          {
            std::ostringstream error_stream;
            std::string output_string =
              "UnstructuredTwoDMeshGeometryBase::setup_boundary_coordinates()";
            error_stream
              << "Something very strange has happened.\n"
              << "The error between the endpoints of the geometric object\n"
              << "and the first and last nodes on the boundary is the same\n"
              << "irrespective of the direction of the coordinate.\n"
              << "This probably means that things are way off.\n"
              << "The errors are " << error << " and " << rev_error << "\n";
            std::cout << error_stream.str();
            throw OomphLibError(
              error_stream.str(), output_string, OOMPH_EXCEPTION_LOCATION);
          }

          // Get the total arclength of the edge
          // last_node_pt->get_coordinates_on_boundary(b, zeta);
          // double zeta_old_range=zeta[0];

          // Store the arclength of the segment (not yet mapped to the
          // boundary coordinates defined by the GeomObject)
          const double zeta_old_range = segment_arclength[is];

          // double zeta_new_range=bound_coord_limits[1]-bound_coord_limits[0];

          // The range to map the zeta values
          double zeta_new_range =
            final_segment_zeta[is] - initial_segment_zeta[is];
          // The initial zeta value for the current segment
          double initial_local_segment_zeta = initial_segment_zeta[is];

#ifdef OOMPH_HAS_MPI

          // If we are working with a distributed mes and the initial
          // and final zeta values for the segments have been
          // established then reset the range to map
          if (this->is_mesh_distributed() &&
              Assigned_segments_initial_zeta_values[b])
          {
            // Re-set the range to map the zeta values
            zeta_new_range =
              std::fabs(final_segment_zeta[is] - initial_segment_zeta[is]);

            // Re-set the initial zeta value of the segment
            initial_local_segment_zeta =
              std::min(initial_segment_zeta[is], final_segment_zeta[is]);
          }

#endif

          // Re-assign boundary coordinate for the case where boundary
          // is represented by polygon
          unsigned use_old = false;
          if (n_vertex == 0) use_old = true;

          // Now scale the coordinates accordingly
          for (std::set<Node*>::iterator it = all_nodes_pt.begin();
               it != all_nodes_pt.end();
               it++)
          {
            // Get the node
            Node* nod_pt = (*it);

            // Get coordinate based on arclength along polygonal repesentation
            nod_pt->get_coordinates_on_boundary(b, zeta);

            if (use_old)
            {
              // Boundary is actually a polygon -- simply rescale
              zeta[0] = initial_local_segment_zeta +
                        (zeta_new_range / zeta_old_range) * (zeta[0]);
            }
            else
            {
              // Scale such that vertex nodes stay where they were
              bool found = false;

              // Loop over vertex nodes
              for (unsigned v = 1; v < n_vertex; v++)
              {
                if ((zeta[0] >= polygonal_vertex_arclength[v - 1].first) &&
                    (zeta[0] <= polygonal_vertex_arclength[v].first))
                {
                  // Increment in intrinsic coordinate along geom object
                  double delta_zeta =
                    (polygonal_vertex_arclength[v].second -
                     polygonal_vertex_arclength[v - 1].second);
                  // Increment in arclength along segment
                  double delta_polyarc =
                    (polygonal_vertex_arclength[v].first -
                     polygonal_vertex_arclength[v - 1].first);

                  // Mapped arclength coordinate
                  double zeta_new =
                    polygonal_vertex_arclength[v - 1].second +
                    delta_zeta *
                      (zeta[0] - polygonal_vertex_arclength[v - 1].first) /
                      delta_polyarc;
                  zeta[0] = zeta_new;

                  // Success!
                  found = true;

                  // Bail out
                  break;
                }
              } // for (v < n_vertex)

              // If we still haven't found it's probably the last point along
              if (!found)
              {
#ifdef PARANOID
                double diff = std::fabs(
                  zeta[0] - polygonal_vertex_arclength[n_vertex - 1].first);
                if (diff >
                    ToleranceForVertexMismatchInPolygons::Tolerable_error)
                {
                  std::ostringstream error_stream;
                  error_stream
                    << "Wasn't able to locate the polygonal arclength exactly\n"
                    << "during re-setup of boundary coordinates and have\n"
                    << "assumed that we're dealing with the final point along\n"
                    << "the curvilinear segment and encountered some roundoff\n"
                    << "However,the difference in the polygonal zeta "
                       "coordinates\n"
                    << "between zeta[0] " << zeta[0]
                    << " and the originallly stored value "
                    << polygonal_vertex_arclength[n_vertex - 1].first << "\n"
                    << "is " << diff
                    << " which exceeds the threshold specified\n"
                    << "in the publically modifiable variable\n"
                    << "ToleranceForVertexMismatchInPolygons::Tolerable_error\n"
                    << "whose current value is: "
                    << ToleranceForVertexMismatchInPolygons::Tolerable_error
                    << "\nPlease check your mesh carefully and increase the\n"
                    << "threshold if you're sure this is appropriate\n";
                  throw OomphLibError(error_stream.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
#endif
                zeta[0] = polygonal_vertex_arclength[n_vertex - 1].second;
              }
            }

            // Assign updated coordinate
            nod_pt->set_coordinates_on_boundary(b, zeta);
          }
        } // if (geom_object_pt != 0)
        else
        {
          // Establish the initial zeta value for this segment
          double z_initial = initial_segment_arclength[is];

          // Only use end points of the whole boundary and pick the
          // bottom left node

          // Set the bottom left coordinate as the first coordinates
          Vector<double> bottom_left_coordinate(2);
          bottom_left_coordinate = first_coordinate;

          // ... do the same with the zeta value
          Vector<double> bottom_left_zeta_coordinate(1);
          bottom_left_zeta_coordinate = first_node_zeta_coordinate;

          // Is the last y-coordinate smaller than y-coordinate of the
          // current bottom left coordinate
          if (last_coordinate[1] < bottom_left_coordinate[1])
          {
            // Re-set the bottom-left coordinate as the last coordinate
            bottom_left_coordinate = last_coordinate;

            // ... do the same with the zeta value
            bottom_left_zeta_coordinate = last_node_zeta_coordinate;
          }
          // The y-coordinates are the same
          else if (last_coordinate[1] == bottom_left_coordinate[1])
          {
            // Then check for the x-coordinate, which is the most-left
            if (last_coordinate[0] < bottom_left_coordinate[0])
            {
              // Re-set the bottom-left coordinate as the last
              // coordinate
              bottom_left_coordinate = last_coordinate;

              // ... do the same with the zeta value
              bottom_left_zeta_coordinate = last_node_zeta_coordinate;
            }
          } // else (The y-coordinates are the same)

          Vector<double> zeta(1, 0.0);

          // Now adjust boundary coordinate so that the bottom left
          // node has a boundary coordinate of 0.0 and that zeta
          // increases away from that point
          zeta = bottom_left_zeta_coordinate;
          const double zeta_ref = zeta[0];

          // Also get the maximum zeta value
          double zeta_max = 0.0;
          for (std::set<Node*>::iterator it = all_nodes_pt.begin();
               it != all_nodes_pt.end();
               it++)
          {
            Node* nod_pt = (*it);
            nod_pt->get_coordinates_on_boundary(b, zeta);

#ifdef OOMPH_HAS_MPI

            // If the mesh is distributed and the initial and final
            // zeta values for the segment have been assigned then
            // check if the segment is inverted, we need to consider
            // that to correctly assgin the zeta values for the segment
            if (this->is_mesh_distributed() &&
                Assigned_segments_initial_zeta_values[b])
            {
              // Is the segment inverted, if that is the case then
              // invert the zeta values
              if (boundary_segment_inverted(b)[is])
              {
                zeta[0] = segment_arclength[is] - zeta[0];
              } // if (boundary_segment_inverted(b)[is])
            }

#endif // #ifdef OOMPH_HAS_MPI

            // Set the zeta value
            zeta[0] += z_initial;
            // Adjust the value based on the bottom-left condition
            zeta[0] -= zeta_ref;
            // If direction is reversed, then take absolute value
            if (zeta[0] < 0.0)
            {
              zeta[0] = std::fabs(zeta[0]);
            }
            // Get the max zeta value (we will use it to scale the
            // values to [0,1])
            if (zeta[0] > zeta_max)
            {
              zeta_max = zeta[0];
            }
            nod_pt->set_coordinates_on_boundary(b, zeta);
          } // Loop through the nodes in the segment (boundary)

#ifdef OOMPH_HAS_MPI
          // After assigning boundary coordinates, BUT before scaling,
          // copy the initial and final segment arclengths so that we
          // know if the values go in increasing or decreasing
          // order. This will be used to identify the correct direction
          // of the segments in the new meshes created in the
          // adaptation method.
          if (this->is_mesh_distributed() &&
              Assigned_segments_initial_zeta_values[b])
          {
            // Get the first face element
            FiniteElement* first_seg_ele_pt = segment_sorted_ele_pt[is].front();

#ifdef PARANOID
            // Check if the face element is nonhalo, it shouldn't, but
            // better check
            if (first_seg_ele_pt->is_halo())
            {
              std::ostringstream error_message;
              std::string output_string = "UnstructuredTwoDMeshGeometryBase::"
                                          "setup_boundary_coordinates()";
              error_message << "The first face element in the (" << is
                            << ")-th segment"
                            << " is halo\n";
              throw OomphLibError(
                error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
            } // if (first_seg_ele_pt->is_halo())
#endif

            // Number of nodes
            const unsigned nnod = first_seg_ele_pt->nnode();

            // Get the first node of the current segment
            Node* first_seg_node_pt = first_seg_ele_pt->node_pt(0);
            if (is_inverted[first_seg_ele_pt])
            {
              first_seg_node_pt = first_seg_ele_pt->node_pt(nnod - 1);
            }

            // Get the last face element
            FiniteElement* last_seg_ele_pt = segment_sorted_ele_pt[is].back();

#ifdef PARANOID
            // Check if the face element is nonhalo, it shouldn't, but
            // better check
            if (last_seg_ele_pt->is_halo())
            {
              std::ostringstream error_message;
              std::string output_string = "UnstructuredTwoDMeshGeometryBase::"
                                          "setup_boundary_coordinates()";
              error_message << "The last face element in the (" << is
                            << ")-th segment"
                            << " is halo\n";
              throw OomphLibError(
                error_message.str(), output_string, OOMPH_EXCEPTION_LOCATION);
            } // if (last_seg_ele_pt->is_halo())
#endif

            // Get the last node of the current segment
            Node* last_seg_node_pt = last_seg_ele_pt->node_pt(nnod - 1);
            if (is_inverted[last_seg_ele_pt])
            {
              last_seg_node_pt = last_seg_ele_pt->node_pt(0);
            }

            // Now get the first and last node boundary coordinate values
            Vector<double> first_seg_arclen(1);
            Vector<double> last_seg_arclen(1);

            first_seg_node_pt->get_coordinates_on_boundary(b, first_seg_arclen);
            last_seg_node_pt->get_coordinates_on_boundary(b, last_seg_arclen);

            // Update the initial and final segments arclength
            boundary_segment_initial_arclength(b)[is] = first_seg_arclen[0];
            boundary_segment_final_arclength(b)[is] = last_seg_arclen[0];

            // Update the initial and final coordinates
            Vector<double> updated_segment_initial_coord(2);
            Vector<double> updated_segment_final_coord(2);
            for (unsigned k = 0; k < 2; k++)
            {
              updated_segment_initial_coord[k] = first_seg_node_pt->x(k);
              updated_segment_final_coord[k] = last_seg_node_pt->x(k);
            }

            // Set the updated initial coordinate
            boundary_segment_initial_coordinate(b)[is] =
              updated_segment_initial_coord;

            // Set the updated final coordinate
            boundary_segment_final_coordinate(b)[is] =
              updated_segment_final_coord;

          } // if (this->is_mesh_distributed() &&
          // Assigned_segments_initial_zeta_values[b])
#endif // #ifdef OOMPH_HAS_MPI

          // The max. value will be incorrect if we are working with
          // distributed meshes where the current boundary has been
          // split by the distribution process. Here correct the
          // maximum value
          if (zeta_max < boundary_arclength)
          {
            zeta_max = boundary_arclength;
          }

          // Scale all surface coordinates so that all the points be on
          // the range [0, 1]
          for (std::set<Node*>::iterator it = all_nodes_pt.begin();
               it != all_nodes_pt.end();
               it++)
          {
            // Get the node
            Node* nod_pt = (*it);

            // Get the boundary coordinate
            nod_pt->get_coordinates_on_boundary(b, zeta);

            // Scale the value of the current node
            zeta[0] /= zeta_max;

            // Set the new scaled value
            nod_pt->set_coordinates_on_boundary(b, zeta);
          } // Loop over the nodes

        } // else if (geom_object_pt != 0)

      } // for (is < nsegments)

      // =================================================================
      // END: Scale the boundary coordinates based on whether the
      // boundary has an associated GeomObj or not
      // =================================================================

      // Cleanup
      for (unsigned e = 0; e < nel; e++)
      {
        delete face_el_pt[e];
        face_el_pt[e] = 0;
      }

    } // if (nel > 0), from the beginning of the method

    // Indicate that boundary coordinate has been set up
    Boundary_coordinate_exists[b] = true;
  }


  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////


  //=================================================
  ///  Helper namespace for BCInfo object used
  /// in the identification of boundary elements.
  //=================================================
  namespace TriangleBoundaryHelper
  {
    /// Structure for Boundary Informations
    struct BCInfo
    {
      /// Face ID
      unsigned Face_id;

      /// Boundary ID
      unsigned Boundary;

      /// Pointer to bulk finite element
      FiniteElement* FE_pt;
    };

  } // namespace TriangleBoundaryHelper

} // namespace oomph

#endif
