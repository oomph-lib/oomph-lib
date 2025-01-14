#ifndef DOMAIN_BOUNDARIES_HEADER
#define DOMAIN_BOUNDARIES_HEADER

#include "generic.h"
#include "meshes/triangle_mesh.h"

namespace oomph
{
  template<class ELEMENT>
  RefineableTriangleMesh<ELEMENT>* create_sector_mesh(double contact_angle,
                                            double radius,
                                            unsigned nsegment,
                                            double element_area)
  {
    enum
    {
      Slip_boundary_id,
      Far_field_boundary_id,
      Free_surface_boundary_id,
    };

    const double x_center = 1.0;
    const double y_center = 0;

    Circle* circle_pt = new Circle(x_center, y_center, radius);

    TriangleMeshCurveSection* polyline_pt = 0;

    Vector<Vector<double>> vertex_coord;
    Vector<double> current_vertex(2);
    Vector<double> starting_vertex(2);

    starting_vertex[0] = x_center;
    starting_vertex[1] = y_center;
    vertex_coord.push_back(starting_vertex);
    current_vertex[0] = x_center;
    current_vertex[1] = y_center + radius;
    vertex_coord.push_back(current_vertex);
    polyline_pt = new TriangleMeshPolyLine(vertex_coord, Slip_boundary_id);

    Vector<TriangleMeshCurveSection*> boundary_polyline_pt;
    boundary_polyline_pt.push_back(polyline_pt);
    vertex_coord.clear();

    const double zeta_start = 0.5 * MathematicalConstants::Pi;
    const double zeta_end = 0.5 * MathematicalConstants::Pi + contact_angle;

    polyline_pt = new TriangleMeshCurviLine(
      circle_pt, zeta_start, zeta_end, nsegment, Far_field_boundary_id);

    boundary_polyline_pt.push_back(polyline_pt);

    current_vertex[0] = x_center - radius * sin(contact_angle);
    current_vertex[1] = y_center + radius * cos(contact_angle);
    vertex_coord.push_back(current_vertex);
    vertex_coord.push_back(starting_vertex);
    polyline_pt =
      new TriangleMeshPolyLine(vertex_coord, Free_surface_boundary_id);

    boundary_polyline_pt.push_back(polyline_pt);

    // Create an interior point: Assuming the contact angle is alpha > 1e-2
    // and the radius is > 0.5
    Vector<double> internal_point_pt(2);
    internal_point_pt[0] = 0.5;
    internal_point_pt[1] = 0;

    // Set up a bool for whether the internal point is fixed
    const bool is_internal_point_fixed = true;

    // True the curve sections into a closed curve
    TriangleMeshClosedCurve* outer_closed_curve_pt =
      new TriangleMeshClosedCurve(
        boundary_polyline_pt, internal_point_pt, is_internal_point_fixed);

    // WARNING: outer_closed_curve_pt's destruction is handled by
    // TriangleMeshParameters. This is unexpected behaviour.
    TriangleMeshParameters triangle_mesh_parameters(outer_closed_curve_pt);
    triangle_mesh_parameters.element_area() = element_area;

    // For some reason the vertices are 1.7e-13 apart!?
    ToleranceForVertexMismatchInPolygons::Tolerable_error = 2e-13;

    return new RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters);
  }
} // namespace oomph
#endif
