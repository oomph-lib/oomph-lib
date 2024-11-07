#ifndef DOMAIN_BOUNDARIES_HEADER
#define DOMAIN_BOUNDARIES_HEADER

#include "generic.h"
#include "meshes/triangle_mesh.h"
#include "parameter_values.h"

namespace oomph
{
  enum
  {
    Slip_boundary_id,
    Far_field_boundary_id,
    Free_surface_boundary_id,
    Inner_boundary_id,
  };


  void create_sector_domain(
    double contact_angle,
    Vector<TriangleMeshCurveSection*>& boundary_polyline_pt)
  {
    const double radius = 1.0;
    const unsigned nsegment = parameters::nsegment;
    const double x_center = 2;
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
  }

  void create_sector_domain_with_interface(
    double contact_angle,
    TriangleMeshClosedCurve*& outer_closed_curve_pt,
    Vector<TriangleMeshOpenCurve*>& internal_open_curves_pt)
  {
    const double radius = 1.0;
    const double inner_radius = parameters::inner_radius;
    const unsigned nsegment = parameters::nsegment;
    const double x_center = 2;
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
    current_vertex[1] = y_center + inner_radius;
    vertex_coord.push_back(current_vertex);

    current_vertex[0] = x_center;
    current_vertex[1] = y_center + radius;
    vertex_coord.push_back(current_vertex);

    polyline_pt = new TriangleMeshPolyLine(vertex_coord, Slip_boundary_id);

    // Setup the sector boundaries
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

    current_vertex[0] = x_center - inner_radius * sin(contact_angle);
    current_vertex[1] = y_center + inner_radius * cos(contact_angle);
    vertex_coord.push_back(current_vertex);

    vertex_coord.push_back(starting_vertex);
    polyline_pt =
      new TriangleMeshPolyLine(vertex_coord, Free_surface_boundary_id);

    boundary_polyline_pt.push_back(polyline_pt);

    // Create an interior point: Assuming the contact angle is alpha > 1e-2
    // and the radius is > 0.5
    Vector<double> internal_point_pt(2);
    internal_point_pt[0] = 2 - 0.5 * sin(1e-2 * MathematicalConstants::Pi);
    internal_point_pt[1] = 0.5 * cos(1e-2 * MathematicalConstants::Pi);

    // Set up a bool for whether the internal point is fixed
    const bool is_internal_point_fixed = true;

    // True the curve sections into a closed curve
    outer_closed_curve_pt = new TriangleMeshClosedCurve(
      boundary_polyline_pt, internal_point_pt, is_internal_point_fixed);

    Circle* inner_circle_pt = new Circle(x_center, y_center, inner_radius);
    Vector<TriangleMeshCurveSection*> interface_curve_pt;
    interface_curve_pt.push_back(new TriangleMeshCurviLine(
      inner_circle_pt, zeta_start, zeta_end, nsegment, Inner_boundary_id));

    interface_curve_pt[0]->connect_initial_vertex_to_polyline(
      dynamic_cast<TriangleMeshPolyLine*>(boundary_polyline_pt[0]), 1);
    interface_curve_pt[0]->connect_final_vertex_to_polyline(
      dynamic_cast<TriangleMeshPolyLine*>(boundary_polyline_pt[2]), 1);

    internal_open_curves_pt.push_back(
      new TriangleMeshOpenCurve(interface_curve_pt));
  }
} // namespace oomph
#endif
