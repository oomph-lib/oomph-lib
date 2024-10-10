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

#ifndef BO_HEIGHT_CONTROL_SINGULAR_AXISYM_DYNAMIC_CAP_PROBLEM_HEADER
#define BO_HEIGHT_CONTROL_SINGULAR_AXISYM_DYNAMIC_CAP_PROBLEM_HEADER

#include <chrono>
#include <algorithm>
#include <numeric>

// OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"

// The mesh
#include "meshes/triangle_mesh.h"

#include "my_error_estimator.h"

#include "axisym_fluid_slip_elements.h"

#include "axisym_navier_stokes_flux_elements.h"

#include "net_flux_elements.h"

#include "height_element.h"

#include "parameters.h"


#include "complex_less.h"

#include "volume_constraint_elements_with_output.h"

#include "my_eigenproblem.h"

#include "free_surface_elements.h"

#include "utility_functions.h"

#include "hijacked_projectable_axisymmteric_Ttaylor_hood_elements.h"

#include "eigensolution_functions.h"
#include "eigensolution_elements.h"
#include "singular_fluid_traction_elements.h"
#include "parameter_values.h"
#include "parameter_functions.h"


namespace oomph
{
  //===========start_of_axisym_dynamic_cap_problem_class========================
  // A class that solves the Axisymmetric Navier--Stokes equations
  // to compute the shape of a dynamic interface in a cylindrical container
  // with imposed contact angle at the boundary.
  //============================================================================
  template<class ELEMENT, class TIMESTEPPER>
  class BoHeightControlSingularAxisymDynamicCapProblem : public MyEigenproblem
  {
  public:
    // Typedef all the different types of elements used for clarity.
    // I am using capitals, even though these are not templates.
    typedef FreeSurfaceElement<ELEMENT> FREE_SURFACE_ELEMENT;
    typedef DebugElasticAxisymmetricVolumeConstraintBoundingElement<ELEMENT>
      VOLUME_COMPUTATION_ELEMENT;
    typedef VolumeConstraintElementWithOutput VOLUME_CONSTRAINT_ELEMENT;
    typedef AxisymmetricNavierStokesSlipElement<ELEMENT> SLIP_ELEMENT;
    typedef DebugImposeImpenetrabilityElement<ELEMENT> NO_PENETRATION_ELEMENT;
    typedef AxisymmetricNavierStokesFluxElement<ELEMENT> FLUX_ELEMENT;
    typedef NetFluxElement NET_FLUX_ELEMENT;
    typedef ElasticPointFluidInterfaceBoundingElement<ELEMENT>
      CONTACT_LINE_ELEMENT;
    typedef HeightElement HEIGHT_ELEMENT;

    // Constructor
    BoHeightControlSingularAxisymDynamicCapProblem(const double& contact_angle,
                                                   const bool& has_restart)
      : Ca(1.0),
        Bo(1.0),
        Re(0.0),
        ReSt(0.0),
        Viscosity_ratio(1.0),
        St(1.0),
        ReInvFr_data_pt(new Data(1)),
        Volume((Global_Physical_Parameters::Initial_fluid_height) /
               2.0), // Initialise the value of the volume
                     // the physical volume divided by 2pi
        Pext(1.0), // Initialise the external pressure to some random value
        Contact_angle(contact_angle), // Initialise the contact angle
        Right_angle(0.5 * MathematicalConstants ::Pi),
        Zero_sigma(0),
        Max_adapt(1), // Initialise the maximum number of adapt steps
        Wall_velocity_data_pt(new Data(1)),
        Height_step_data_pt(ReInvFr_data_pt),
        Contact_line_solid_node_pt(0),
        Inner_corner_solid_node_pt(0),
        Z2_error_estimator_pt(0),
        Corner_error_estimator_pt(0),
        Free_surface_mesh_pt(0),
        Slip_boundary_mesh_pt(0),
        No_penetration_boundary_mesh_pt(0),
        Flux_mesh_pt(0),
        Volume_computation_mesh_pt(0),
        Contact_angle_mesh_pt(0),
        Eigensolution_slip_mesh_pt(0),
        Singularity_scaling_mesh_pt(0),
        Pressure_contribution_mesh_1_pt(0),
        Pressure_contribution_mesh_2_pt(0),
        Pressure_contribution_geom_mesh_1_pt(0),
        Pressure_contribution_geom_mesh_2_pt(0),
        Volume_constraint_mesh_pt(0),
        Net_flux_mesh_pt(0),
        Height_mesh_pt(0),
        Backup_volume_constraint_lagrange_multiplier(1.0),
        Backup_net_flux_lagrange_multiplier(1.0),
        Backup_point_kinematic_lagrange_multiplier(1.0),
        problem_type(Normal_problem),
        Is_steady(true),
        Using_contact_angle_error_estimator(false)
    {
      //======================================================================
      // Set up the rest of the parameters
      //======================================================================

      this->Start_time = std::chrono::high_resolution_clock::now();

      // Create time stepper
      const bool is_timestepper_adaptive = true;
      this->add_time_stepper_pt(new TIMESTEPPER(is_timestepper_adaptive));
      this->Maximum_dt = Mesh_Control_Parameters::Max_timestep;
      // Set the constituive law
      this->Constitutive_law_pt =
        new GeneralisedHookean(&Mesh_Control_Parameters::Nu);
      // Set the maximum number of newton iterations
      this->max_newton_iterations() =
        Newton_Solver_Parameters::Max_newton_iterations;
      // Set the maximum residual before assuming newton is not converging
      this->max_residuals() = Newton_Solver_Parameters::Max_residual;

      this->newton_solver_tolerance() =
        Newton_Solver_Parameters::Newton_solver_tolerance;
      // Suppress warning for restarting
      this->Suppress_warning_about_actions_before_read_unstructured_meshes =
        true;

      this->Desired_newton_iterations_ds = 3;
      this->Desired_proportion_of_arc_length = 0.01;

      // Create DocInfo object (allows checking if output directory exists)
      this->doc_info().number() = 0;

      //======================================================================
      // Create the refineable bulk mesh
      //======================================================================
      // Create the bulk mesh and its elements
      if (has_restart)
      {
        create_simple_bulk_mesh();
      }
      else
      {
        create_bulk_mesh();
      }

      // Output boundary and mesh initial mesh for information
      this->Bulk_mesh_pt->output_boundaries(this->doc_info().directory() +
                                            "/boundaries.dat");
      this->Bulk_mesh_pt->output(this->doc_info().directory() + "/mesh.dat");

      //======================================================================
      // Create the error estimators
      //======================================================================
      set_contact_line_node_pt();
      create_error_estimators();

      //======================================================================
      // Set up external pressure
      //======================================================================
      // Create a Data object whose single value stores the
      // external pressure
      External_pressure_data_pt = new Data(1);

      // Set external pressure
      External_pressure_data_pt->set_value(0, Pext);

      // Regard the external pressure is an unknown and add
      // it to the problem's global data so it gets included
      // in the equation numbering. Note that, at the moment,
      // there's no equation that determines its value!
      add_global_data(External_pressure_data_pt);

      // Setup Eotvos number
      ReInvFr_data_pt->pin(0);

      // Setup wall velocity
      Wall_velocity_data_pt->set_value(0, Slip_Parameters::wall_velocity);
      Slip_Parameters::wall_velocity_pt = Wall_velocity_data_pt->value_pt(0);
      Wall_velocity_data_pt->pin(0);

      //======================================================================
      // Create the empty non-refineable meshes and add to problem
      //======================================================================
      // Allocate storage the sub-meshes
      if (problem_type == Normal_problem)
      {
        Free_surface_mesh_pt = new Mesh;
        Slip_boundary_mesh_pt = new Mesh;
        Flux_mesh_pt = new Mesh;
        Volume_computation_mesh_pt = new Mesh;
        Contact_angle_mesh_pt = new Mesh;
        Volume_constraint_mesh_pt = new Mesh;
        Net_flux_mesh_pt = new Mesh;
        No_penetration_boundary_mesh_pt = new Mesh;
        Height_mesh_pt = new Mesh;
      }

      if (Net_flux_mesh_pt)
      {
        this->add_sub_mesh(Net_flux_mesh_pt);
      }

      this->add_sub_mesh(Bulk_mesh_pt);


      // Now need to add all the meshes
      if (Free_surface_mesh_pt)
      {
        this->add_sub_mesh(Free_surface_mesh_pt);
      }
      if (Slip_boundary_mesh_pt)
      {
        this->add_sub_mesh(Slip_boundary_mesh_pt);
      }
      if (No_penetration_boundary_mesh_pt)
      {
        this->add_sub_mesh(No_penetration_boundary_mesh_pt);
      }
      if (Flux_mesh_pt)
      {
        this->add_sub_mesh(Flux_mesh_pt);
      }
      if (Volume_computation_mesh_pt)
      {
        this->add_sub_mesh(Volume_computation_mesh_pt);
      }
      if (Contact_angle_mesh_pt)
      {
        this->add_sub_mesh(Contact_angle_mesh_pt);
      }
      if (Volume_constraint_mesh_pt)
      {
        this->add_sub_mesh(Volume_constraint_mesh_pt);
      }
      if (Height_mesh_pt)
      {
        this->add_sub_mesh(Height_mesh_pt);
      }

      Eigensolution_slip_mesh_pt = new Mesh;
      this->add_sub_mesh(Eigensolution_slip_mesh_pt);
      Singularity_scaling_mesh_pt = new Mesh;
      this->add_sub_mesh(Singularity_scaling_mesh_pt);
      Pressure_contribution_mesh_1_pt = new Mesh;
      this->add_sub_mesh(Pressure_contribution_mesh_1_pt);
      Pressure_contribution_mesh_2_pt = new Mesh;
      this->add_sub_mesh(Pressure_contribution_mesh_2_pt);

      // Add the global data
      this->add_global_data(ReInvFr_data_pt);
      this->add_global_data(Wall_velocity_data_pt);

      // and build the global mesh.
      // We can either build the global mesh here or create an empty mesh and
      // "rebuild" the global mesh in actions after adapt
      this->mesh_pt() = new Mesh;

      //======================================================================
      // Call actions after adapt to create the non-refineable elements and
      // setup the remainder of the problem
      //======================================================================
      actions_after_adapt();

    } // end_of_constructor


  private:
    // Create the bulk mesh and its elements
    void create_bulk_mesh()
    {
      oomph_info << "Creating the initial mesh" << std::endl;

      // Create the Outer_boundary_polyline_pt
      // If we have a 90 degree contact angle
      if (Contact_angle == 0.5 * MathematicalConstants::Pi)
      {
        // create a rectangular domain.
        create_rectangle_domain();
      }
      else
      {
        // Use a circular meniscus as a first approximation.
        create_circular_domain();
      }

      // Set the Outer_boundary_polyline_pt tolerances

      // Set a measure of the maximum local curvature before refining. Default
      // = 0.08. If `d` is the argument, the radius of curvature is = 0.5/d *
      // sqrt(d^2 + 0.5^2)
      Outer_boundary_polyline_pt->set_polyline_refinement_tolerance(
        Mesh_Control_Parameters::Polyline_refinement_tolerence);

      // Set a measure of the minimum local curvature before unrefining.
      // Default = 0.04. If `d` is the argument, the radius of curvature is =
      // 0.5/d * sqrt(d^2 + 0.5^2)
      Outer_boundary_polyline_pt->set_polyline_unrefinement_tolerance(
        Mesh_Control_Parameters::Polyline_unrefinement_tolerence);

      // Now build the mesh, based on the boundaries specified by
      //---------------------------------------------------------
      // polygons just created
      //----------------------

      // Convert to "closed curve" objects
      TriangleMeshClosedCurve* outer_closed_curve_pt =
        Outer_boundary_polyline_pt;

      // Use the TriangleMeshParameter object for gathering all
      // the necessary arguments for the TriangleMesh object
      TriangleMeshParameters triangle_mesh_parameters(outer_closed_curve_pt);

      // Define the maximum element area
      triangle_mesh_parameters.element_area() =
        Mesh_Control_Parameters::Max_element_size;

      // Construct mesh
      Bulk_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(
        triangle_mesh_parameters, this->time_stepper_pt());
      // Bulk_mesh_pt->set_print_level_timings_adaptation(3);

      Bulk_mesh_pt->max_element_size() =
        Mesh_Control_Parameters::Max_element_size;
      Bulk_mesh_pt->min_element_size() =
        Mesh_Control_Parameters::Min_element_size;
      Bulk_mesh_pt->min_permitted_angle() =
        Mesh_Control_Parameters::Min_permitted_angle;
      Bulk_mesh_pt->max_keep_unrefined() = 400;

      refine_mesh_for_weak_contact_angle_constraint();
    }

    void create_simple_bulk_mesh()
    {
      oomph_info << "Creating the simple initial mesh" << std::endl;

      create_rectangle_domain();

      Outer_boundary_polyline_pt->set_polyline_refinement_tolerance(
        Mesh_Control_Parameters::Polyline_refinement_tolerence);
      Outer_boundary_polyline_pt->set_polyline_unrefinement_tolerance(
        Mesh_Control_Parameters::Polyline_unrefinement_tolerence);

      TriangleMeshClosedCurve* outer_closed_curve_pt =
        Outer_boundary_polyline_pt;

      TriangleMeshParameters triangle_mesh_parameters(outer_closed_curve_pt);

      triangle_mesh_parameters.element_area() =
        Mesh_Control_Parameters::Max_element_size;

      Bulk_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(
        triangle_mesh_parameters, this->time_stepper_pt());

      Bulk_mesh_pt->max_element_size() =
        Mesh_Control_Parameters::Max_element_size;
      Bulk_mesh_pt->min_element_size() =
        Mesh_Control_Parameters::Min_element_size;
      Bulk_mesh_pt->min_permitted_angle() =
        Mesh_Control_Parameters::Min_permitted_angle;
      Bulk_mesh_pt->max_keep_unrefined() = 400;
    }

    // Create the Outer_boundary_polyline_pt via one of the three methods
    // below. Used in the creation of the bulk elements.
    void create_rectangle_domain()
    {
      // Halfwidth of domain
      double half_width = 1.0;

      // Domain height
      double domain_height = Volume / pow(half_width, 2.0) * 2.0;

      // Build the boundary segments for outer boundary, consisting of
      //--------------------------------------------------------------
      // four separate polylines
      //------------------------
      Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);

      // Each polyline only has two vertices -- provide storage for their
      // coordinates
      Vector<Vector<double>> vertex_coord(2);
      for (unsigned i = 0; i < 2; i++)
      {
        vertex_coord[i].resize(2);
      }

      // First polyline: Free_surface_boundary_id
      vertex_coord[0][0] = 0.0;
      vertex_coord[0][1] = 0.0;
      vertex_coord[1][0] = half_width;
      vertex_coord[1][1] = 0.0;

      // Build the 1st boundary polyline
      boundary_polyline_pt[0] =
        new TriangleMeshPolyLine(vertex_coord, Free_surface_boundary_id);

      // Second boundary polyline: Outer wall with slip
      vertex_coord[0][0] = vertex_coord[1][0];
      vertex_coord[0][1] = vertex_coord[1][1];
      vertex_coord[1][0] = half_width;
      vertex_coord[1][1] = domain_height;

      // Build the 2nd boundary polyline
      boundary_polyline_pt[1] =
        new TriangleMeshPolyLine(vertex_coord, Outer_boundary_with_slip_id);

      // Third boundary polyline: Outflow
      vertex_coord[0][0] = vertex_coord[1][0];
      vertex_coord[0][1] = vertex_coord[1][1];
      vertex_coord[1][0] = 0.0;
      vertex_coord[1][1] = domain_height;

      // Build the 3rd boundary polyline
      boundary_polyline_pt[2] =
        new TriangleMeshPolyLine(vertex_coord, Upper_boundary_id);

      // Fourth boundary polyline: Bottom wall
      vertex_coord[0][0] = vertex_coord[1][0];
      vertex_coord[0][1] = vertex_coord[1][1];
      vertex_coord[1][0] = 0.0;
      vertex_coord[1][1] = 0.0;

      // Build the 4th boundary polyline
      boundary_polyline_pt[3] =
        new TriangleMeshPolyLine(vertex_coord, Inner_boundary_id);

      // Set max length for lower boundary
      boundary_polyline_pt[0]->set_maximum_length(
        Mesh_Control_Parameters::Max_free_surface_polyline_length);

      boundary_polyline_pt[1]->set_maximum_length(
        Mesh_Control_Parameters::Max_slip_polyline_length);

      // Create the triangle mesh polygon for outer boundary
      Outer_boundary_polyline_pt =
        new TriangleMeshPolygon(boundary_polyline_pt);
    }

    void create_circular_domain()
    {
      // Halfwidth of domain
      double half_width = 1.0;

      // Domain height
      double domain_height = Volume / pow(half_width, 2.0) * 2.0;

      // Number of points to use for the free surface polyline
      const unsigned npoints =
        Mesh_Control_Parameters::initial_number_of_free_surface_points;

      double radius = 1.0 / (cos(Contact_angle));
      double zeta_step =
        (0.5 * MathematicalConstants::Pi - Contact_angle) / double(npoints - 1);

      // Shift surface to ensure the volume is conserved
      double shift =
        2.0 / 3.0 *
        ((1 - 2 * ((radius) < 0)) * pow(pow(radius, 2.0) - 1, 3.0 / 2.0) -
         pow(radius, 3.0));

      double x_center = 0;
      double y_center = shift;

      Circle* circle_pt = new Circle(x_center, y_center, radius);

      // Intrinsic coordinate along GeomObject defining the bubble
      Vector<double> zeta(1);
      // Position vector to GeomObject defining the bubble
      Vector<double> coord(2);

      // Build the boundary segments for outer boundary, consisting of
      //--------------------------------------------------------------
      // four separate polylines
      //------------------------
      Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);

      Vector<Vector<double>> free_surface_vertex_coord;
      Vector<double> current_vertex(2);
      for (unsigned ipoint = 0; ipoint < npoints; ipoint++)
      {
        zeta[0] = 0.5 * MathematicalConstants::Pi - zeta_step * double(ipoint);
        circle_pt->position(zeta, current_vertex);
        free_surface_vertex_coord.push_back(current_vertex);
      }

      // Build the 1st boundary polyline
      boundary_polyline_pt[0] = new TriangleMeshPolyLine(
        free_surface_vertex_coord, Free_surface_boundary_id);

      // The rest of the polylines only have two vertices -- provide storage
      // for their coordinates
      Vector<Vector<double>> vertex_coord(2);
      for (unsigned ipoint = 0; ipoint < 2; ipoint++)
      {
        vertex_coord[ipoint].resize(2);
      }

      // Second boundary polyline: Outer wall with slip
      vertex_coord[0] = free_surface_vertex_coord.back();
      vertex_coord[1][0] = half_width;
      vertex_coord[1][1] = domain_height;

      // Build the 2nd boundary polyline
      boundary_polyline_pt[1] =
        new TriangleMeshPolyLine(vertex_coord, Outer_boundary_with_slip_id);

      // Fourth boundary polyline: Upper boundary
      vertex_coord[0][0] = vertex_coord[1][0];
      vertex_coord[0][1] = vertex_coord[1][1];
      vertex_coord[1][0] = 0.0;
      vertex_coord[1][1] = domain_height;

      // Build the 4rd boundary polyline
      boundary_polyline_pt[2] =
        new TriangleMeshPolyLine(vertex_coord, Upper_boundary_id);

      // Fifth boundary polyline: Inner wall
      vertex_coord[0][0] = vertex_coord[1][0];
      vertex_coord[0][1] = vertex_coord[1][1];
      vertex_coord[1] = free_surface_vertex_coord.front();

      // Build the 5th boundary polyline
      boundary_polyline_pt[3] =
        new TriangleMeshPolyLine(vertex_coord, Inner_boundary_id);

      // Set max length for lower boundary
      boundary_polyline_pt[0]->set_maximum_length(
        Mesh_Control_Parameters::Max_free_surface_polyline_length);

      boundary_polyline_pt[1]->set_maximum_length(
        Mesh_Control_Parameters::Max_slip_polyline_length);

      // Create the triangle mesh polygon for outer boundary
      Outer_boundary_polyline_pt =
        new TriangleMeshPolygon(boundary_polyline_pt);
    }

    // Create the free surface elements
    void create_free_surface_elements()
    {
      // Loop over the free surface boundary and create the "interface
      // elements
      unsigned b = Free_surface_boundary_id;

      // How many bulk fluid elements are adjacent to boundary b?
      unsigned n_element = Bulk_mesh_pt->nboundary_element(b);

      // Loop over the bulk fluid elements adjacent to boundary b?
      for (unsigned e = 0; e < n_element; e++)
      {
        // Get pointer to the bulk fluid element that is
        // adjacent to boundary b
        ELEMENT* bulk_elem_pt =
          dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b, e));

        // Find the index of the face of element e along boundary b
        int face_index = Bulk_mesh_pt->face_index_at_boundary(b, e);

        // Create new element
        FREE_SURFACE_ELEMENT* el_pt =
          new FREE_SURFACE_ELEMENT(bulk_elem_pt,
                                   face_index,
                                   this->time_stepper_pt(),
                                   Lagrange_id::Kinematic);

        // Add the appropriate boundary number
        el_pt->set_boundary_number_in_bulk_mesh(b);

        // Add the capillary number
        el_pt->ca_pt() = &Ca;
        el_pt->st_pt() = &St;

        // Add the external pressure data
        el_pt->set_external_pressure_data(External_pressure_data_pt);

        // Add it to the mesh
        Free_surface_mesh_pt->add_element_pt(el_pt);
      }
    }
    // Create the volume constraint elements
    void create_volume_constraint_elements()
    {
      // Build the single volume constraint element
      if (Volume_constraint_mesh_pt)
      {
        VOLUME_CONSTRAINT_ELEMENT* vol_constraint_element =
          new VOLUME_CONSTRAINT_ELEMENT(&Volume, External_pressure_data_pt, 0);
        Volume_constraint_mesh_pt->add_element_pt(vol_constraint_element);

        if (Volume_computation_mesh_pt)
        {
          // Loop over all boundaries
          for (unsigned b = 0; b < Bulk_mesh_pt->nboundary(); b++)
          {
            // How many bulk fluid elements are adjacent to boundary b?
            unsigned n_element = Bulk_mesh_pt->nboundary_element(b);

            // Loop over the bulk fluid elements adjacent to boundary b?
            for (unsigned e = 0; e < n_element; e++)
            {
              // Get pointer to the bulk fluid element that is
              // adjacent to boundary b
              ELEMENT* bulk_elem_pt =
                dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b, e));

              // Find the index of the face of element e along boundary b
              int face_index = Bulk_mesh_pt->face_index_at_boundary(b, e);

              // Create new element
              VOLUME_COMPUTATION_ELEMENT* el_pt =
                new VOLUME_COMPUTATION_ELEMENT(bulk_elem_pt, face_index);

              // Set the "master" volume control element
              el_pt->set_volume_constraint_element(vol_constraint_element);

              // Add it to the mesh
              Volume_computation_mesh_pt->add_element_pt(el_pt);
            }
          }
        }
      }
    } // end_of_create_volume_constraint_elements


    // Create the contact angle element
    void create_contact_angle_element()
    {
      // Find the element and node at the end of the free surface which meets
      // the wall

      // Inialise storage for bounding element
      FluidInterfaceBoundingElement* el_pt = 0;

      // Here we have no guarantee of order so we need to loop over all
      // surface elements to find the one that is next to the outer boundary
      unsigned n_free_surface = Free_surface_mesh_pt->nelement();
      for (unsigned e = 0; e < n_free_surface; e++)
      {
        // Locally cache the element pointer
        FREE_SURFACE_ELEMENT* bulk_el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(e));

        // Read out number of nodes in the element
        unsigned n_node = bulk_el_pt->nnode();

        // Is the "left" hand node on the boundary
        if (bulk_el_pt->node_pt(0)->is_on_boundary(Outer_boundary_with_slip_id))
        {
          // Create bounding element on "left" hand face, with the normal
          // pointing into the fluid
          el_pt = bulk_el_pt->make_bounding_element(-1);

          // Exit loop
          break;
        }

        // Is the "right" hand node on the boundary
        if (bulk_el_pt->node_pt(n_node - 1)
              ->is_on_boundary(Outer_boundary_with_slip_id))
        {
          // Create bounding element on "right" hand face, with the normal
          // pointing out of the fluid
          el_pt = bulk_el_pt->make_bounding_element(1);

          // Exit loop
          break;
        }
      }

      // Set the contact angle function
      el_pt->set_contact_angle(
        &Contact_angle, Global_Physical_Parameters::Use_strong_imposition);

      // Set the capillary number
      el_pt->ca_pt() = &Ca;

      // Set the wall normal of the external boundary
      el_pt->wall_unit_normal_fct_pt() =
        &Global_Physical_Parameters::wall_unit_normal_fct;

      // Add the element to the mesh
      Contact_angle_mesh_pt->add_element_pt(el_pt);

      // Here we have no guarantee of order so we need to loop over all
      // surface elements to find the one that is next to the outer boundary
      for (unsigned e = 0; e < n_free_surface; e++)
      {
        // Locally cache the element pointer
        FREE_SURFACE_ELEMENT* bulk_el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(e));

        // Read out number of nodes in the element
        unsigned n_node = bulk_el_pt->nnode();

        // Is the "left" hand node on the boundary
        if (bulk_el_pt->node_pt(0)->is_on_boundary(Inner_boundary_id))
        {
          // Create bounding element on "left" hand face, with the normal
          // pointing into the fluid
          el_pt = bulk_el_pt->make_bounding_element(-1);

          // Exit loop
          break;
        }

        // Is the "right" hand node on the boundary
        if (bulk_el_pt->node_pt(n_node - 1)->is_on_boundary(Inner_boundary_id))
        {
          // Create bounding element on "right" hand face, with the normal
          // pointing out of the fluid
          el_pt = bulk_el_pt->make_bounding_element(1);

          // Exit loop
          break;
        }
      }

      // Set the contact angle function
      el_pt->set_contact_angle(
        &Right_angle, Global_Physical_Parameters::Use_strong_imposition);

      // Set the capillary number
      el_pt->ca_pt() = &Ca;

      // Set sigma
      el_pt->sigma_pt() = &Zero_sigma;

      // Set the wall normal of the external boundary
      el_pt->wall_unit_normal_fct_pt() =
        &Global_Physical_Parameters::wall_unit_normal_fct;

      // Add the element to the mesh
      Contact_angle_mesh_pt->add_element_pt(el_pt);

    } // end_of_create_contact_angle_element


    // Create the slip boundary elements
    void create_slip_elements(const unsigned& b,
                              Mesh* const& bulk_mesh_pt,
                              Mesh* const& surface_mesh_pt)
    {
      // How many bulk elements are adjacent to boundary b?
      unsigned n_element = bulk_mesh_pt->nboundary_element(b);

      // Loop over the bulk elements adjacent to boundary b?
      for (unsigned e = 0; e < n_element; e++)
      {
        // Get pointer to the bulk element that is adjacent to boundary b
        ELEMENT* bulk_elem_pt =
          dynamic_cast<ELEMENT*>(bulk_mesh_pt->boundary_element_pt(b, e));

        // What is the index of the face of element e along boundary b
        int face_index = bulk_mesh_pt->face_index_at_boundary(b, e);
        SLIP_ELEMENT* slip_element_pt = 0;

        // Build the corresponding slip element
        slip_element_pt = new SLIP_ELEMENT(bulk_elem_pt, face_index);

        // Set the pointer to the prescribed slip function
        slip_element_pt->slip_fct_pt() = &Slip_Parameters::prescribed_slip_fct;

        slip_element_pt->wall_velocity_fct_pt() =
          &Slip_Parameters::prescribed_wall_velocity_fct;

        slip_element_pt->add_external_data(Wall_velocity_data_pt);

        // Add the prescribed-flux element to the surface mesh
        surface_mesh_pt->add_element_pt(slip_element_pt);

      } // end of loop over bulk elements adjacent to boundary b
    }

    // Create the flux elements
    void create_height_elements()
    {
      HEIGHT_ELEMENT* height_el_pt = new HEIGHT_ELEMENT(
        Inner_corner_solid_node_pt, Contact_line_solid_node_pt);

      height_el_pt->set_parameter_data_pt(Height_step_data_pt);

      Height_mesh_pt->add_element_pt(height_el_pt);
    }

    void create_slip_eigen_elements(Mesh* const& bulk_mesh_pt)
    {
      oomph_info << "create_slip_eigen_elements" << std::endl;

      // Loop over the free surface boundary and create the "interface elements
      unsigned b = Outer_boundary_with_slip_id;

      // How many bulk fluid elements are adjacent to boundary b?
      unsigned n_element = bulk_mesh_pt->nboundary_element(b);

      // Loop over the bulk fluid elements adjacent to boundary b?
      for (unsigned e = 0; e < n_element; e++)
      {
        // Get pointer to the bulk fluid element that is
        // adjacent to boundary b
        ELEMENT* bulk_elem_pt =
          dynamic_cast<ELEMENT*>(bulk_mesh_pt->boundary_element_pt(b, e));

        if (bulk_elem_pt->is_augmented())
        {
          // Find the index of the face of element e along boundary b
          int face_index = bulk_mesh_pt->face_index_at_boundary(b, e);

          // Create new element
          SingularNavierStokesTractionElement<ELEMENT>* el_pt =
            new SingularNavierStokesTractionElement<ELEMENT>(
              bulk_elem_pt,
              face_index,
              Singularity_scaling_mesh_pt->element_pt(0)->internal_data_pt(0));

          el_pt->traction_fct_pt() = &parameters::eigensolution_slip_fct;

          // Add it to the mesh
          Eigensolution_slip_mesh_pt->add_element_pt(el_pt);
        }
      }
    }

    void create_singularity_scaling_elements()
    {
      oomph_info << "create_singularity_scaling_elements" << std::endl;
      SingularNavierStokesSolutionElement<ELEMENT>* el_pt =
        new SingularNavierStokesSolutionElement<ELEMENT>;

      // Set the pointer to the velocity singular function for this
      // element, defined in parameters namespace
      el_pt->velocity_singular_fct_pt() = &parameters::velocity_singular_fct;

      // Set the pointer to the gradient of the velocity singular
      // function for this element, defined in parameters namespace
      el_pt->grad_velocity_singular_fct_pt() =
        &parameters::grad_velocity_singular_fct;

      // Set the pointer to the first pressure singular function for this
      // element, defined in parameters namespace
      el_pt->pressure_singular_fct_pt() = &parameters::pressure_singular_fct;

      // The singular function satisfies the Stokes equation
      el_pt->singular_function_satisfies_stokes_equation() = false;

      // el_pt->pin_c();
      el_pt->set_c(0.0);

      // Add element to the mesh
      Singularity_scaling_mesh_pt->add_element_pt(el_pt);
    }

    // Create the no penetration boundary elements
    void create_no_penetration_elements(const unsigned& b,
                                        Mesh* const& bulk_mesh_pt,
                                        Mesh* const& surface_mesh_pt)
    {
      // How many bulk elements are adjacent to boundary b?
      unsigned n_element = bulk_mesh_pt->nboundary_element(b);

      // Loop over the bulk elements adjacent to boundary b?
      for (unsigned e = 0; e < n_element; e++)
      {
        // Get pointer to the bulk element that is adjacent to boundary b
        ELEMENT* bulk_elem_pt =
          dynamic_cast<ELEMENT*>(bulk_mesh_pt->boundary_element_pt(b, e));

        // What is the index of the face of element e along boundary b
        int face_index = bulk_mesh_pt->face_index_at_boundary(b, e);

        // Build the corresponding slip element
        NO_PENETRATION_ELEMENT* no_penetration_element_pt =
          new NO_PENETRATION_ELEMENT(
            bulk_elem_pt, face_index, Lagrange_id::No_penetration);

        // Add the prescribed-flux element to the surface mesh
        surface_mesh_pt->add_element_pt(no_penetration_element_pt);

      } // end of loop over bulk elements adjacent to boundary b
    }

    // Create the flux elements
    void create_flux_elements()
    {
      if (Net_flux_mesh_pt)
      {
        // Build master element
        NET_FLUX_ELEMENT* net_flux_el_pt =
          new NET_FLUX_ELEMENT(&Flux_Parameters::flux_fct, time_pt());

        // Add NetFluxControlElement to its mesh
        Net_flux_mesh_pt->add_element_pt(net_flux_el_pt);

        if (Flux_mesh_pt)
        {
          // Loop over the free surface boundary and create the interface
          // elements
          unsigned b = Upper_boundary_id;

          // How many bulk fluid elements are adjacent to boundary b?
          unsigned n_element = Bulk_mesh_pt->nboundary_element(b);

          // Loop over the bulk fluid elements adjacent to boundary b?
          for (unsigned e = 0; e < n_element; e++)
          {
            // Get pointer to the bulk fluid element that is
            // adjacent to boundary b
            ELEMENT* bulk_elem_pt =
              dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b, e));

            // Find the index of the face of element e along boundary b
            int face_index = Bulk_mesh_pt->face_index_at_boundary(b, e);

            // Create new element
            FLUX_ELEMENT* el_pt = new FLUX_ELEMENT(
              bulk_elem_pt,
              face_index,
              net_flux_el_pt->get_lagrange_multiplier_data_pt());

            // Add it to the mesh
            Flux_mesh_pt->add_element_pt(el_pt);
          }
        }
      }
    }

    void setup_mesh_interaction()
    {
      oomph_info << "setup_mesh_interaction" << std::endl;
      // Find corner bulk element
      unsigned element_index;
      unsigned node_index;
      find_corner_bulk_element_and_node(Outer_boundary_with_slip_id,
                                        Free_surface_boundary_id,
                                        element_index,
                                        node_index);
      ELEMENT* corner_bulk_element_pt =
        dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(
          Outer_boundary_with_slip_id, element_index));

      // Tell the CEquationElement object about its associated Navier-Stokes
      // element, the component of the velocity whose derivative will be
      // computed in the residual, the local coordinate in the Navier_Stokes
      // element at which the residual will be computed, and the direction of
      // the derivative
      Vector<double> s1_pt(2);
      s1_pt[0] = 0.0;
      s1_pt[1] = 0.0;
      unsigned* const direction_pt = new unsigned(0);
      SingularNavierStokesSolutionElement<ELEMENT>* singular_el_pt =
        dynamic_cast<SingularNavierStokesSolutionElement<ELEMENT>*>(
          Singularity_scaling_mesh_pt->element_pt(0));
      singular_el_pt->set_wrapped_navier_stokes_element_pt(
        corner_bulk_element_pt, s1_pt, direction_pt);
      Vector<double> x(2, 0.0);
      corner_bulk_element_pt->get_x(s1_pt, x);
      oomph_info << "First singular element point, ";
      oomph_info << "x: " << x[0] << ", ";
      oomph_info << "y: " << x[1] << ", ";
      oomph_info << std::endl;

      // Loop over the augmented bulk elements
      unsigned n_aug_bulk = Augmented_bulk_element_number.size();
      for (unsigned e = 0; e < n_aug_bulk; e++)
      {
        // Augment elements
        // Upcast from GeneralisedElement to the present element
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
          Bulk_mesh_pt->element_pt(Augmented_bulk_element_number[e]));

        // Set the pointer to the element that determines the amplitude
        // of the singular fct
        el_pt->add_c_equation_element_pt(singular_el_pt);
      }
    }

    void setup_pressure_contribution_elements()
    {
      oomph_info << "setup_pressure_contribution_elements" << std::endl;
      // Locate zeta for the first element
      Vector<double> zeta(1, parameters::small_r);

      GeomObject* sub_geom_object_pt = 0;
      Vector<double> s(1);

      Pressure_contribution_geom_mesh_1_pt->locate_zeta(
        zeta, sub_geom_object_pt, s);

      dynamic_cast<PressureContributionElement<ELEMENT>*>(sub_geom_object_pt)
        ->set_evaluation_point(s);

      zeta[0] = 1 - parameters::small_r;

      Pressure_contribution_geom_mesh_2_pt->locate_zeta(
        zeta, sub_geom_object_pt, s);

      dynamic_cast<PressureContributionElement<ELEMENT>*>(sub_geom_object_pt)
        ->set_evaluation_point(s);
    }

    void create_pressure_contribution_1_elements()
    {
      oomph_info << "create_pressure_contribution_1_elements" << std::endl;
      const bool add_to_residual = true;

      // Loop over boundary elements
      const unsigned n_boundary_element =
        Bulk_mesh_pt->nboundary_element(Outer_boundary_with_slip_id);
      for (unsigned element_index = 0; element_index < n_boundary_element;
           element_index++)
      {
        // Find corner bulk element
        ELEMENT* bulk_element_pt =
          dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(
            Outer_boundary_with_slip_id, element_index));

        // Find the index of the face of element e along boundary b
        int face_index = Bulk_mesh_pt->face_index_at_boundary(
          Outer_boundary_with_slip_id, element_index);

        PressureContributionElement<ELEMENT>* el_pt =
          new PressureContributionElement<ELEMENT>(
            bulk_element_pt, face_index, add_to_residual);

        el_pt->set_pressure_data_pt(
          Singularity_scaling_mesh_pt->element_pt(0)->internal_data_pt(0));

        el_pt->set_boundary_number_in_bulk_mesh(Outer_boundary_with_slip_id);

        Pressure_contribution_mesh_1_pt->add_element_pt(el_pt);
      }
    }

    void create_pressure_contribution_2_elements()
    {
      oomph_info << "create_pressure_contribution_2_elements" << std::endl;
      const bool add_to_residual = false;

      // Loop over boundary elements
      const unsigned n_boundary_element =
        Bulk_mesh_pt->nboundary_element(Free_surface_boundary_id);
      for (unsigned element_index = 0; element_index < n_boundary_element;
           element_index++)
      {
        // Find corner bulk element
        ELEMENT* bulk_element_pt =
          dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(
            Free_surface_boundary_id, element_index));

        // Find the index of the face of element e along boundary b
        int face_index = Bulk_mesh_pt->face_index_at_boundary(
          Free_surface_boundary_id, element_index);

        PressureContributionElement<ELEMENT>* el_pt =
          new PressureContributionElement<ELEMENT>(
            bulk_element_pt, face_index, add_to_residual);

        el_pt->set_pressure_data_pt(
          Singularity_scaling_mesh_pt->element_pt(0)->internal_data_pt(0));

        el_pt->set_boundary_number_in_bulk_mesh(Free_surface_boundary_id);

        Pressure_contribution_mesh_2_pt->add_element_pt(el_pt);
      }
    }

    void create_mesh_as_geom_object()
    {
      Pressure_contribution_geom_mesh_1_pt =
        new MeshAsGeomObject(Pressure_contribution_mesh_1_pt);

      Pressure_contribution_geom_mesh_2_pt =
        new MeshAsGeomObject(Pressure_contribution_mesh_2_pt);
    }


    //============================================================================
    // Usage functions
    //============================================================================
  public:
    void update_triangulateio()
    {
      dynamic_cast<TriangleMesh<ELEMENT>*>(Bulk_mesh_pt)
        ->update_triangulateio();
    }

    TriangulateIO triangulateio_representation_of_mesh()
    {
      dynamic_cast<TriangleMesh<ELEMENT>*>(Bulk_mesh_pt)
        ->update_triangulateio();
      return dynamic_cast<TriangleMesh<ELEMENT>*>(Bulk_mesh_pt)
        ->triangulateio_representation();
    }

    void check_mass_matrix()
    {
      oomph_info << "check_mass_matrix" << std::endl;
      {
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0));
        unsigned n = el_pt->ndof();
        Vector<double> residuals(n, 0.0);
        DenseMatrix<double> jacobian(n, n, 0.0);
        DenseMatrix<double> mass_matrix(n, n, 0.0);
        el_pt->get_jacobian_and_mass_matrix(residuals, jacobian, mass_matrix);
      }
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(0));
        unsigned n = el_pt->ndof();
        Vector<double> residuals(n, 0.0);
        DenseMatrix<double> jacobian(n, n, 0.0);
        DenseMatrix<double> mass_matrix(n, n, 0.0);
        el_pt->get_jacobian_and_mass_matrix(residuals, jacobian, mass_matrix);
      }
      {
        SLIP_ELEMENT* el_pt =
          dynamic_cast<SLIP_ELEMENT*>(Slip_boundary_mesh_pt->element_pt(0));
        unsigned n = el_pt->ndof();
        Vector<double> residuals(n, 0.0);
        DenseMatrix<double> jacobian(n, n, 0.0);
        DenseMatrix<double> mass_matrix(n, n, 0.0);
        el_pt->get_jacobian_and_mass_matrix(residuals, jacobian, mass_matrix);
      }
      if (No_penetration_boundary_mesh_pt)
      {
        NO_PENETRATION_ELEMENT* el_pt = dynamic_cast<NO_PENETRATION_ELEMENT*>(
          No_penetration_boundary_mesh_pt->element_pt(0));
        unsigned n = el_pt->ndof();
        Vector<double> residuals(n, 0.0);
        DenseMatrix<double> jacobian(n, n, 0.0);
        DenseMatrix<double> mass_matrix(n, n, 0.0);
        el_pt->get_jacobian_and_mass_matrix(residuals, jacobian, mass_matrix);
      }
      {
        if (Flux_mesh_pt)
        {
          FLUX_ELEMENT* el_pt =
            dynamic_cast<FLUX_ELEMENT*>(Flux_mesh_pt->element_pt(0));
          unsigned n = el_pt->ndof();
          Vector<double> residuals(n, 0.0);
          DenseMatrix<double> jacobian(n, n, 0.0);
          DenseMatrix<double> mass_matrix(n, n, 0.0);
          el_pt->get_jacobian_and_mass_matrix(residuals, jacobian, mass_matrix);
        }
      }
      {
        VOLUME_COMPUTATION_ELEMENT* el_pt =
          dynamic_cast<VOLUME_COMPUTATION_ELEMENT*>(
            Volume_computation_mesh_pt->element_pt(0));
        unsigned n = el_pt->ndof();
        Vector<double> residuals(n, 0.0);
        DenseMatrix<double> jacobian(n, n, 0.0);
        DenseMatrix<double> mass_matrix(n, n, 0.0);
        el_pt->get_jacobian_and_mass_matrix(residuals, jacobian, mass_matrix);
      }
      {
        VOLUME_CONSTRAINT_ELEMENT* el_pt =
          dynamic_cast<VOLUME_CONSTRAINT_ELEMENT*>(
            Volume_constraint_mesh_pt->element_pt(0));
        unsigned n = el_pt->ndof();
        Vector<double> residuals(n, 0.0);
        DenseMatrix<double> jacobian(n, n, 0.0);
        DenseMatrix<double> mass_matrix(n, n, 0.0);
        el_pt->get_jacobian_and_mass_matrix(residuals, jacobian, mass_matrix);
      }
      {
        if (Net_flux_mesh_pt)
        {
          NET_FLUX_ELEMENT* el_pt =
            dynamic_cast<NET_FLUX_ELEMENT*>(Net_flux_mesh_pt->element_pt(0));
          unsigned n = el_pt->ndof();
          Vector<double> residuals(n, 0.0);
          DenseMatrix<double> jacobian(n, n, 0.0);
          DenseMatrix<double> mass_matrix(n, n, 0.0);
          el_pt->get_jacobian_and_mass_matrix(residuals, jacobian, mass_matrix);
        }
      }
    }


    // Refine the mesh for the weak contact angle constraint
    void refine_mesh_for_weak_contact_angle_constraint()
    {
      oomph_info << "Refining the mesh about the contact line" << std::endl;

      double max_error = 1e6;
      // double old_max_error = max_error + 1.0;
      double min_error = 0.0;
      const unsigned n_adapt =
        Mesh_Control_Parameters::Max_number_of_adapts_for_refinement;
      // double tol = 1e-3;
      unsigned i_adapt = 0;
      // Call the adapt function until the maximum of the estimated error is
      // below 1e0 ( = Bulk_mesh->max_permitted_error)
      while (max_error > 1e0 && i_adapt < n_adapt)
      //&& abs(old_max_error - max_error) < tol)
      {
        set_error_estimator();

        // old_max_error = max_error;
        compute_error_estimate(max_error, min_error);
        oomph_info << "Number of adaptions: " << i_adapt
                   << ", Max error: " << max_error << std::endl;

        unsigned n_elements = Bulk_mesh_pt->nelement();
        Vector<double> elemental_error(n_elements);
        Mesh* fluid_mesh_pt = dynamic_cast<Mesh*>(Bulk_mesh_pt);
        Bulk_mesh_pt->spatial_error_estimator_pt()->get_element_errors(
          fluid_mesh_pt, elemental_error);
        Bulk_mesh_pt->adapt(elemental_error);
        Bulk_mesh_pt->output("RESLT/mesh" + to_string(i_adapt) + ".dat");
        i_adapt++;
      }
    }

    // Solve for the initial state
    void solve_for_initial_state()
    {
      oomph_info << "Solving for initial state" << std::endl;

      // Ensure the problem is static
      if (!this->Is_steady)
      {
        make_steady();
      }

      // Solve the steady problem
      this->steady_newton_solve(Max_adapt);
    }

    // A custom adaptive steady newton solve function.
    // Uses the function is_adaption_needed() to determine if we should adapt
    // rather than just adapting until max_adapt or there is no elements that
    // need refining/unrefining.
    int steady_newton_solve_adapt_if_needed(const unsigned& max_adapt)
    {
      if (!this->Is_steady)
      {
        OomphLibError("Problem is not steady. Call make_steady before steady "
                      "newton solve.",
                      "steady_newton_solve_adapt_if_needed",
                      OOMPH_EXCEPTION_LOCATION);
      }
      // Loop count
      unsigned n = 0;
      bool local_is_adaption_needed = false;
      do
      {
        // If adaption is needed, adapt
        if (local_is_adaption_needed)
        {
          adapt();
          n++;
        }
        // Solve steady problem
        steady_newton_solve(0);
        local_is_adaption_needed = is_adaption_needed();
        // Increment loop count
      } while (n < max_adapt && local_is_adaption_needed);
      int is_solved = 0;
      if (!local_is_adaption_needed)
      {
        is_solved = 1;
      }
      return is_solved;
    }

    void pin_wall_velocity()
    {
      Wall_velocity_data_pt->pin(0);
    }

    void unpin_wall_velocity()
    {
      Wall_velocity_data_pt->unpin(0);
    }

    void pin_bond_number()
    {
      ReInvFr_data_pt->pin(0);
    }

    void unpin_bond_number()
    {
      ReInvFr_data_pt->unpin(0);
    }

    void set_height_step_parameter_to_bond_number()
    {
      Height_step_data_pt = ReInvFr_data_pt;
    }

    void set_height_step_parameter_to_wall_velocity()
    {
      Height_step_data_pt = Wall_velocity_data_pt;
    }

    // Takes a continuation step using a fixed height drop and solving for the
    // parameter.
    double height_step_solve(const double& ds)
    {
      pre_height_solve(ds);

      // this->debug_jacobian();
      this->steady_newton_solve();

      post_height_solve();
      return ds;
    }

    void pre_height_solve(const double& ds)
    {
      dynamic_cast<HEIGHT_ELEMENT*>(Height_mesh_pt->element_pt(0))
        ->step_height(ds);
      dynamic_cast<HEIGHT_ELEMENT*>(Height_mesh_pt->element_pt(0))
        ->pin_height();
      Height_step_data_pt->unpin(0);

      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }
    void post_height_solve()
    {
      dynamic_cast<HEIGHT_ELEMENT*>(Height_mesh_pt->element_pt(0))
        ->unpin_height();
      Height_step_data_pt->pin(0);

      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }


    // Check whether we need to adapt.
    //
    bool is_adaption_needed()
    {
      // Check if the Z2 error estimator is ok.
      double local_max_z2_error = 0;
      double local_min_z2_error = 0;
      get_z2_error(local_max_z2_error, local_min_z2_error);

      // If the Z2 error is not within tolerance
      if (local_max_z2_error > Mesh_Control_Parameters::Max_permitted_z2_error)
      {
        // Then adapt is needed
        oomph_info << "Adapt is needed due to Z2 error" << std::endl;
        return true;
      }

      // if the error estimator flag is 1 or 2
      // Check if the contact angle error is ok.
      if (Contact_angle_mesh_pt)
      {
        if (Mesh_Control_Parameters::Error_estimator_flag > 0)
        {
          double expected_contact_angle = 0.0;
          double actual_contact_angle = 0.0;
          dynamic_cast<CONTACT_LINE_ELEMENT*>(
            Contact_angle_mesh_pt->element_pt(0))
            ->calculate_contact_angle(expected_contact_angle,
                                      actual_contact_angle);

          oomph_info << "Contact angle error: "
                     << abs(expected_contact_angle - actual_contact_angle) *
                          180.0 / MathematicalConstants::Pi
                     << std::endl;
          // If the contact angle error is not within tolerance
          if (abs(expected_contact_angle - actual_contact_angle) * 180.0 /
                MathematicalConstants::Pi >
              1.0)
          {
            // If we are using the contact angle error estimator
            if (Using_contact_angle_error_estimator)
            {
              // We need to increase the resolution also
              Mesh_Control_Parameters::min_element_length *= 0.5;
            }

            // Then adapt is needed
            oomph_info << "Adapt is needed due to contact angle error" << std::endl;
            return true;
          }
        }
      }

      // if the error estimator flag is 2
      // Check if the inner angle error is ok.
      if (Contact_angle_mesh_pt)
      {
        if (Mesh_Control_Parameters::Error_estimator_flag > 0)
        {
          double expected_contact_angle = 0.0;
          double actual_contact_angle = 0.0;
          dynamic_cast<CONTACT_LINE_ELEMENT*>(
            Contact_angle_mesh_pt->element_pt(1))
            ->calculate_contact_angle(expected_contact_angle,
                                      actual_contact_angle);
          oomph_info << "Inner corner error: "
                     << abs(expected_contact_angle - actual_contact_angle) *
                          180.0 / MathematicalConstants::Pi
                     << std::endl;
          // If the contact angle error is not within tolerance
          if (abs(expected_contact_angle - actual_contact_angle) * 180.0 /
                MathematicalConstants::Pi >
              1.0)
          {
            // If we are using the contact angle error estimator
            if (Using_contact_angle_error_estimator)
            {
              // We need to increase the resolution also
              Mesh_Control_Parameters::inner_min_element_length *= 0.5;
            }

            // Then adapt is needed
            oomph_info << "Adapt is needed due to inner angle error" << std::endl;
            return true;
          }
        }
      }

      // If the max Z2 error is much smaller than the permitted, then adapt to
      // unrefine
      // if (local_max_z2_error /
      // Mesh_Control_Parameters::Max_permitted_z2_error <
      //     1e-1)
      // {
      //   // Then adapt is needed
      //   oomph_info << "Adapt is needed due to over-refining of the Z2
      //   error"
      //              << std::endl;
      //   return true;
      // }


      oomph_info << "No further adaption is needed" << std::endl;
      return false;
    }

    double get_contact_angle_error()
    {
      double expected_contact_angle = 0.0;
      double actual_contact_angle = 0.0;
      dynamic_cast<CONTACT_LINE_ELEMENT*>(Contact_angle_mesh_pt->element_pt(0))
        ->calculate_contact_angle(expected_contact_angle, actual_contact_angle);
      return abs(expected_contact_angle - actual_contact_angle);
    }

    void pin_contact_line()
    {
      oomph_info << "Pin contact line" << std::endl;

      Contact_line_solid_node_pt->pin_position(0);
      Contact_line_solid_node_pt->pin_position(1);
      Contact_line_solid_node_pt->pin(0);
      Contact_line_solid_node_pt->pin(1);
      Contact_line_solid_node_pt->pin(2);

      fix_contact_line_lagrange_multiplier();
      fix_cl_no_penetration_lagrange_multiplier();

      // Rebuild the global mesh
      this->rebuild_global_mesh();

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void pin_centre()
    {
      oomph_info << "Pin centre" << std::endl;

      Inner_corner_solid_node_pt->pin_position(0);
      Inner_corner_solid_node_pt->pin_position(1);
      Inner_corner_solid_node_pt->pin(0);
      Inner_corner_solid_node_pt->pin(1);
      Inner_corner_solid_node_pt->pin(2);

      fix_centre_lagrange_multiplier();

      // Rebuild the global mesh
      this->rebuild_global_mesh();

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void pin_centre_element()
    {
      oomph_info << "Pin centre element" << std::endl;
      unsigned n_element;
      n_element = Free_surface_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 3; m++)
        {
          if (el_pt->node_pt(m)->is_on_boundary(Inner_boundary_id))
          {
            for (unsigned p = 0; p < 3; p++)
            {
              el_pt->pin_lagrange_multiplier(p);
              // Pin RC and RS
              dynamic_cast<SolidNode*>(el_pt->node_pt(p))->pin_position(0);
              dynamic_cast<SolidNode*>(el_pt->node_pt(p))->pin_position(1);

              el_pt->node_pt(p)->pin(0);
              el_pt->node_pt(p)->pin(1);
              el_pt->node_pt(p)->pin(2);
            }
          }
        }
      }

      // Rebuild the global mesh
      this->rebuild_global_mesh();

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void pin_centre_half()
    {
      oomph_info << "Pin position" << std::endl;
      unsigned n_element;
      n_element = Free_surface_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 3; m++)
        {
          if (el_pt->node_pt(m)->x(0) <= 0.5)
          {
            el_pt->pin_lagrange_multiplier(m);
            // Pin RC and RS
            dynamic_cast<SolidNode*>(el_pt->node_pt(m))->pin_position(0);
            dynamic_cast<SolidNode*>(el_pt->node_pt(m))->pin_position(1);

            el_pt->node_pt(m)->pin(0);
            el_pt->node_pt(m)->pin(1);
            el_pt->node_pt(m)->pin(2);
          }
        }
      }

      // Rebuild the global mesh
      this->rebuild_global_mesh();

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void free_contact_line()
    {
      // Node* node_pt = 0;
      // find_corner_node(
      //  Outer_boundary_with_slip_id, Free_surface_boundary_id, node_pt);

      // dynamic_cast<SolidNode*>(node_pt)->unpin_position(0);
      // dynamic_cast<SolidNode*>(node_pt)->unpin_position(1);
      Contact_line_solid_node_pt->unpin_position(0);
      Contact_line_solid_node_pt->unpin_position(1);
    }

    void pin_kinematic_lagrange_multiplier()
    {
      // Pin all kinematic lagrange multipliers
      const unsigned n_el = Free_surface_mesh_pt->nelement();
      for (unsigned i_el = 0; i_el < n_el; i_el++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(i_el));
        const unsigned n_nod = el_pt->nnode();
        for (unsigned i_nod = 0; i_nod < n_nod; i_nod++)
        {
          // Get boundary node
          el_pt->pin_lagrange_multiplier(i_nod);
        }
      }

      // Set up the equation numbering so we are ready to solve the problem.
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void free_kinematic_lagrange_multiplier()
    {
      const unsigned n_el = Free_surface_mesh_pt->nelement();
      for (unsigned i_el = 0; i_el < n_el; i_el++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(i_el));
        el_pt->free();
      }
    }

    void set_kinematic_lagrange_multiplier(const double& value)
    {
      const unsigned n_el = Free_surface_mesh_pt->nelement();
      for (unsigned i_el = 0; i_el < n_el; i_el++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(i_el));
        const unsigned n_nod = el_pt->nnode();
        for (unsigned i_nod = 0; i_nod < n_nod; i_nod++)
        {
          el_pt->set_lagrange_multiplier(i_nod, value);
        }
      }
    }

    void pin_no_penetration_lagrange_multiplier(const double& value)
    {
      if (No_penetration_boundary_mesh_pt)
      {
        const unsigned n_el = No_penetration_boundary_mesh_pt->nelement();
        for (unsigned i_el = 0; i_el < n_el; i_el++)
        {
          NO_PENETRATION_ELEMENT* no_penetration_element_pt =
            dynamic_cast<NO_PENETRATION_ELEMENT*>(
              No_penetration_boundary_mesh_pt->element_pt(i_el));
          const unsigned n_nod = no_penetration_element_pt->nnode();
          for (unsigned i_nod = 0; i_nod < n_nod; i_nod++)
          {
            no_penetration_element_pt->pin_lagrange_multiplier(i_nod);
          }
        }
      }
    }

    void fix_contact_line_lagrange_multiplier()
    {
      // Pin all kinematic lagrange multipliers
      const unsigned n_el = Free_surface_mesh_pt->nelement();
      for (unsigned i_el = 0; i_el < n_el; i_el++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(i_el));
        const unsigned n_nod = el_pt->nnode();
        for (unsigned i_nod = 0; i_nod < n_nod; i_nod++)
        {
          if (el_pt->node_pt(i_nod)->is_on_boundary(
                Outer_boundary_with_slip_id))
          {
            // Get boundary node
            el_pt->pin_lagrange_multiplier(i_nod);
          }
        }
      }
    }

    void fix_centre_lagrange_multiplier()
    {
      // Pin all kinematic lagrange multipliers
      const unsigned n_el = Free_surface_mesh_pt->nelement();
      for (unsigned i_el = 0; i_el < n_el; i_el++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(i_el));
        const unsigned n_nod = el_pt->nnode();
        for (unsigned i_nod = 0; i_nod < n_nod; i_nod++)
        {
          if (el_pt->node_pt(i_nod)->is_on_boundary(Inner_boundary_id))
          {
            // Get boundary node
            el_pt->pin_lagrange_multiplier(i_nod);
          }
        }
      }
    }

    void fix_cl_no_penetration_lagrange_multiplier()
    {
      if (No_penetration_boundary_mesh_pt)
      {
        const unsigned n_el = No_penetration_boundary_mesh_pt->nelement();
        for (unsigned i_el = 0; i_el < n_el; i_el++)
        {
          NO_PENETRATION_ELEMENT* no_penetration_element_pt =
            dynamic_cast<NO_PENETRATION_ELEMENT*>(
              No_penetration_boundary_mesh_pt->element_pt(i_el));
          const unsigned n_nod = no_penetration_element_pt->nnode();
          for (unsigned i_nod = 0; i_nod < n_nod; i_nod++)
          {
            if (no_penetration_element_pt->node_pt(i_nod)->is_on_boundary(
                  Free_surface_boundary_id))
            {
              no_penetration_element_pt->pin_lagrange_multiplier(i_nod);
            }
          }
        }
      }
    }

    void find_corner_bulk_element_and_node(const unsigned& boundary_1_id,
                                           const unsigned& boundary_2_id,
                                           unsigned& element_index,
                                           unsigned& node_index)
    {
      unsigned n_boundary_element =
        Bulk_mesh_pt->nboundary_element(boundary_1_id);
      for (unsigned e = 0; e < n_boundary_element; e++)
      {
        // Locally cache the element pointer
        FiniteElement* bulk_el_pt =
          Bulk_mesh_pt->boundary_element_pt(boundary_1_id, e);

        // Read out number of nodes in the element
        unsigned n_node = bulk_el_pt->nnode();
        for (unsigned i_node = 0; i_node < n_node; i_node++)
        {
          // If the node is on the free surface boundary as well then ...
          if (bulk_el_pt->node_pt(i_node)->is_on_boundary(boundary_2_id) &&
              bulk_el_pt->node_pt(i_node)->is_on_boundary(boundary_1_id))
          {
            // set the output arguments,
            element_index = e;
            node_index = i_node;

            // Return to exit both loops and end function
            return;
          }
        }
      }
      // If not found, issue warning and return anyway
      oomph_info << "Warning: No corner node found!" << std::endl;

      return;
    }

    void restrict_horizontal_displacement()
    {
      unsigned n_node = Bulk_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        static_cast<SolidNode*>(Bulk_mesh_pt->node_pt(n))->pin_position(0);
      }
    }

    void restrict_vertical_displacement()
    {
      unsigned n_node = Bulk_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        static_cast<SolidNode*>(Bulk_mesh_pt->node_pt(n))->pin_position(1);
      }
    }

    // Ensure the problem is dynamic and time step forward using 'n_tsteps'
    // steps of size 'dt'
    void timestep(const double& dt, const double& ft)
    {
      // Ensure the problem is dynamic
      if (this->Is_steady)
      {
        make_unsteady();
      }

      // Set up variables for the adaption procedure
      unsigned steps_between_adapt =
        Mesh_Control_Parameters::interval_between_adapts;
      unsigned local_max_adapt = 0;

      // Set up variables for unsteady_newton_solve
      double local_dt = dt;
      const double temporal_tolerance =
        Mesh_Control_Parameters::Temporal_tolerance;
      bool first_timestep = false;
      bool shift = true;
      unsigned it = 0;
      double final_time_tolerance = 1e-6;
      while (this->time_pt()->time() < ft - final_time_tolerance)
      {
        oomph_info << "t: " << this->time() << std::endl;

        // If (the step counter + 1) is a multiple of the number of steps
        // allowed between adapts, then ...
        if (it % steps_between_adapt == 0)
        {
          // ... Set the local_max_adapt to the global one ...
          local_max_adapt = this->Max_adapt;
          // if (it == 0)
          //{
          //  local_max_adapt = 5;
          //}
        }
        else
        {
          // ..otherwise don't allow adapting.
          local_max_adapt = 0;
        }

        // Call unsteady newton solver
        try
        {
          if (Mesh_Control_Parameters::Use_adaptive_timestepping)
          {
            local_dt =
              this->doubly_adaptive_unsteady_newton_solve(local_dt,
                                                          temporal_tolerance,
                                                          local_max_adapt,
                                                          first_timestep,
                                                          shift);
          }
          else
          {
            this->unsteady_newton_solve(
              local_dt, local_max_adapt, first_timestep, shift);
          }
        }
        catch (std::runtime_error& err)
        {
          oomph_info << "Caught error" << std::endl;
          break;
        }

        // Document the solution
        if ((it + 1) % Doc_Parameters::interval_between_doc == 0)
        {
          this->doc_solution();
        }

        // Dump the solution
        if ((it + 1) % Restart_Parameters::interval_between_dump == 0)
        {
          this->create_restart_file();
        }
        it++;
      }
    }

    // Make the problem static. Remove the flux elements and add the volume
    // ones, by calling actions before and after adapt with the Is_steady flag
    // set to true
    void make_steady()
    {
      oomph_info << "make_steady" << std::endl;

      actions_before_adapt();

      this->Is_steady = true;

      actions_after_adapt();
    }

    // Make the problem static. Add the flux elements and remove the volume
    // ones, by calling actions before and after adapt with the Is_steady flag
    // set to false
    void make_unsteady()
    {
      oomph_info << "make_unsteady" << std::endl;
      actions_before_adapt();

      this->Is_steady = false;

      actions_after_adapt();
    }

    // Set the bond number represented by ReInvFr
    void set_bond_number(const double& bond_number)
    {
      this->Bo = bond_number;

      this->ReInvFr_data_pt->set_value(0, this->Bo / this->Ca);
    }

    // Set the bond number represented by ReInvFr
    double* reynolds_number_inverse_froude_number_pt()
    {
      return this->ReInvFr_data_pt->value_pt(0);
    }

    // Set the Capillary without changing the Bond number
    void set_capillary_number(const double& capillary_number)
    {
      // Set the new capillary number
      this->Ca = capillary_number;

      this->ReInvFr_data_pt->set_value(0, this->Bo / this->Ca);
    }

    // Set the Reynolds number without changing the Bond number or Strouhal
    // number
    void set_reynolds_number(const double& reynolds_number)
    {
      // Set new Reynolds number
      this->Re = reynolds_number;

      // this->ReInvFr = invFr * this->Re;

      this->ReSt = this->St * this->Re;
    }

    // Set the Strouhal number
    void set_strouhal_number(const double& strouhal_number)
    {
      this->St = strouhal_number;

      this->ReSt = this->St * this->Re;
    }

    void set_viscosity_ratio(const double& viscosity_ratio)
    {
      this->Viscosity_ratio = viscosity_ratio;
    }

    // Set the contact angle, not this is overridden by the contact angle
    // function.
    void set_contact_angle(const double& contact_angle)
    {
      Contact_angle = contact_angle;
    }

    // Get the contact angle pointer. Note: accessing this pointer after the
    // problem has been deleted with cause undefined behaviour.
    double* get_contact_angle_pt()
    {
      return &Contact_angle;
    }

    // Set the maximum number of adapts when timestepping or solving the
    // static problem
    void set_max_adapt(const unsigned& max_adapt)
    {
      Max_adapt = max_adapt;
    }

    void set_directory(const std::string& dir_name)
    {
      this->doc_info().set_directory(dir_name);
    }

    std::string get_directory()
    {
      return this->doc_info().directory();
    }

    void set_doc_number(const unsigned& number)
    {
      this->doc_info().number() = number;
    }

    unsigned get_doc_number()
    {
      return this->doc_info().number();
    }

    Mesh* bulk_mesh_pt()
    {
      return Bulk_mesh_pt;
    }

    Mesh* free_surface_mesh_pt()
    {
      return Free_surface_mesh_pt;
    }

    Mesh* slip_surface_mesh_pt()
    {
      return Slip_boundary_mesh_pt;
    }

    // Function to print and debug the Jacobian and residual.
    // Compare specified Jacobian to a numerical one
    void debug_residuals()
    {
      DoubleVector residuals;

      get_residuals(residuals);

      oomph_info << "Residuals" << std::endl;
      for (unsigned i = ndof() - 10; i < ndof(); i++)
      {
        oomph_info << "i: " << i;
        printf(", res: %10.5g", residuals[i]);
        oomph_info << std::endl;
      }
    }

    void debug_jacobian()
    {
      oomph_info << "debug_jacobian" << std::endl;
      if (!this->Is_steady)
      {
        oomph_info << "WARNING: The problem is not steady! Expect differences "
                      "in the jacobian."
                   << std::endl;
      }

      DoubleVector residuals;
      DenseDoubleMatrix jacobian;
      DoubleVector residualsFD;
      DenseDoubleMatrix jacobianFD(ndof());

      get_jacobian(residuals, jacobian);
      jacobian.sparse_indexed_output(
        this->doc_info().directory() + "/jacJ.dat", 16, true);
      get_fd_jacobian(residualsFD, jacobianFD);
      jacobianFD.sparse_indexed_output(
        this->doc_info().directory() + "/jacfdJ.dat", 16, true);

      bool jacobians_are_different = compare_matrices(jacobian, jacobianFD);

      if (!jacobians_are_different)
      {
        oomph_info << "Computed Jacobian matches finite differenced Jacobian"
                   << std::endl;
      }
      else
      {
        oomph_info << "WARNING: Computed Jacobian is different to the finite "
                      "differenced Jacobian"
                   << std::endl;

        std::ofstream output_stream;
        output_stream.open(this->doc_info().directory() + "/dofs.txt");
        this->describe_dofs(output_stream);
        output_stream.close();
      }
    }

    void debug_mass_matrix()
    {
      oomph_info << "debug_mass_matrix" << std::endl;

      DoubleVector residuals;
      CRDoubleMatrix jacobian;
      CRDoubleMatrix mass_matrix;
      get_eigenproblem_matrices(mass_matrix, jacobian);

      std::ofstream output_stream;
      output_stream.open(this->doc_info().directory() + "/massJM.dat");
      mass_matrix.sparse_indexed_output(output_stream, 16, true);
      output_stream.close();
      output_stream.open(this->doc_info().directory() + "/jacJM.dat");
      jacobian.sparse_indexed_output(output_stream, 16, true);
      output_stream.close();

      DoubleVector residualsFD;
      DenseDoubleMatrix jacobianFD(ndof());
      DenseDoubleMatrix mass_matrixFD(ndof());
      get_fd_mass_matrix(residualsFD, jacobianFD, mass_matrixFD);

      output_stream.open(this->doc_info().directory() + "/massfdJM.dat");
      mass_matrixFD.sparse_indexed_output(output_stream, 16, true);
      output_stream.close();
      output_stream.open(this->doc_info().directory() + "/jacfdJM.dat");
      jacobianFD.sparse_indexed_output(output_stream, 16, true);
      output_stream.close();

      bool jacobians_are_different = compare_matrices(jacobian, jacobianFD);
      if (!jacobians_are_different)
      {
        oomph_info << "Computed Jacabian matches finite differenced Jacobian"
                   << std::endl;
      }
      else
      {
        oomph_info
          << "Computed Jacobian does not match finite differenced Jacobian"
          << std::endl;
      }
      bool mass_matrices_are_different =
        compare_matrices(mass_matrix, mass_matrixFD);
      if (!mass_matrices_are_different)
      {
        oomph_info
          << "Computed Mass matrix matches finite differenced Mass matrix"
          << std::endl;
      }
      else
      {
        oomph_info << "Computed Mass matrix does not match finite differenced "
                      "Mass matrix"
                   << std::endl;
      }

      // Output dofs if different
      if (mass_matrices_are_different || jacobians_are_different)
      {
        std::ofstream output_stream;
        output_stream.open(this->doc_info().directory() + "/dofs.txt");
        this->describe_dofs(output_stream);
        output_stream.close();
      }
    }

    //================================================================
    /// Get the full Mass by finite differencing
    //================================================================
    void get_fd_mass_matrix(DoubleVector& residuals,
                            DenseMatrix<double>& jacobian,
                            DenseMatrix<double>& mass_matrix)
    {
#ifdef OOMPH_HAS_MPI

      if (Problem_has_been_distributed)
      {
        OomphLibWarning("This is unlikely to work with a distributed problem",
                        OOMPH_EXCEPTION_LOCATION);
      }
#endif

      time_stepper_pt()->make_steady();

      DoubleVector tempRes;
      get_fd_jacobian(residuals, jacobian);

      time_stepper_pt()->undo_make_steady();

      const double dt = time_stepper_pt()->time_pt()->dt();
      DoubleVector residuals_unsteady;
      DenseMatrix<double> jacobian_unsteady;
      get_fd_jacobian(residuals_unsteady, jacobian_unsteady);

      mass_matrix = jacobian;
      for (unsigned i = 0; i < mass_matrix.nrow(); i++)
      {
        for (unsigned j = 0; j < mass_matrix.ncol(); j++)
        {
          /// The 2/3 is due to the BDF<2> scheme.
          mass_matrix(j, i) =
            (2.0 / 3.0) * dt * (-jacobian_unsteady(j, i) + jacobian(j, i));
        }
      }
    }

    // Open the trace files
    // Do this before documenting solution. Opening once per simulation rather
    // than at each time step. Closed in the destructor. This is public so
    // that they can be opened after adjusting problem parameters. If this is
    // the first time, then add the variable headers.
    void open_trace_files(bool first_time)
    {
      std::string filename;
      if (first_time)
      {
        // Open trace file
        filename = this->doc_info().directory() + "/trace.dat";
        Trace_file.open(filename);
        Trace_file.precision(16);
        add_header_to_trace_file();

        // Open trace file
        filename = this->doc_info().directory() + "/volume_trace.dat";

        Volume_trace_file.open(filename);
        Volume_trace_file << "doc_number ";
        Volume_trace_file << "prescribed_volume ";
        Volume_trace_file << "lagrange_multiplier" << std::endl;

        // Open trace file
        filename = this->doc_info().directory() + "/flux_trace.dat";
        Flux_trace_file.open(filename);
        Flux_trace_file << "doc_number ";
        Flux_trace_file << "time ";
        Flux_trace_file << "prescribed_flux ";
        Flux_trace_file << "lagrange_multiplier" << std::endl;

        // Open trace file
        filename = this->doc_info().directory() + "/contact_angle_trace.dat";
        Contact_angle_trace_file.open(filename);
        Contact_angle_trace_file << "doc_number ";
        Contact_angle_trace_file << "x ";
        Contact_angle_trace_file << "y ";
        Contact_angle_trace_file << "prescribed_angle ";
        Contact_angle_trace_file << "computed_angle ";
        Contact_angle_trace_file << "lagrange_multiplier";
        Contact_angle_trace_file << std::endl;

        filename = this->doc_info().directory() + "/inner_angle_trace.dat";
        Inner_angle_trace_file.open(filename);
        Inner_angle_trace_file << "doc_number ";
        Inner_angle_trace_file << "x ";
        Inner_angle_trace_file << "y ";
        Inner_angle_trace_file << "prescribed_angle ";
        Inner_angle_trace_file << "computed_angle ";
        Inner_angle_trace_file << "lagrange_multiplier";
        Inner_angle_trace_file << std::endl;
      }
      else
      {
        filename = this->doc_info().directory() + "/trace.dat";
        Trace_file.open(filename, std::ios_base::app);
        Trace_file.precision(16);

        filename = this->doc_info().directory() + "/volume_trace.dat";
        Volume_trace_file.open(filename, std::ios_base::app);

        Volume_trace_file << "doc_number ";
        filename = this->doc_info().directory() + "/flux_trace.dat";
        Flux_trace_file.open(filename, std::ios_base::app);

        filename = this->doc_info().directory() + "/contact_angle_trace.dat";
        Contact_angle_trace_file.open(filename, std::ios_base::app);

        filename = this->doc_info().directory() + "/inner_angle_trace.dat";
        Inner_angle_trace_file.open(filename, std::ios_base::app);
      }
    }

    void close_trace_files()
    {
      // Open trace file
      Trace_file.close();
      Volume_trace_file.close();
      Flux_trace_file.close();
      Contact_angle_trace_file.close();
      Inner_angle_trace_file.close();
    }

    // Add a header to the dat trace file
    void add_header_to_trace_file()
    {
      Trace_file << "doc_number ";
      Trace_file << "time ";
      Trace_file << "contact_angle ";
      Trace_file << "Bo ";
      Trace_file << "Ca ";
      Trace_file << "Re ";
      Trace_file << "St ";
      Trace_file << "ReSt ";
      Trace_file << "ReInvFr ";
      Trace_file << "wall_velocity ";
      Trace_file << "ExtP ";
      Trace_file << "h0 ";
      Trace_file << "h1 ";
      Trace_file << "H ";
      Trace_file << "ndof ";
      Trace_file << "min_element_length ";
      Trace_file << "inner_min_element_length ";
      Trace_file << "velocity_norm ";
      Trace_file << "deformation_norm ";
      Trace_file << "max_err ";
      Trace_file << "min_err ";
      Trace_file << "node_position_error_temporal_error ";
      Trace_file << "velocity_value_error_temporal_error ";
      Trace_file << "newton_iterations ";
      Trace_file << "sec_since_start";
      Trace_file << std::endl;
    }

    // Document the solution
    // Outputs the solution in the bulk and on the surfaces.
    // Uses the Z2 error estimator to compute an approximation to the error on
    // each element
    void doc_solution()
    {
      oomph_info << "doc_solution" << std::endl;
      int local_rank = 0;
#ifdef OOMPH_HAS_MPI
      OomphCommunicator* comm_pt = MPI_Helpers::communicator_pt();

      // Ensure each rank has completed previous computations
      MPI_Barrier(comm_pt->mpi_comm());
      local_rank = comm_pt->my_rank();
#endif

      // Set up error estimates within each bulk element
      ErrorEstimator* err_est_pt = Bulk_mesh_pt->spatial_error_estimator_pt();
      Bulk_mesh_pt->spatial_error_estimator_pt() = Z2_error_estimator_pt;
      double max_err = 0;
      double min_err = 0;
      compute_error_estimate(max_err, min_err);
      Bulk_mesh_pt->spatial_error_estimator_pt() = err_est_pt;

      if (local_rank == 0)
      {
        oomph_info << "Doc Number: " << this->doc_info().number() << std::endl;

        // Output stream
        std::ofstream output_stream;
        std::string filename;

        // Number of plot points
        unsigned npts = Plot_Parameters::Bulk_element_number_of_plot_points;
        unsigned npts_surface =
          Plot_Parameters::Surface_element_number_of_plot_points;

        // Output bulk domain
        if (Doc_Parameters::Doc_bulk)
        {
          filename = this->doc_info().directory() + "/soln" +
                     to_string(this->doc_info().number()) + ".dat";
          output_stream.open(filename);
          Bulk_mesh_pt->output(output_stream, npts);
          output_stream.close();
        }

        if (Doc_Parameters::Doc_free_surface && Free_surface_mesh_pt)
        {
          // Output free surface
          filename = this->doc_info().directory() + "/free_surface" +
                     to_string(this->doc_info().number()) + ".dat";
          output_stream.open(filename);
          output_stream << "x ";
          output_stream << "y ";
          output_stream << "u ";
          output_stream << "v ";
          output_stream << "w ";
          output_stream << "p ";
          output_stream << "lagrange_multiplier ";
          output_stream << std::endl;
          Free_surface_mesh_pt->output(output_stream, npts);
          output_stream.close();
        }

        // If the problem is static ...
        // ... output the volume ...
        if (Doc_Parameters::Doc_volume && Volume_constraint_mesh_pt)
        {
          Volume_trace_file << this->doc_info().number() << " ";
          dynamic_cast<VOLUME_CONSTRAINT_ELEMENT*>(
            Volume_constraint_mesh_pt->element_pt(0))
            ->output(Volume_trace_file);
        }
        // ... output the flux on the upper boundary
        if (Doc_Parameters::Doc_flux)
        {
          if (Flux_mesh_pt)
          {
            filename = this->doc_info().directory() + "/flux_surface" +
                       to_string(this->doc_info().number()) + ".dat";
            output_stream.open(filename);
            doc_flux(output_stream, npts_surface);
            output_stream.close();
          }

          if (Net_flux_mesh_pt)
          {
            Flux_trace_file << this->doc_info().number() << " ";
            dynamic_cast<NET_FLUX_ELEMENT*>(Net_flux_mesh_pt->element_pt(0))
              ->output(Flux_trace_file);
          }
        }

        // Output the slip on the slip boundary

        if (Doc_Parameters::Doc_slip && Slip_boundary_mesh_pt)
        {
          filename = this->doc_info().directory() + "/slip_surface" +
                     to_string(this->doc_info().number()) + ".dat";
          output_stream.open(filename);
          doc_slip(output_stream, npts_surface);
          output_stream.close();
        }

        // Output the no-penetration
        if (Doc_Parameters::Doc_no_penetration &&
            No_penetration_boundary_mesh_pt)
        {
          filename = this->doc_info().directory() + "/no_penetration_surface" +
                     to_string(this->doc_info().number()) + ".dat";
          output_stream.open(filename);
          output_stream << "x ";
          output_stream << "y ";
          output_stream << "u ";
          output_stream << "v ";
          output_stream << "w ";
          output_stream << "lagrange_multiplier ";
          output_stream << std::endl;
          No_penetration_boundary_mesh_pt->output(output_stream);
          output_stream.close();
        }

        if (Doc_Parameters::Doc_contact_angle && Contact_angle_mesh_pt)
        {
          // Output the contact angle
          Contact_angle_trace_file << this->doc_info().number() << " ";
          doc_contact_angle(Contact_angle_trace_file);

          Inner_angle_trace_file << this->doc_info().number() << " ";
          doc_inner_angle(Inner_angle_trace_file);
        }

        if (Doc_Parameters::Doc_trace)
        {
          // Output trace which includes all parameters and global data
          doc_trace(max_err, min_err);
        }

        // Document height drop to std::cout
        dynamic_cast<HEIGHT_ELEMENT*>(Height_mesh_pt->element_pt(0))
          ->output(std::cout);

        // Bump up counter
        this->doc_info().number()++;
      }

#ifdef OOMPH_HAS_MPI
      // Ensure each rank doesn't start further computations, until
      // documenting has completed.
      MPI_Barrier(comm_pt->mpi_comm());
#endif

    } // end_of_doc_solution

    // Doc adaptivity targets
    void doc_adaptivity_targets(std::ostream& outfile)
    {
      Bulk_mesh_pt->doc_adaptivity_targets(outfile);
    }

    // Document the flux on the upper boundary
    void doc_flux(std::ostream& out, const unsigned& n_plot_points)
    {
      // Output dat file header
      out << "x y u v w" << std::endl;

      Flux_mesh_pt->output(out, n_plot_points);
    }

    // Document the slip on the slip boundary
    void doc_slip(std::ostream& out, const unsigned& n_plot_points)
    {
      // Output dat file header
      out << "x y l_x l_y l_z n_x n_y u v w" << std::endl;

      Slip_boundary_mesh_pt->output(out, n_plot_points);
    }


    // Document the slip on the slip boundary
    void doc_contact_angle(std::ostream& out)
    {
      // Output contact angle
      dynamic_cast<CONTACT_LINE_ELEMENT*>(Contact_angle_mesh_pt->element_pt(0))
        ->output(out);
    }

    // Document the slip on the slip boundary
    void doc_inner_angle(std::ostream& out)
    {
      // Output contact angle
      dynamic_cast<CONTACT_LINE_ELEMENT*>(Contact_angle_mesh_pt->element_pt(1))
        ->output(out);
    }

    // Doc the trace
    void doc_trace(const double& max_err, const double& min_err)
    {
      Node* left_hand_node_pt = 0;
      find_interface_end_node(false, left_hand_node_pt);
      Node* right_hand_node_pt = 0;
      find_interface_end_node(true, right_hand_node_pt);


      // Document the contact angle (in degrees),
      Trace_file << this->doc_info().number() << " ";
      Trace_file << this->time_pt()->time() << " ";
      Trace_file << Contact_angle * 180.0 / MathematicalConstants::Pi << " ";
      // the parameters,
      Trace_file << Bo << " ";
      Trace_file << Ca << " ";
      Trace_file << Re << " ";
      Trace_file << St << " ";
      Trace_file << ReSt << " ";
      Trace_file << ReInvFr_data_pt->value(0) << " ";
      Trace_file << Wall_velocity_data_pt->value(0) << " ";
      // the external pressure,
      Trace_file << External_pressure_data_pt->value(0) << " ";
      // the height of the interface at the centre of the container,
      Trace_file << left_hand_node_pt->x(1) << " ";
      // the height at the wall,
      Trace_file << right_hand_node_pt->x(1) << " ";
      // the height difference,
      Trace_file << right_hand_node_pt->x(1) - left_hand_node_pt->x(1) << " ";
      // the number of degrees of freedom
      Trace_file << this->ndof() << " ";
      // the element length of the small elements about the contact angle
      Trace_file << Mesh_Control_Parameters::min_element_length << " ";
      // the element length of the small elements about the contact angle
      Trace_file << Mesh_Control_Parameters::inner_min_element_length << " ";
      // the global velocity norm
      Trace_file << this->global_velocity_norm() << " ";
      // the global deformation norm
      Trace_file << this->global_deformation_norm() << " ";
      // the estimate of the max error,
      Trace_file << max_err << " ";
      // the estimate of the min error.
      Trace_file << min_err << " ";
      // and the estimate of the temporal error.
      Trace_file << node_position_error_temporal_error() << " ";
      Trace_file << velocity_value_error_temporal_error() << " ";
      Trace_file << this->Nnewton_iter_taken << " ";

      // std::chrono::high_resolution_clock::time_point current_time =
      //  std::chrono::high_resolution_clock::now();
      // std::chrono::duration<double> diff_time =
      //  std::chrono::duration_cast<std::chrono::seconds>(current_time -
      //                                                   this->Start_time);
      // Trace_file << diff_time.count();

      // End of line
      Trace_file << std::endl;
    }

    virtual void read(std::ifstream& restart_file, bool& unsteady_restart)
    {
      this->flush_sub_meshes();
      this->flush_global_data();

      this->add_global_data(External_pressure_data_pt);
      this->add_sub_mesh(Net_flux_mesh_pt);
      this->add_sub_mesh(Bulk_mesh_pt);
      this->add_sub_mesh(Free_surface_mesh_pt);
      this->add_sub_mesh(Slip_boundary_mesh_pt);
      this->add_sub_mesh(No_penetration_boundary_mesh_pt);
      this->add_sub_mesh(Flux_mesh_pt);
      this->add_sub_mesh(Volume_computation_mesh_pt);
      this->add_sub_mesh(Contact_angle_mesh_pt);
      this->add_sub_mesh(Volume_constraint_mesh_pt);
      this->rebuild_global_mesh();
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;

      Problem::read(restart_file, unsteady_restart);

      this->add_sub_mesh(Height_mesh_pt);
      this->add_sub_mesh(Eigensolution_slip_mesh_pt);
      this->add_sub_mesh(Singularity_scaling_mesh_pt);
      this->add_sub_mesh(Pressure_contribution_mesh_1_pt);
      this->add_sub_mesh(Pressure_contribution_mesh_2_pt);
      this->add_global_data(ReInvFr_data_pt);
      this->add_global_data(Wall_velocity_data_pt);
      this->rebuild_global_mesh();
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    virtual void dump(std::ofstream& dump_file)
    {
      this->flush_sub_meshes();
      this->flush_global_data();

      this->add_global_data(External_pressure_data_pt);
      this->add_sub_mesh(Net_flux_mesh_pt);
      this->add_sub_mesh(Bulk_mesh_pt);
      this->add_sub_mesh(Free_surface_mesh_pt);
      this->add_sub_mesh(Slip_boundary_mesh_pt);
      this->add_sub_mesh(No_penetration_boundary_mesh_pt);
      this->add_sub_mesh(Flux_mesh_pt);
      this->add_sub_mesh(Volume_computation_mesh_pt);
      this->add_sub_mesh(Contact_angle_mesh_pt);
      this->add_sub_mesh(Volume_constraint_mesh_pt);
      this->rebuild_global_mesh();
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;

      Problem::dump(dump_file);

      this->add_sub_mesh(Height_mesh_pt);
      this->add_sub_mesh(Eigensolution_slip_mesh_pt);
      this->add_sub_mesh(Singularity_scaling_mesh_pt);
      this->add_sub_mesh(Pressure_contribution_mesh_1_pt);
      this->add_sub_mesh(Pressure_contribution_mesh_2_pt);
      this->add_global_data(ReInvFr_data_pt);
      this->add_global_data(Wall_velocity_data_pt);
      this->rebuild_global_mesh();
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    // Create a restart file
    // Does not include all the problem and system parameters
    void create_restart_file()
    {
      // Save current solution
      std::ofstream dump_filestream;
      std::string restart_filename = this->doc_info().directory() + "/restart" +
                                to_string(this->doc_info().number()) + ".dat";
      dump_filestream.open(restart_filename);
      dump_filestream.precision(16);
      // actions_before_adapt();
      dump(dump_filestream);
      // actions_after_adapt();
      dump_filestream.close();

      std::string triangulateio_filename =
        this->doc_info().directory() + "/triangulateio" +
        to_string(this->doc_info().number()) + ".dat";
      dump_filestream.open(triangulateio_filename);
      Bulk_mesh_pt->dump_triangulateio(dump_filestream);
      dump_filestream.close();
    }

    double get_centre_point_z()
    {
      return Inner_corner_solid_node_pt->x(1);
    }

  private:
    // Set the elements internal variable and function pointers
    void setup_bulk_elements()
    {
      // Loop over the elements to set the consitutive law and jacobian
      unsigned n_bulk = Bulk_mesh_pt->nelement();
      for (unsigned e = 0; e < n_bulk; e++)
      {
        // Upcast from GeneralisedElement to the present element
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

        // Set the Reynolds number
        el_pt->re_pt() = &Re;

        // Set the Womersley number
        el_pt->re_st_pt() = &ReSt;

        // Set viscosity ratio
        el_pt->viscosity_ratio_pt() = &Viscosity_ratio;

        // Set the product of the Reynolds number and the inverse of the
        // Froude number
        el_pt->re_invfr_pt() = ReInvFr_data_pt->value_pt(0);

        // Set the direction of gravity
        el_pt->g_pt() = &Global_Physical_Parameters::G;

        // Set the constitutive law
        el_pt->constitutive_law_pt() = Constitutive_law_pt;
        el_pt->disable_inertia();
        // el_pt->enable_evaluate_jacobian_by_fd();

        el_pt->time_stepper_pt() = this->time_stepper_pt();

        el_pt->add_external_data(ReInvFr_data_pt);

        // Augmented elements close to the corner
        // Check distance from
        // s centre is centre of mass of a uniform triangle, so (1/3,1/3) for
        // Triangle[(1,0),(0,1),(0,0)]
        Vector<double> s_centre(2, 1.0 / 3.0);
        Vector<double> element_centre_x(2, 0.0);
        el_pt->get_x(s_centre, element_centre_x);
        double dist = 0;
        for (unsigned i = 0; i < 2; i++)
        {
          dist += pow(
            element_centre_x[i] - Contact_line_solid_node_pt->position(i), 2.0);
        }
        dist = pow(dist, 0.5);

        // If the distance to the corner is within the "inner" region, ...
        const double inner_radius = 0.3;
        if (dist < inner_radius)
        {
          el_pt->augment();
          Augmented_bulk_element_number.push_back(e);
        }
      }
      oomph_info << Augmented_bulk_element_number.size()
                 << " augmented elements" << std::endl;
    }

    // Find interface end node
    void find_interface_end_node(const bool& node_on_outer_boundary,
                                 Node*& node_pt)
    {
      unsigned boundary_id;
      if (node_on_outer_boundary)
      {
        boundary_id = Outer_boundary_with_slip_id;
      }
      else
      {
        boundary_id = Inner_boundary_id;
      }

      unsigned n_nod = Bulk_mesh_pt->nboundary_node(Free_surface_boundary_id);
      for (unsigned inod = 0; inod < n_nod; inod++)
      {
        // Get boundary node
        Node* nod_pt =
          Bulk_mesh_pt->boundary_node_pt(Free_surface_boundary_id, inod);

        if (nod_pt->is_on_boundary(boundary_id))
        {
          node_pt = nod_pt;

          return;
        }
      }
    }

    // Find corner node and return whether it has been found
    bool find_corner_node(const unsigned& first_boundary_id,
                          const unsigned& second_boundary_id,
                          SolidNode*& node_pt)
    {
      unsigned n_nod = Bulk_mesh_pt->nboundary_node(first_boundary_id);
      for (unsigned inod = 0; inod < n_nod; inod++)
      {
        // Get boundary node
        Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(first_boundary_id, inod);

        if (nod_pt->is_on_boundary(second_boundary_id))
        {
          node_pt = dynamic_cast<SolidNode*>(nod_pt);

          return true;
        }
      }
      return false;
    }

  public:
    void get_z2_error(double& max_err, double& min_err)
    {
      ErrorEstimator* err_est_pt = Bulk_mesh_pt->spatial_error_estimator_pt();
      Bulk_mesh_pt->spatial_error_estimator_pt() = Z2_error_estimator_pt;
      compute_error_estimate(max_err, min_err);
      Bulk_mesh_pt->spatial_error_estimator_pt() = err_est_pt;
    }

  private:
    // Compute error estimates and save to elements
    void compute_error_estimate(double& max_err, double& min_err)
    {
      // Get error estimator
      ErrorEstimator* err_est_pt = Bulk_mesh_pt->spatial_error_estimator_pt();

      // Get/output error estimates
      unsigned n_elements = Bulk_mesh_pt->nelement();
      Vector<double> elemental_error(n_elements);

      // We need a dynamic cast, get_element_errors from the Bulk_mesh_pt
      // Dynamic cast is used because get_element_errors require a Mesh* ans
      // not a SolidMesh*
      Mesh* fluid_mesh_pt = dynamic_cast<Mesh*>(Bulk_mesh_pt);
      set_contact_line_node_pt();
      err_est_pt->get_element_errors(fluid_mesh_pt, elemental_error);

      // Set errors for post-processing and find extrema
      max_err = 0.0;
      min_err = 1e6;
      for (unsigned e = 0; e < n_elements; e++)
      {
        dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))
          ->set_error(elemental_error[e]);

        max_err = std::max(max_err, elemental_error[e]);
        min_err = std::min(min_err, elemental_error[e]);
      }

      oomph_info << "Max error is " << max_err << std::endl;
      oomph_info << "Min error is " << min_err << std::endl;
    }

    // Set the boundary conditions
    void set_boundary_conditions()
    {
      // pin_total_velocity_eqn();

      pin_velocity_on_boundary(w_index, Free_surface_boundary_id);

      pin_velocity_on_boundary(w_index, Outer_boundary_with_slip_id);
      pin_position_on_boundary(0, Outer_boundary_with_slip_id);
      pin_no_penetration_lagrange_multiplier(0.0e-3);
      pin_velocity_on_boundary(u_index, Outer_boundary_with_slip_id);

      pin_velocity_on_boundary(u_index, Upper_boundary_id);
      set_velocity_on_boundary(u_index, Upper_boundary_id, 0.0);
      pin_velocity_on_boundary(w_index, Upper_boundary_id);
      pin_position_on_boundary(1, Upper_boundary_id);

      pin_velocity_on_boundary(u_index, Inner_boundary_id);
      set_velocity_on_boundary(u_index, Inner_boundary_id, 0.0);
      pin_velocity_on_boundary(w_index, Inner_boundary_id);
      pin_position_on_boundary(0, Inner_boundary_id);

      // pin_horizontal_displacement();
      pin_azimuthal_velocity();

      if (!this->Is_steady)
      {
        // Fix the volume constraint. This fixes the external pressure so we
        // don't have to pin an internal pressure.
        pin_volume_constraint();
        unpin_interior_pressure();
      }
      else
      {
        // Pin the flux constraint
        if (Net_flux_mesh_pt)
        {
          pin_flux_constraint();
        }

        // Pin an interior pressure
        pin_interior_pressure();
      }

      if (Slip_Parameters::slip_length == 0 &&
          std::abs(Wall_velocity_data_pt->value(0)) >= 1e-8)
      {
        pin_velocity_on_boundary(v_index, Outer_boundary_with_slip_id);
        pin_contact_line();
      }
    }

    void set_dirichlet_bc_on_boundary(const unsigned b)
    {
      // local boolean
      const bool pin_bc = false;

      // Loop over the nodes on the boundary
      unsigned n_boundary_element = Bulk_mesh_pt->nboundary_element(b);
      for (unsigned n = 0; n < n_boundary_element; n++)
      {
        // Upcast from GeneralsedElement to the present element
        ELEMENT* el_pt =
          dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b, n));

        const unsigned n_node = el_pt->nnode();
        for (unsigned i = 0; i < n_node; i++)
        {
          Node* node_pt = el_pt->node_pt(i);
          if (node_pt->is_on_boundary(b))
          {
            Vector<double> x(2);
            const unsigned DIM = 2;
            for (unsigned d = 0; d < DIM; d++)
            {
              x[d] = node_pt->x(d);
            }

            // Vector<double> u = parameters::velocity_singular_fct(x);
            Vector<double> u = parameters::boundary_condition_velocity_fct(x);

            if (pin_bc)
            {
              // oomph_info << u[0] << ", " << u[1] << std::endl;

              // el_pt->pin_momentum_eqn(i, 0);
              // el_pt->pin_momentum_eqn(i, 1);
              el_pt->set_velocity_value(i, 0, u[0]);
              el_pt->set_velocity_value(i, 1, u[1]);
            }
            else
            {
              // Its imposed value is 1.0
              el_pt->set_velocity_dirichlet_value_on_node(i, 0, u[0]);
              el_pt->set_velocity_dirichlet_value_on_node(i, 1, u[1]);

              // The x velocity component at node j is subject to Dirichlet BC
              // Its imposed value is 0.0 which is already the value by default
              // in the vector Imposed_velocity_values_at_node. So there is no
              // need to call the function set_velocity_dirichlet_value_on_node
              el_pt->impose_velocity_dirichlet_bc_on_node(i, 0);

              // The y velocity component at node j is subject to Dirichlet BC
              // Its imposed value is 0.0 which is already the value by default
              // in the vector Imposed_velocity_values_at_node. So there is no
              // need to call the function set_velocity_dirichlet_value_on_node
              el_pt->impose_velocity_dirichlet_bc_on_node(i, 1);
            }
          }
        }
      }
    }

  public:
    void perturb_free_surface()
    {
      const unsigned n_node =
        Bulk_mesh_pt->nboundary_node(Free_surface_boundary_id);
      for (unsigned n = 0; n < n_node; n++)
      {
        double r =
          Bulk_mesh_pt->boundary_node_pt(Free_surface_boundary_id, n)->x(0);
        const double amplitude = 0.01;
        Bulk_mesh_pt->boundary_node_pt(Free_surface_boundary_id, n)->x(1) +=
          amplitude * (-7.0 / 60.0 + pow(r, 2.0) / 2.0 - pow(r, 3.0) / 3.0);
      }
    }

    void perturb_free_surface2()
    {
      const unsigned n_node =
        Bulk_mesh_pt->nboundary_node(Free_surface_boundary_id);
      for (unsigned n = 0; n < n_node; n++)
      {
        double r =
          Bulk_mesh_pt->boundary_node_pt(Free_surface_boundary_id, n)->x(0);
        if (r > 0.5)
        {
          const double amplitude = 0.01;
          Bulk_mesh_pt->boundary_node_pt(Free_surface_boundary_id, n)->x(1) +=
            amplitude * (r - 0.5) * (r - 0.5) * (1.0 - r);
        }
      }
    }

    void perturb_displacement(const double& amplitude)
    {
      const unsigned n_node = Bulk_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        double r = Bulk_mesh_pt->node_pt(n)->x(0);
        double z = Bulk_mesh_pt->node_pt(n)->x(1);
        // Constrains for the perturbation
        // Volume is conserved
        // y'(0) = 0
        // y'(1) = 0
        // Bulk_mesh_pt->node_pt(n)->x(1) +=
        //  amplitude * (-7.0 / 60.0 + pow(r, 2.0) / 2.0 - pow(r, 3.0) / 3.0)
        //  * (3.5 - z);
        Bulk_mesh_pt->node_pt(n)->x(1) +=
          amplitude * ((tanh(7.0 * (r - 0.4)) - 0.646608127717086) /
                       1.992181886660662 * (3.5 - z));
      }
    }

    void perturb_fluid()
    {
      const unsigned n_node = Bulk_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        double r = Bulk_mesh_pt->node_pt(n)->x(0) - 0.5;
        double z = Bulk_mesh_pt->node_pt(n)->x(1) - 3.5 / 2.0;
        Bulk_mesh_pt->node_pt(n)->set_value(u_index, z);
        Bulk_mesh_pt->node_pt(n)->set_value(v_index, -r);
      }
    }

  private:
    void unpin_all_boundaries()
    {
      unsigned n_bound = Bulk_mesh_pt->nboundary();
      for (unsigned b = 0; b < n_bound; b++)
      {
        // Loop over the nodes on the boundary
        unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
        for (unsigned n = 0; n < n_boundary_node; n++)
        {
          Bulk_mesh_pt->boundary_node_pt(b, n)->unpin(u_index);
          Bulk_mesh_pt->boundary_node_pt(b, n)->unpin(v_index);
          Bulk_mesh_pt->boundary_node_pt(b, n)->unpin(w_index);
          dynamic_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(b, n))
            ->unpin_position(0);
          dynamic_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(b, n))
            ->unpin_position(1);
        }
      }
    }

  public:
    void pin_volume_constraint()
    {
      External_pressure_data_pt->pin(0);
    }

    void pin_flux_constraint()
    {
      if (Net_flux_mesh_pt)
      {
        dynamic_cast<NET_FLUX_ELEMENT*>(Net_flux_mesh_pt->element_pt(0))
          ->internal_data_pt(0)
          ->pin(0);
      }
    }

  private:
  public:
    void pin_velocity_on_boundary(const unsigned& j, const unsigned& b)
    {
      // Loop over the nodes on the boundary
      unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        // No conditons on the free surface. The dynamic and kinematic
        // conditions are handled by the interface elements
        Bulk_mesh_pt->boundary_node_pt(b, n)->pin(j);
      }
    }

    void set_velocity_on_boundary(const unsigned& j,
                                  const unsigned& b,
                                  const double& value)
    {
      // Loop over the nodes on the boundary
      unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        // No conditons on the free surface. The dynamic and kinematic
        // conditions are handled by the interface elements
        Bulk_mesh_pt->boundary_node_pt(b, n)->set_value(j, value);
      }
    }

  private:
    void set_velocity_on_upper_boundary_to_parabola()
    {
      unsigned n_boundary_node =
        this->Bulk_mesh_pt->nboundary_node(Upper_boundary_id);
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        Node* node_pt =
          this->Bulk_mesh_pt->boundary_node_pt(Upper_boundary_id, n);
        double x = node_pt->x(0);
        double flux = 0;
        Flux_Parameters::flux_fct(this->time(), flux);
        double U = 2 * (Wall_velocity_data_pt->value(0) + flux) * x * x -
                   Wall_velocity_data_pt->value(0) - 2 * flux;
        node_pt->set_value(v_index, U);
      }
    }

    void pin_total_velocity_eqn()
    {
      // Loop over the elements to set the consitutive law and jacobian
      unsigned n_bulk = Bulk_mesh_pt->nelement();
      for (unsigned e = 0; e < n_bulk; e++)
      {
        // Upcast from GeneralisedElement to the present element
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
        for (unsigned n = 0; n < 6; n++)
        {
          for (unsigned i = 0; i < 3; i++)
          {
            el_pt->pin_total_velocity_eqn(n, i);
          }
        }
      }
    }

    void pin_position_on_boundary(const unsigned& j, const unsigned& b)
    {
      // Loop over the nodes on the boundary
      unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        // No conditons on the free surface. The dynamic and kinematic
        // conditions are handled by the interface elements
        dynamic_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(b, n))
          ->pin_position(j);
      }
    }

  public:
    void pin_fluid()
    {
      oomph_info << "pin_fluid" << std::endl;
      unsigned n_element = Bulk_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(n))->pin();
      }

      pin_kinematic_lagrange_multiplier(0.0);

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void pin_solid_boundaries()
    {
      for (unsigned j = 0; j < 2; j++)
      {
        for (unsigned b = 0; b < 4; b++)
        {
          pin_position_on_boundary(j, b);
        }
      }
      pin_volume_constraint();

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    };

    void pin_fs_solid_boundary()
    {
      for (unsigned j = 0; j < 2; j++)
      {
        for (unsigned b = 0; b < 4; b++)
        {
          if (b == Free_surface_boundary_id)
          {
            pin_position_on_boundary(j, b);
          }
        }
      }
      pin_volume_constraint();

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    };

    void pin_solid_boundaries(const unsigned& b)
    {
      for (unsigned j = 0; j < 2; j++)
      {
        pin_position_on_boundary(j, b);
      }
      if (b == Free_surface_boundary_id)
      {
        pin_volume_constraint();
      }

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    };

    void pin_fluid_boundary(const unsigned& b)
    {
      for (unsigned j = 0; j < 3; j++)
      {
        pin_velocity_on_boundary(j, b);
      }
      if (b == Free_surface_boundary_id)
      {
        pin_kinematic_lagrange_multiplier(0.0);
      }

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    };

    void pin_solid_boundaries_except_upper()
    {
      for (unsigned j = 0; j < 2; j++)
      {
        for (unsigned b = 0; b < 4; b++)
        {
          if (b != Upper_boundary_id)
          {
            pin_position_on_boundary(j, b);
          }
        }
      }
      pin_volume_constraint();

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    };

  private:
    void pin_azimuthal_velocity()
    {
      unsigned n_node = Bulk_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        Bulk_mesh_pt->node_pt(n)->pin(w_index);
        Bulk_mesh_pt->node_pt(n)->set_value(w_index, 0.0);
      }
    }

  public:
    void pin_horizontal_displacement()
    {
      unsigned n_node = Bulk_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        dynamic_cast<SolidNode*>(Bulk_mesh_pt->node_pt(n))->pin_position(0);
      }

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void pin_vertical_displacement()
    {
      unsigned n_node = Bulk_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        dynamic_cast<SolidNode*>(Bulk_mesh_pt->node_pt(n))->pin_position(1);
      }
    }

    void pin_all_pressure()
    {
      unsigned n_element = Bulk_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(n));
        el_pt->pin_pressure();
      }
    }

    void pin_interior_pressure()
    {
      SolidNode* node_pt = 0;
      find_corner_node(Outer_boundary_with_slip_id, Upper_boundary_id, node_pt);
      node_pt->pin(3);
      node_pt->set_value(3, 0.0);
      // Inner_corner_solid_node_pt->pin(3);
      // Inner_corner_solid_node_pt->set_value(3, 0.0);
    }

  private:
    void unpin_interior_pressure()
    {
      SolidNode* node_pt = 0;
      find_corner_node(Outer_boundary_with_slip_id, Upper_boundary_id, node_pt);
      node_pt->unpin(3);
      // Inner_corner_solid_node_pt->unpin(3);
    }

  public:
    void pin_solid()
    {
      oomph_info << "pin_solid" << std::endl;
      // Pin all solid node positions
      unsigned n_node = Bulk_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        dynamic_cast<SolidNode*>(Bulk_mesh_pt->node_pt(n))->pin_position(0);
        dynamic_cast<SolidNode*>(Bulk_mesh_pt->node_pt(n))->pin_position(1);
      }

      pin_volume_constraint();
      pin_kinematic_lagrange_multiplier(0.0);
    }

  private:
    // Set the contact line node pointer
    void set_contact_line_node_pt()
    {
      find_corner_node(Outer_boundary_with_slip_id,
                       Free_surface_boundary_id,
                       Contact_line_solid_node_pt);
      find_corner_node(Inner_boundary_id,
                       Free_surface_boundary_id,
                       Inner_corner_solid_node_pt);
    }

    // Create both the Z2 and corner error estimators
    // This is called in the constructor, then we can swap between the two as
    // needed.
    void create_error_estimators()
    {
      oomph_info << "create_error_estimators" << std::endl;
      // Create the Z2 error estimator
      Z2_error_estimator_pt = new Z2ErrorEstimator;
      // And set a reference flux norm to avoid high errors for the static
      // solution.
      Z2_error_estimator_pt->reference_flux_norm() = 0.01;

      // Create the corner estimator to refine sufficently around the contact
      // line and inner corner.
      create_corner_error_estimator();
    }

    // Create the corner estimator to refine sufficently around the contact
    // line and inner corner.
    void create_corner_error_estimator()
    {
      // If one already exists
      if (Corner_error_estimator_pt)
      {
        // Delete it
        delete Corner_error_estimator_pt;
      }

      // We need to know which is the corner node
      set_contact_line_node_pt();
      // Create error estimator with the mesh control parameters.
      Corner_error_estimator_pt = new CornerErrorEstimator(
        Contact_line_solid_node_pt,
        Inner_corner_solid_node_pt,
        &Mesh_Control_Parameters::min_element_length,
        &Mesh_Control_Parameters::inner_min_element_length,
        Mesh_Control_Parameters::element_length_ratio);
    }

    bool is_almost_static()
    {
      double velocity_norm = 0;
      velocity_norm = global_velocity_norm();

      if (std::abs(Wall_velocity_data_pt->value(0)) < 1e-8 ||
          velocity_norm < 1e-8)
      {
        return true;
      }
      else
      {
        return false;
      }
    }

    // Set the error estimator
    void set_error_estimator()
    {
      oomph_info << "set_error_estimator()" << std::endl;

      double velocity_norm = 0;
      velocity_norm = global_velocity_norm();

      if (std::abs(Wall_velocity_data_pt->value(0)) < 1e-8 ||
          velocity_norm < 1e-8)
      {
        Using_contact_angle_error_estimator = true;
        create_corner_error_estimator();

        // ... use the corner error estimator
        Bulk_mesh_pt->spatial_error_estimator_pt() = Corner_error_estimator_pt;

        // Set the refinement tolerances
        Bulk_mesh_pt->min_permitted_error() =
          Mesh_Control_Parameters::Min_permitted_mesh_residual;
        Bulk_mesh_pt->max_permitted_error() =
          Mesh_Control_Parameters::Max_permitted_mesh_residual;
      }
      else
      {
        Using_contact_angle_error_estimator = false;
        // ... use the Z2 error estimator.
        Bulk_mesh_pt->spatial_error_estimator_pt() = Z2_error_estimator_pt;

        // Set the refinement tolerances
        Bulk_mesh_pt->min_permitted_error() =
          Mesh_Control_Parameters::Min_permitted_z2_error;
        Bulk_mesh_pt->max_permitted_error() =
          Mesh_Control_Parameters::Max_permitted_z2_error;
      }
      if (Using_contact_angle_error_estimator)
      {
        std::cout << "Using contact angle error estimator" << std::endl;
      }
      else
      {
        std::cout << "Using Z2 error estimator" << std::endl;
      }
    }

    // Set an approximation to the initial pressure gradient -- not used
    void set_initial_pressure()
    {
      const unsigned n_element = this->Bulk_mesh_pt->nelement();
      for (unsigned i_element = 0; i_element < n_element; i_element++)
      {
        // Upcast from GeneralisedElement to the present element
        ELEMENT* el_pt =
          dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i_element));

        // In these 9-node elements, the 4-th node is always positioned at
        // s1=0,s2=0 where s1, s2 are the local coordinates (and the "first"
        // node is the 0-th)
        const double eulerian_z_pos_middle_node = el_pt->node_pt(4)->x(1);

        // Determine the value of the pressure at this node
        const double p_val_at_middle_node = Global_Physical_Parameters::G[1] *
                                            this->ReInvFr_data_pt->value(0) *
                                            (eulerian_z_pos_middle_node - 3.5);

        // Specify the pressure analytically
        el_pt->fix_pressure(0, p_val_at_middle_node);
        el_pt->fix_pressure(2, 0.0);
        el_pt->fix_pressure(1, 0.0);

        el_pt->free_pressure(0);
        el_pt->free_pressure(1);
        el_pt->free_pressure(2);
      }
    }

    /// Update the problem specs before next timestep:
    /// Set Dirchlet boundary conditions from exact solution.
    void actions_before_implicit_timestep()
    {
      set_velocity_on_upper_boundary_to_parabola();
    }

    // Actions before adapt
    void actions_before_adapt()
    {
      oomph_info << "actions_before_adapt" << std::endl;

      unaugment_elements();

      // Reset boundary conditions
      unpin_interior_pressure();
      unpin_all_boundaries();

      // Reset error estimator for bulk mesh
      set_error_estimator();

      //======================================================================
      // Backup Lagrange multipliers
      //======================================================================
      this->backup_lagrange_multipliers();

      //======================================================================
      // Delete all non-refineable elements
      //======================================================================
      delete_non_refineable_elements();

      //======================================================================
      // Rebuild the global mesh
      //======================================================================
      this->rebuild_global_mesh();
    }

    void unaugment_elements()
    {
      const unsigned n_aug_bulk = Augmented_bulk_element_number.size();
      for (unsigned e = 0; e < n_aug_bulk; e++)
      {
        // Augment elements
        // Upcast from GeneralisedElement to the present element
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
          Bulk_mesh_pt->element_pt(Augmented_bulk_element_number[e]));
        el_pt->unaugment();
      }
      Augmented_bulk_element_number.clear();
    }

    void delete_non_refineable_elements()
    {
      delete_elements(No_penetration_boundary_mesh_pt);
      delete_elements(Flux_mesh_pt);
      delete_elements(Net_flux_mesh_pt);
      delete_elements(Slip_boundary_mesh_pt);
      delete_elements(Contact_angle_mesh_pt);
      delete_elements(Volume_constraint_mesh_pt);
      delete_elements(Volume_computation_mesh_pt);
      delete_elements(Free_surface_mesh_pt);
      if (Height_mesh_pt)
      {
        delete_elements(Height_mesh_pt);
      }

      delete_elements(Eigensolution_slip_mesh_pt);
      delete_elements(Singularity_scaling_mesh_pt);
      delete_elements(Pressure_contribution_mesh_1_pt);
      delete_elements(Pressure_contribution_mesh_2_pt);
      delete Pressure_contribution_geom_mesh_1_pt;
      delete Pressure_contribution_geom_mesh_2_pt;
    }

    // Actions before adapt
    void actions_after_adapt()
    {
      oomph_info << "actions_after_adapt" << std::endl;
      //======================================================================
      // Setup the remaining parts of the bulk mesh
      //======================================================================
      setup_refineable_elements();

      //======================================================================
      // Create non-refineable elements
      //======================================================================
      create_non_refineable_elements();

      //======================================================================
      // Restore the backed up mesh if there is one
      //======================================================================
      this->restore_lagrange_multipliers();

      this->rebuild_global_mesh();
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
      //======================================================================
      // Set the boundary conditions
      //======================================================================
      set_boundary_conditions();

      // Rebuild the global mesh
      this->rebuild_global_mesh();

      // Set the timestepper to steady/unsteady
      if (Is_steady)
      {
        time_stepper_pt()->make_steady();
      }
      else
      {
        time_stepper_pt()->undo_make_steady();
      }

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    // void actions_after_newton_solve()
    // {
    //   Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
    //   if (Free_surface_mesh_pt)
    //   {
    //     set_kinematic_lagrange_multiplier(0.0);
    //   }
    // }

    void setup_refineable_elements()
    {
      // Reset the Contact_line_solid_node_pt and Inner_corner_solid_node_pt
      find_corner_node(Outer_boundary_with_slip_id,
                       Free_surface_boundary_id,
                       Contact_line_solid_node_pt);
      find_corner_node(Inner_boundary_id,
                       Free_surface_boundary_id,
                       Inner_corner_solid_node_pt);

      // Reset error estimator for bulk mesh
      set_error_estimator();

      // Set the variable and function pointers
      setup_bulk_elements();

      // Reset lagrangian coordinates
      reset_lagrange();
    }

  public:
    void reset_lagrange()
    {
      Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
      if (Free_surface_mesh_pt)
      {
        set_kinematic_lagrange_multiplier(0.0);
      }
    }

  private:
    void create_non_refineable_elements()
    {
      // Create the elements of the other meshes
      // Create the free surface elements
      if (Free_surface_mesh_pt)
      {
        create_free_surface_elements();
      }

      // Create the volume constraint elements
      create_volume_constraint_elements();

      // Need to make the bounding element
      if (Contact_angle_mesh_pt)
      {
        create_contact_angle_element();
      }

      // Need to make the slip boundary elements
      if (Slip_boundary_mesh_pt)
      {
        create_slip_elements(
          Outer_boundary_with_slip_id, Bulk_mesh_pt, Slip_boundary_mesh_pt);
      }

      // Create the no penetration elements on the slip boundary
      if (No_penetration_boundary_mesh_pt)
      {
        create_no_penetration_elements(Outer_boundary_with_slip_id,
                                       Bulk_mesh_pt,
                                       No_penetration_boundary_mesh_pt);
      }

      // Create the flux elements
      create_flux_elements();

      // Create the height elements
      create_height_elements();

      // Create the other meshes
      if (!Augmented_bulk_element_number.empty())
      {
        parameters::x_centre_node_pt = Contact_line_solid_node_pt;
        std::cout << "Make augmented elements" << std::endl;
        create_singularity_scaling_elements();
        create_pressure_contribution_1_elements();
        create_pressure_contribution_2_elements();

        create_slip_eigen_elements(Bulk_mesh_pt);

        // Setup the mesh interaction between the bulk and singularity meshes
        setup_mesh_interaction();
        create_mesh_as_geom_object();
        setup_pressure_contribution_elements();
      }
    }

    // Backup and restore the lagrange multipliers
    void backup_lagrange_multipliers()
    {
      if (No_penetration_boundary_mesh_pt)
      {
        Backed_up_surface_mesh_pt.push_back(
          new BackupMeshForProjection<TElement<1, 3>>(
            No_penetration_boundary_mesh_pt,
            Outer_boundary_with_slip_id,
            Lagrange_id::No_penetration));
        Backup_mesh_index_vector.push_back(Backup_mesh_id::No_penetration);
      }

      if (Free_surface_mesh_pt)
      {
        Backed_up_surface_mesh_pt.push_back(
          new BackupMeshForProjection<TElement<1, 3>>(Free_surface_mesh_pt,
                                                      Free_surface_boundary_id,
                                                      Lagrange_id::Kinematic));
        Backup_mesh_index_vector.push_back(Backup_mesh_id::Free_surface);
      }

      if (Net_flux_mesh_pt)
      {
        Backup_net_flux_lagrange_multiplier =
          dynamic_cast<NET_FLUX_ELEMENT*>(Net_flux_mesh_pt->element_pt(0))
            ->internal_data_pt(0)
            ->value(0);
      }
    }

    void restore_lagrange_multipliers()
    {
      if (Net_flux_mesh_pt)
      {
        dynamic_cast<NET_FLUX_ELEMENT*>(Net_flux_mesh_pt->element_pt(0))
          ->internal_data_pt(0)
          ->set_value(0, Backup_net_flux_lagrange_multiplier);
      }

      if (!Backup_mesh_index_vector.empty())
      {
        for (unsigned n = 0; n < Backup_mesh_index_vector.size(); n++)
        {
          switch (Backup_mesh_index_vector.back())
          {
            case Backup_mesh_id::No_penetration:
              if (No_penetration_boundary_mesh_pt)
              {
                Backed_up_surface_mesh_pt.back()->project_onto_new_mesh(
                  this->No_penetration_boundary_mesh_pt);
              }
              break;
            case Backup_mesh_id::Free_surface:
              Backed_up_surface_mesh_pt.back()->project_onto_new_mesh(
                this->Free_surface_mesh_pt);
              break;
          }
          Backed_up_surface_mesh_pt.pop_back();
          Backup_mesh_index_vector.pop_back();
        }
      }
    }

    double node_position_error_temporal_error()
    {
      // Set the position error to zero
      double node_position_error = 0.0;

      // Find out how many nodes there are in the problem
      const unsigned n_node = Bulk_mesh_pt->nnode();

      // Find number of dimensions of the node
      const unsigned n_dim = Bulk_mesh_pt->node_pt(0)->ndim();

      // Loop over the nodes and calculate the errors in the positions
      for (unsigned n = 0; n < n_node; n++)
      {
        // Loop over the dimensions
        for (unsigned i = 0; i < n_dim; i++)
        {
          // Get position error
          double error =
            Bulk_mesh_pt->node_pt(n)
              ->position_time_stepper_pt()
              ->temporal_error_in_position(Bulk_mesh_pt->node_pt(n), i);

          // Add the square of the individual error to the position error
          node_position_error += error * error;
        }
      }

      // Now scale the errors by the number of dimensions.
      node_position_error /= n_dim;

      // Now scale the errors by the number of nodes.
      node_position_error /= n_node;

      return sqrt(node_position_error);
    }

    double velocity_value_error_temporal_error()
    {
      // Set the velocity value error to zero
      double velocity_value_error = 0.0;

      // Find out how many nodes there are in the problem
      const unsigned n_node = Bulk_mesh_pt->nnode();

      // Find number of dimensions of the node
      const unsigned n_dim = Bulk_mesh_pt->node_pt(0)->ndim();

      // Loop over the nodes and calculate the errors in the positions
      for (unsigned n = 0; n < n_node; n++)
      {
        // Loop over the dimensions
        for (unsigned i = 0; i < n_dim + 1; i++)
        {
          // Get solution error
          double error =
            Bulk_mesh_pt->node_pt(n)
              ->time_stepper_pt()
              ->temporal_error_in_value(Bulk_mesh_pt->node_pt(n), i);

          // Add the square of the individual error to the position error
          velocity_value_error += error * error;
        }
      }

      // Now scale the errors by the number of directions.
      velocity_value_error /= n_dim + 1;

      // Now scale the errors by the number of nodes.
      velocity_value_error /= n_node;

      // Return the global error
      return sqrt(velocity_value_error);
    }

    double global_temporal_error_norm()
    {
      // Return the sum of errors
      return (node_position_error_temporal_error() +
              velocity_value_error_temporal_error());
    }

    double global_velocity_norm()
    {
      // Set the velocity value error to zero
      double velocity_value_norm = 0.0;

      // Find out how many nodes there are in the problem
      const unsigned n_node = Bulk_mesh_pt->nnode();

      // Find number of dimensions of the node
      const unsigned n_dim = Bulk_mesh_pt->node_pt(0)->ndim();

      // Loop over the nodes and calculate the errors in the positions
      for (unsigned n = 0; n < n_node; n++)
      {
        // Loop over the dimensions
        for (unsigned i = 0; i < n_dim + 1; i++)
        {
          // Get solution error
          double value = Bulk_mesh_pt->node_pt(n)->value(i);

          // Add the square of the individual error to the position error
          velocity_value_norm += value * value;
        }
      }

      // Now scale the norm by the number of directions.
      velocity_value_norm /= n_dim + 1;

      // Now scale the norm by the number of nodes.
      velocity_value_norm /= n_node;

      // Return the global error
      return sqrt(velocity_value_norm);
    }

    double global_deformation_norm()
    {
      // Set the velocity value error to zero
      double deformation_norm = 0.0;

      // Find out how many nodes there are in the problem
      const unsigned n_node = Bulk_mesh_pt->nnode();

      // Find number of dimensions of the node
      const unsigned n_dim = Bulk_mesh_pt->node_pt(0)->ndim();

      // Loop over the nodes and calculate the errors in the positions
      for (unsigned n = 0; n < n_node; n++)
      {
        // Loop over the dimensions
        for (unsigned i = 0; i < n_dim; i++)
        {
          // Get solution error
          double value = dynamic_cast<SolidNode*>(Bulk_mesh_pt->node_pt(n))
                           ->lagrangian_position(i);

          // Add the square of the individual value to the position value
          deformation_norm += value * value;
        }
      }

      // Now scale the norm by the number of directions.
      deformation_norm /= n_dim;

      // Now scale the norm by the number of nodes.
      deformation_norm /= n_node;

      // Return the global error
      return sqrt(deformation_norm);
    }


    //============================================================================
    // Destruction functions
    //============================================================================
  public:
    // Destructor: clean up memory allocated by the object
    ~BoHeightControlSingularAxisymDynamicCapProblem()
    {
      // Close the trace file
      close_trace_files();

      // Delete all non-refineable elements
      delete_non_refineable_elements();

      // Delete all pointers created with "new" in reverse order
      delete Height_mesh_pt;
      // Delete all pointers created with "new" in reverse order
      if (Net_flux_mesh_pt)
      {
        delete Net_flux_mesh_pt;
      }
      if (Flux_mesh_pt)
      {
        delete Flux_mesh_pt;
      }
      if (No_penetration_boundary_mesh_pt)
      {
        delete No_penetration_boundary_mesh_pt;
      }
      if (Slip_boundary_mesh_pt)
      {
        delete Slip_boundary_mesh_pt;
      }
      if (Contact_angle_mesh_pt)
      {
        delete Contact_angle_mesh_pt;
      }
      if (Volume_computation_mesh_pt)
      {
        delete Volume_computation_mesh_pt;
      }
      if (Volume_constraint_mesh_pt)
      {
        delete Volume_constraint_mesh_pt;
      }
      if (Free_surface_mesh_pt)
      {
        delete Free_surface_mesh_pt;
      }
      if (Height_mesh_pt)
      {
        delete_elements(Height_mesh_pt);
      }

      // Next delete the external data
      delete External_pressure_data_pt;

      delete Bulk_mesh_pt;

      delete Constitutive_law_pt;
    }

  private:
    // Delete elements of the passed mesh
    void delete_elements(Mesh* mesh_pt)
    {
      if (mesh_pt)
      {
        // Delete the slip elements
        unsigned n_element = mesh_pt->nelement();
        for (unsigned e = 0; e < n_element; e++)
        {
          delete mesh_pt->element_pt(e);
        }
        // Now flush the storage
        mesh_pt->flush_element_and_node_storage();
      }
    }

  private:
    //============================================================================
    // Class Variables
    //============================================================================
    // Dimensionless parameters
    double Ca;
    double Bo;
    double Re;
    double ReSt;
    double Viscosity_ratio;
    double St;
    Data* ReInvFr_data_pt;

    // The prescribed volume of the fluid if static
    double Volume;

    // The external pressure
    double Pext;

    // The contact angle
    double Contact_angle;
    double Right_angle;
    double Zero_sigma;

    // Maximum number of adapt steps
    unsigned Max_adapt;

    // Constitutive law used to determine the mesh deformation
    ConstitutiveLaw* Constitutive_law_pt;

    // Data object whose single value stores the external pressure
    Data* External_pressure_data_pt;

    Data* Wall_velocity_data_pt;
    Data* Height_step_data_pt;

    // Data object whose single value stores the continuation parameter
    Data* Parameter_data_pt;

    // Contact line node pointer
    SolidNode* Contact_line_solid_node_pt;
    SolidNode* Inner_corner_solid_node_pt;

    // Error Estimator pointers
    Z2ErrorEstimator* Z2_error_estimator_pt;
    ErrorEstimator* Corner_error_estimator_pt;

    // Trace file
    std::ofstream Trace_file;
    std::ofstream Volume_trace_file;
    std::ofstream Flux_trace_file;
    std::ofstream Contact_angle_trace_file;
    std::ofstream Inner_angle_trace_file;

    std::chrono::high_resolution_clock::time_point Start_time;

    // Storage for the bulk mesh
    RefineableSolidTriangleMesh<ELEMENT>* Bulk_mesh_pt;


  private:
    // Storage for the free surface mesh
    Mesh* Free_surface_mesh_pt;

    // Storage for the slip elements
    Mesh* Slip_boundary_mesh_pt;

    // Storage for the no penetration elements
    Mesh* No_penetration_boundary_mesh_pt;

    // Storage for the flux elements
    Mesh* Flux_mesh_pt;

    // Storage for the elements that compute the enclosed volume
    Mesh* Volume_computation_mesh_pt;

    // Storage for the element bounding the free surface
    Mesh* Contact_angle_mesh_pt;

    // Storage for the augmented element
    Mesh* Eigensolution_slip_mesh_pt;
    Mesh* Singularity_scaling_mesh_pt;
    Mesh* Pressure_contribution_mesh_1_pt;
    Mesh* Pressure_contribution_mesh_2_pt;
    MeshAsGeomObject* Pressure_contribution_geom_mesh_1_pt;
    MeshAsGeomObject* Pressure_contribution_geom_mesh_2_pt;

    // Storage for the volume constraint
    Mesh* Volume_constraint_mesh_pt;

    // Storage for the net flux element
    Mesh* Net_flux_mesh_pt;

    // Storage for the height element
    Mesh* Height_mesh_pt;

    // Backup mesh
    Vector<BackupMeshForProjection<TElement<1, 3>>*> Backed_up_surface_mesh_pt;

    // List of augmented element numbers
    Vector<unsigned> Augmented_bulk_element_number;

    // Backup mesh id
    enum struct Backup_mesh_id
    {
      No_penetration,
      Free_surface,
    };

    enum Lagrange_id
    {
      No_penetration,
      Kinematic,
    };

    // Backup mesh indices
    Vector<Backup_mesh_id> Backup_mesh_index_vector;

    // Triangle mesh polygon for outer boundary
    TriangleMeshPolygon* Outer_boundary_polyline_pt;

    // Backup lagrange multipliers
    double Backup_volume_constraint_lagrange_multiplier;
    double Backup_net_flux_lagrange_multiplier;
    double Backup_point_kinematic_lagrange_multiplier;

    // Enumeration of domain boundaries
    enum
    {
      Upper_boundary_id,
      Outer_boundary_with_slip_id,
      Free_surface_boundary_id,
      Inner_boundary_id,
    };

    // Enumeration of velocity indices
    enum
    {
      u_index,
      v_index,
      w_index
    };

    // Enumeration of the different problems
    enum Problem_type
    {
      Normal_problem,
      Bulk_only_problem
    };

    Problem_type problem_type;


    // Static problem state Boolean
    bool Is_steady;
    bool Using_contact_angle_error_estimator;
  };
} // namespace oomph
#endif
