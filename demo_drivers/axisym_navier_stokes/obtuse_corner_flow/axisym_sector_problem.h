#ifndef AXISYM_SECTOR_PROBLEM_HEADER
#define AXISYM_SECTOR_PROBLEM_HEADER

#include "generic.h"
#include "solid.h"
//#include "constitutive.h"
#include "navier_stokes.h"
#include "fluid_interface.h"
#include "meshes/triangle_mesh.h"

#include "free_surface_element.h"
#include "far_field_element.h"
#include "domain_boundaries.h"
#include "parameters.h"
#include "parameter_values.h"
#include "parameter_functions.h"
#include "my_error_estimator.h"

#include "fluid_slip_elements.h"

namespace oomph
{
  // Problem class
  template<class ELEMENT>
  class AxisymSectorProblem : public Problem
  {
  public:
    // Constructor
    AxisymSectorProblem() : Contact_line_solid_node_pt(0)
    {
      // Create and add the timestepper
      add_time_stepper_pt(new BDF<2>);

      // Create an empty mesh
      add_bulk_mesh();
      add_non_adaptive_sub_meshes();
      build_global_mesh();

      // Assign doc info pointer
      Doc_info_pt = new DocInfo("RESLT");
      Doc_info_pt->number() = 0;

      // From here should be actions_after_adapt

      // Omega
      setup_bulk_elements();

      // S2
      // create_far_field_elements(Bulk_mesh_pt);

      // S0
      create_no_penetration_elements(Bulk_mesh_pt);
      create_slip_elements(Bulk_mesh_pt);

      // L0

      rebuild_global_mesh();
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;

      // Remove solid mesh equations and azimuthal velocity
      pin_solid_dofs();

      // Set the boundary conditions
      set_boundary_conditions();

      // Setup all the equation numbering and look-up schemes
      rebuild_global_mesh();
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;

      // Document the initial conditions
      doc_solution();
    }

  private:
    RefineableSolidTriangleMesh<ELEMENT>* Bulk_mesh_pt;
    SolidNode* Contact_line_solid_node_pt;

    Mesh* No_penetration_boundary_mesh_pt;
    Mesh* Far_field_mesh_pt;
    Mesh* Slip_boundary_mesh_pt;
    Mesh* Contact_angle_mesh_pt;

    DocInfo* Doc_info_pt;

    // Enumeration of Lagrangian identities
    enum Lagrange_id
    {
      Free_surface,
      No_penetration,
      Far_field,
    };

  public:
    void add_non_adaptive_sub_meshes()
    {
      No_penetration_boundary_mesh_pt = new Mesh;
      add_sub_mesh(No_penetration_boundary_mesh_pt);
      Far_field_mesh_pt = new Mesh;
      add_sub_mesh(Far_field_mesh_pt);
      Slip_boundary_mesh_pt = new Mesh;
      add_sub_mesh(Slip_boundary_mesh_pt);
      Contact_angle_mesh_pt = new Mesh;
      add_sub_mesh(Contact_angle_mesh_pt);
    }

    void setup_bulk_elements()
    {
      set_contact_line_node_pt(Bulk_mesh_pt);

      unsigned n_bulk = Bulk_mesh_pt->nelement();
      for (unsigned e = 0; e < n_bulk; e++)
      {
        // Upcast from GeneralisedElement to the present element
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

        // Set the Reynolds number
        el_pt->re_pt() = &parameters::reynolds_number;

        // Set the Reynolds Strouhal number
        el_pt->re_st_pt() = &parameters::strouhal_reynolds_number;
      }
    }

    // Destructor
    ~AxisymSectorProblem()
    {
      const unsigned n_sub_mesh = nsub_mesh();
      for (unsigned i = 0; i < n_sub_mesh; i++)
      {
        delete_elements(this->mesh_pt(i));
        delete this->mesh_pt(i);
        this->mesh_pt(i) = 0;
      }
      delete Doc_info_pt;
      Doc_info_pt = 0;
    }

    // Delete the created elements
    void delete_elements(Mesh* local_mesh_pt)
    {
      // Delete the slip elements
      unsigned n_element = local_mesh_pt->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        delete local_mesh_pt->element_pt(e);
        local_mesh_pt->element_pt(e) = 0;
      }
      // Now flush the storage
      local_mesh_pt->flush_element_and_node_storage();
    }

    void doc_solution();

    void debug_residuals()
    {
      DoubleVector res;
      oomph_info << "get_residuals" << endl;
      get_residuals(res);
      res.output("res.dat");
    }

    void debug_jacobian()
    {
      oomph_info << "output dofs" << endl;
      std::ofstream output_stream("dofs.txt");
      describe_dofs(output_stream);
      output_stream.close();

      CRDoubleMatrix jac;
      DoubleVector res;
      oomph_info << "get_jacobian" << endl;
      get_jacobian(res, jac);
      jac.sparse_indexed_output("jac.dat", 15, 1);

      unsigned n_dof = ndof();
      DenseDoubleMatrix fd_jac(n_dof);
      oomph_info << "get_fd_jacobian" << endl;
      get_fd_jacobian(res, fd_jac);
      fd_jac.sparse_indexed_output("fd_jac.dat", 15, 1);
    }
    void make_steady()
    {
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(8));
      el_pt->fix_pressure(0, 0.0);
      Vector<double> s(2, 0.0);
      Vector<double> x(2, 0.0);
      el_pt->get_x(s, x);
      oomph_info << "Pressure point, ";
      oomph_info << "x: " << x[0] << ", ";
      oomph_info << "y: " << x[1] << ", ";
      oomph_info << endl;
    }

    void make_unsteady()
    {
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(8));
      el_pt->free_pressure(0);
      Vector<double> s(2, 0.0);
      Vector<double> x(2, 0.0);
      el_pt->get_x(s, x);
      oomph_info << "Pressure point, ";
      oomph_info << "x: " << x[0] << ", ";
      oomph_info << "y: " << x[1] << ", ";
      oomph_info << endl;
    }

  private:
    void add_bulk_mesh();
    void create_slip_elements(Mesh* const& bulk_mesh_pt);
    void create_no_penetration_elements(Mesh* const& bulk_mesh_pt);
    void create_contact_angle_elements(Mesh* const& bulk_mesh_pt,
                                       Mesh* const& surface_mesh_pt);
    void create_far_field_elements(Mesh* const& bulk_mesh_pt);

  public:
    void pin_far_field_elements();

  private:
    void refine_mesh_for_weak_contact_angle_constraint();
    void set_contact_line_node_pt(const Mesh* const bulk_mesh_pt);
    void compute_error_estimate(double& max_err, double& min_err);
    void find_corner_bulk_element(const unsigned& boundary_1_id,
                                  const unsigned& boundary_2_id,
                                  unsigned& element_index);
    void find_corner_bulk_node(const unsigned& boundary_1_id,
                               const unsigned& boundary_2_id,
                               unsigned& element_index,
                               unsigned& node_index);

    void setup_pressure_contribution_elements();

  public:
    void pin_solid_dofs()
    {
      oomph_info << "pin_solid_dofs" << endl;

      // Set all azimuthal velocities to zero
      unsigned n_node = Bulk_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        Bulk_mesh_pt->node_pt(n)->pin_position(0);
        Bulk_mesh_pt->node_pt(n)->pin_position(1);
      }


      /* We don't have a contact angle element
      // Fix the extra kinematic lagrange_multiplier of the contact angle
      // point
      dynamic_cast<PointFluidInterfaceBoundingElement*>(
        Contact_angle_mesh_pt->element_pt(0))
        ->fix_lagrange_multiplier(0.0);
      */

      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void unpin_solid_dofs()
    {
      oomph_info << "pin_solid_dofs" << endl;

      // Set all azimuthal velocities to zero
      unsigned n_node = Bulk_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        Bulk_mesh_pt->node_pt(n)->unpin_position(0);
        Bulk_mesh_pt->node_pt(n)->unpin_position(1);
      }


      /* We don't have a contact angle element
      // Fix the extra kinematic lagrange_multiplier of the contact angle
      // point
      dynamic_cast<PointFluidInterfaceBoundingElement*>(
        Contact_angle_mesh_pt->element_pt(0))
        ->fix_lagrange_multiplier(0.0);
      */

      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

  private:
    void set_boundary_conditions()
    {
      oomph_info << "set_boundary_conditions" << endl;

      // Set fluid boundary conditions
      // set_dirichlet_bc_on_boundary(Far_field_boundary_id);
      // set_dirichlet_bc_on_boundary(Slip_boundary_id);
      // set_dirichlet_bc_on_boundary(Free_surface_boundary_id);

      // pin_free_surface_contact_line_point();

      // Fix end points far field boundary condition lagrange_multipliers
      // pin_far_field_lagrange_multiplier_end_points();
    }

    void pin_no_penetration_contact_line_point()
    {
      // Fix end points for the no pen boundary condition lagrange_multipliers
      const unsigned n_el = No_penetration_boundary_mesh_pt->nelement();
      const double zero_value = 0.0;
      for (unsigned i_el = 0; i_el < n_el; i_el++)
      {
        ImposeImpenetrabilityElement<ELEMENT>* el_pt =
          dynamic_cast<ImposeImpenetrabilityElement<ELEMENT>*>(
            No_penetration_boundary_mesh_pt->element_pt(i_el));
        const unsigned n_nod = el_pt->nnode();
        for (unsigned i_nod = 0; i_nod < n_nod; i_nod++)
        {
          // Get boundary node
          const Node* const node_pt = el_pt->node_pt(i_nod);
          // If node is on either of the other boundaries
          if (node_pt->is_on_boundary(Slip_boundary_id) &&
              node_pt->is_on_boundary(Free_surface_boundary_id))
          {
            // then fix the lagrange multiplier to zero
            oomph_info << "Fix lagrange_multiplier" << std::endl;
            el_pt->fix_lagrange_multiplier(i_nod, 0.0);
          }
        }
      }
    }


    void pin_far_field_lagrange_multiplier_end_points()
    {
      // Non-augmented
      {
        const unsigned n_el = Far_field_mesh_pt->nelement();
        const double zero_value = 0.0;
        for (unsigned i_el = 0; i_el < n_el; i_el++)
        {
          FarFieldElement<ELEMENT>* el_pt =
            dynamic_cast<FarFieldElement<ELEMENT>*>(
              Far_field_mesh_pt->element_pt(i_el));
          const unsigned n_nod = el_pt->nnode();
          for (unsigned i_nod = 0; i_nod < n_nod; i_nod++)
          {
            // Get boundary node
            const Node* const node_pt = el_pt->node_pt(i_nod);
            // If node is on either of the other boundaries
            if (node_pt->is_on_boundary(Slip_boundary_id) ||
                node_pt->is_on_boundary(Free_surface_boundary_id))
            {
              // then fix the lagrange multiplier to zero
              oomph_info << "Fix lagrange_multiplier" << std::endl;
              el_pt->fix_lagrange_multiplier(i_nod, 0, zero_value);
              el_pt->fix_lagrange_multiplier(i_nod, 1, zero_value);
            }
          }
        }
      }
    }
  };

  template<class ELEMENT>
  void AxisymSectorProblem<ELEMENT>::add_bulk_mesh()
  {
    oomph_info << "add_bulk_mesh" << endl;

    // Generate triangle_mesh_parameters
    // ---------------------------------

    // Setup the sector boundaries
    Vector<TriangleMeshCurveSection*> boundary_polyline_pt;
    create_sector_domain(parameters::contact_angle, boundary_polyline_pt);

    // Create an interior point: Assuming the contact angle is alpha > 1e-2
    // and the radius is > 0.5
    Vector<double> internal_point_pt(2);
    internal_point_pt[0] = 2 - 0.5 * sin(1e-2 * MathematicalConstants::Pi);
    internal_point_pt[1] = 0.5 * cos(1e-2 * MathematicalConstants::Pi);

    // Set up a bool for whether the internal point is fixed
    const bool is_internal_point_fixed = true;

    // True the curve sections into a closed curve
    TriangleMeshClosedCurve* outer_closed_curve_pt =
      new TriangleMeshClosedCurve(
        boundary_polyline_pt, internal_point_pt, is_internal_point_fixed);

    // Create the Triangle mesh generator parameters
    // WARNING: outer_closed_curve_pt's destruction is handled by
    // TriangleMeshParameters. This is unexpected behaviour.
    TriangleMeshParameters triangle_mesh_parameters(outer_closed_curve_pt);

    triangle_mesh_parameters.element_area() = parameters::element_area;

    // ---------------------------------

    // Generate the mesh using the template ELEMENT
    Bulk_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(
      triangle_mesh_parameters, this->time_stepper_pt());

    // Add mesh to problem
    add_sub_mesh(Bulk_mesh_pt);

    // refine_mesh_for_weak_contact_angle_constraint();
  }

  template<class ELEMENT>
  void AxisymSectorProblem<ELEMENT>::create_slip_elements(
    Mesh* const& bulk_mesh_pt)
  {
    oomph_info << "create_slip_elements" << endl;

    unsigned b = Slip_boundary_id;

    // How many bulk elements are adjacent to boundary b?
    unsigned n_element = bulk_mesh_pt->nboundary_element(b);

    // Loop over the bulk elements adjacent to boundary b?
    for (unsigned e = 0; e < n_element; e++)
    {
      // Get pointer to the bulk element that is adjacent to boundary
      // b
      ELEMENT* bulk_elem_pt =
        dynamic_cast<ELEMENT*>(bulk_mesh_pt->boundary_element_pt(b, e));

      // What is the index of the face of element e along boundary b
      int face_index = bulk_mesh_pt->face_index_at_boundary(b, e);
      NavierStokesSlipElement<ELEMENT>* slip_element_pt = 0;

      // Build the corresponding slip element
      slip_element_pt =
        new NavierStokesSlipElement<ELEMENT>(bulk_elem_pt, face_index);

      // Set the pointer to the prescribed slip function
      slip_element_pt->slip_fct_pt() = &parameters::prescribed_slip_fct;

      slip_element_pt->wall_velocity_fct_pt() =
        &parameters::prescribed_wall_velocity_fct;

      // Add the prescribed-flux element to the surface mesh
      Slip_boundary_mesh_pt->add_element_pt(slip_element_pt);
    }
  }

  template<class ELEMENT>
  void AxisymSectorProblem<ELEMENT>::create_no_penetration_elements(
    Mesh* const& bulk_mesh_pt)
  {
    oomph_info << "create_no_penetration_elements" << endl;

    unsigned b = Slip_boundary_id;

    // How many bulk elements are adjacent to boundary b?
    unsigned n_element = bulk_mesh_pt->nboundary_element(b);

    // Loop over the bulk elements adjacent to boundary b?
    for (unsigned e = 0; e < n_element; e++)
    {
      // Get pointer to the bulk element that is adjacent to boundary
      // b
      ELEMENT* bulk_elem_pt =
        dynamic_cast<ELEMENT*>(bulk_mesh_pt->boundary_element_pt(b, e));

      // What is the index of the face of element e along boundary b
      int face_index = bulk_mesh_pt->face_index_at_boundary(b, e);

      // Build the corresponding slip element
      ImposeImpenetrabilityElement<ELEMENT>* no_penetration_element_pt =
        new ImposeImpenetrabilityElement<ELEMENT>(
          bulk_elem_pt, face_index, Lagrange_id::No_penetration);

      // Add the prescribed-flux element to the surface mesh
      No_penetration_boundary_mesh_pt->add_element_pt(
        no_penetration_element_pt);
    }
  }

  template<class ELEMENT>
  void AxisymSectorProblem<ELEMENT>::create_contact_angle_elements(
    Mesh* const& bulk_mesh_pt, Mesh* const& surface_mesh_pt)
  {
    oomph_info << "create_contact_angle_elements" << endl;

    /* WHAT I WANT
    const unsigned surface_boundary_id = 0;
    const unsigned n_free_surface =
      surface_mesh_pt->nboundary_element(surface_boundary_id);
    for (unsigned e = 0; e < n_free_surface; e++)
    {
      // Get pointer to the bulk fluid element that is
      // adjacent to boundary b
      MyFreeSurfaceElement* surface_elem_pt =
        dynamic_cast<MyFreeSurfaceElement*>(
          surface_mesh_pt->boundary_element_pt(surface_boundary_id, e));

      // Find the index of the face of element e along boundary b
      int face_index =
        surface_mesh_pt->face_index_at_boundary(surface_boundary_id, e);

      // Locally cache the element pointer
      ElasticPointFluidInterfaceBoundingElement<MyFreeSurfaceElement> el_pt =
        new ElasticPointFluidInterfaceBoundingElement<MyFreeSurfaceElement>(
          surface_elem_pt, face_index);

      // Set the contact angle function
      const bool use_strong_imposition = false;
      el_pt->set_contact_angle_fct(&parameters::contact_angle_fct_pt,
                                   use_strong_imposition);

      // Set the capillary number
      el_pt->ca_pt() = &parameters::capillary_number;

      // Set the wall normal of the external boundary
      el_pt->wall_unit_normal_fct_pt() = &parameters::wall_unit_normal_fct;

      // Add the element to the mesh
      Contact_angle_mesh_pt->add_element_pt(el_pt);
    }
    */

    /* WHAT I HAVE
    // Inialise storage for bounding element
    FluidInterfaceBoundingElement* el_pt = 0;

    // Here we have no guarantee of order so we need to loop over all
    // surface elements to find the one that is next to the outer boundary
    unsigned n_free_surface = surface_mesh_pt->nelement();
    for (unsigned e = 0; e < n_free_surface; e++)
    {
      // Locally cache the element pointer
      // ElasticLineFluidInterfaceElement<ELEMENT>* bulk_el_pt =
      //   dynamic_cast<ElasticLineFluidInterfaceElement<ELEMENT>*>(
      //     surface_mesh_pt->element_pt(e));
      MyFreeSurfaceElement<ELEMENT>* bulk_el_pt =
        dynamic_cast<MyFreeSurfaceElement<ELEMENT>*>(
          surface_mesh_pt->element_pt(e));

      // Read out number of nodes in the element
      unsigned n_node = bulk_el_pt->nnode();

      // Is the "left" hand node on the boundary
      if (bulk_el_pt->node_pt(0)->is_on_boundary(Slip_boundary_id))
      {
        // Create bounding element on "left" hand face, with the normal
        // pointing into the fluid
        el_pt = bulk_el_pt->make_bounding_element(-1);

        // Exit loop
        break;
      }

      // Is the "right" hand node on the boundary
      if (bulk_el_pt->node_pt(n_node - 1)->is_on_boundary(Slip_boundary_id))
      {
        // Create bounding element on "right" hand face, with the normal
        // pointing out of the fluid
        el_pt = bulk_el_pt->make_bounding_element(1);

        // Exit loop
        break;
      }
    }

    // Set the contact angle function
    const bool use_strong_imposition = false;
    el_pt->set_contact_angle_fct(&parameters::contact_angle_fct_pt,
                                 use_strong_imposition);

    // Set the capillary number
    el_pt->ca_pt() = &parameters::capillary_number;

    // Set the wall normal of the external boundary
    el_pt->wall_unit_normal_fct_pt() = &parameters::wall_unit_normal_fct;

    // Add the element to the mesh
    Contact_angle_mesh_pt->add_element_pt(el_pt);
    */
  }

  template<class ELEMENT>
  void AxisymSectorProblem<
    ELEMENT>::refine_mesh_for_weak_contact_angle_constraint()
  {
    oomph_info << "Refining the mesh about the contact line" << endl;

    // Set the refinement tolerances
    Bulk_mesh_pt->max_element_size() = 0.5 * pow(2e-1, 2.0);
    Bulk_mesh_pt->min_element_size() =
      0.5 * pow(Mesh_Control_Parameters::min_element_length, 2.0);
    Bulk_mesh_pt->min_permitted_angle() = 15;
    Bulk_mesh_pt->min_permitted_error() = 1.0e-2;
    Bulk_mesh_pt->max_permitted_error() = 1.0;

    double max_error = 1e5;
    double min_error = 0;
    const unsigned max_n_adapt = 20;
    // Call the adapt function until the maximum of the estimated error is
    // below 1e0 ( = Bulk_mesh->max_permitted_error)
    for (unsigned i_adapt = 0; i_adapt < max_n_adapt; i_adapt++)
    {
      // Check refinement

      // Find the contact line node
      set_contact_line_node_pt(Bulk_mesh_pt);

      // Get error estimator
      Bulk_mesh_pt->spatial_error_estimator_pt() =
        new ContactlineErrorEstimator(
          Contact_line_solid_node_pt,
          Mesh_Control_Parameters::min_element_length,
          Mesh_Control_Parameters::element_length_ratio);
      compute_error_estimate(max_error, min_error);
      oomph_info << "Number of adaptions: " << i_adapt
                 << ", Max error: " << max_error << endl;

      // If we haven't refined enough ...
      if (max_error > 1e0)
      {
        //... refine mesh,
        this->adapt();
      }
      else
      {
        // ... else, end loop.
        break;
      }
    }
  }

  template<class ELEMENT>
  void AxisymSectorProblem<ELEMENT>::set_contact_line_node_pt(
    const Mesh* const bulk_mesh_pt)
  {
    oomph_info << "set_contact_line_node_pt" << endl;

    unsigned element_index;
    unsigned node_index;

    find_corner_bulk_node(
      Slip_boundary_id, Free_surface_boundary_id, element_index, node_index);

    Contact_line_solid_node_pt = dynamic_cast<SolidNode*>(
      bulk_mesh_pt->boundary_element_pt(Slip_boundary_id, element_index)
        ->node_pt(node_index));
  }

  template<class ELEMENT>
  void AxisymSectorProblem<ELEMENT>::compute_error_estimate(double& max_err,
                                                            double& min_err)
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

  template<class ELEMENT>
  void AxisymSectorProblem<ELEMENT>::create_far_field_elements(
    Mesh* const& bulk_mesh_pt)
  {
    oomph_info << "create_far_field_elements" << endl;

    // Loop over the free surface boundary and create the "interface elements
    unsigned b = Far_field_boundary_id;

    // How many bulk fluid elements are adjacent to boundary b?
    unsigned n_element = bulk_mesh_pt->nboundary_element(b);

    // Loop over the bulk fluid elements adjacent to boundary b?
    for (unsigned e = 0; e < n_element; e++)
    {
      // Get pointer to the bulk fluid element that is
      // adjacent to boundary b
      ELEMENT* bulk_elem_pt =
        dynamic_cast<ELEMENT*>(bulk_mesh_pt->boundary_element_pt(b, e));

      if (!bulk_elem_pt->is_augmented())
      {
        // Find the index of the face of element e along boundary b
        int face_index = bulk_mesh_pt->face_index_at_boundary(b, e);

        // Create new element
        FarFieldElement<ELEMENT>* el_pt = new FarFieldElement<ELEMENT>(
          bulk_elem_pt, face_index, Lagrange_id::Far_field);

        // Add it to the mesh
        Far_field_mesh_pt->add_element_pt(el_pt);
      }
    }
  }

  template<class ELEMENT>
  void AxisymSectorProblem<ELEMENT>::pin_far_field_elements()
  {
    const unsigned n_el = Far_field_mesh_pt->nelement();
    for (unsigned i_el = 0; i_el < n_el; i_el++)
    {
      FarFieldElement<ELEMENT>* el_pt = dynamic_cast<FarFieldElement<ELEMENT>*>(
        Far_field_mesh_pt->element_pt(i_el));
      const unsigned n_nod = el_pt->nnode();
      for (unsigned i_nod = 0; i_nod < n_nod; i_nod++)
      {
        const unsigned n_eq = 2;
        for (unsigned i_eq = 0; i_eq < n_eq; i_eq++)
        {
          el_pt->pin_lagrange_multiplier(i_nod, i_eq);
        }
      }
    }
  }

  template<class ELEMENT>
  void AxisymSectorProblem<ELEMENT>::find_corner_bulk_element(
    const unsigned& boundary_1_id,
    const unsigned& boundary_2_id,
    unsigned& element_index)
  {
    unsigned dummy_node_index;
    find_corner_bulk_node(
      boundary_1_id, boundary_2_id, element_index, dummy_node_index);
  }

  template<class ELEMENT>
  void AxisymSectorProblem<ELEMENT>::find_corner_bulk_node(
    const unsigned& boundary_1_id,
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
        if (bulk_el_pt->node_pt(i_node)->is_on_boundary(boundary_2_id))
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
    oomph_info << "Warning: No corner node found!" << endl;

    return;
  }

  template<class ELEMENT>
  void AxisymSectorProblem<ELEMENT>::doc_solution()
  {
    unsigned doc_number = Doc_info_pt->number();

    oomph_info << "Doc Number: " << doc_number << endl;

    std::ofstream output_stream;
    char filename[100];
    unsigned npts = 3;

    if (doc_number == 0)
    {
      sprintf(filename, "%s/dofs.txt", Doc_info_pt->directory().c_str());
      output_stream.open(filename);
      describe_dofs(output_stream);
      output_stream.close();
    }

    sprintf(filename,
            "%s/soln%i.dat",
            Doc_info_pt->directory().c_str(),
            Doc_info_pt->number());
    output_stream.open(filename);
    output_stream.precision(15);
    Bulk_mesh_pt->output(output_stream, npts);
    output_stream.close();

    //    sprintf(filename,
    //            "%s/pressure_1_%i.csv",
    //            Doc_info_pt->directory().c_str(),
    //            Doc_info_pt->number());
    //    output_stream.open(filename);
    //    output_stream << "x,y,p," << endl;
    //    Pressure_contribution_mesh_1_pt->output(output_stream, npts);
    //    output_stream.close();
    //
    //    sprintf(filename,
    //            "%s/pressure_2_%i.csv",
    //            Doc_info_pt->directory().c_str(),
    //            Doc_info_pt->number());
    //    output_stream.open(filename);
    //    output_stream << "x,y,p," << endl;
    //    Pressure_contribution_mesh_2_pt->output(output_stream, npts);
    //    output_stream.close();
    //
    //     sprintf(filename,
    //             "%s/slip_surface%i.csv",
    //             Doc_info_pt->directory().c_str(),
    //             Doc_info_pt->number());
    //     output_stream.open(filename);
    //     output_stream << "x,y,l_x,l_y,n_x,n_y,u,v" << endl;
    //     Slip_boundary_mesh_pt->output(output_stream, npts);
    //     output_stream.close();
    //
    //     // ImposeImpenetrabilityElement hasn't implemented an output
    //     function
    //     // sprintf(filename,
    //     //         "%s/no_penetration_surface%i.csv",
    //     //         Doc_info_pt->directory().c_str(),
    //     //         Doc_info_pt->number());
    //     // output_stream.open(filename);
    //     // output_stream << "x,y,u,v,p,lagrange_multiplier" << endl;
    //     // No_penetration_boundary_mesh_pt->output(output_stream, 3);
    //     // output_stream.close();
    //
    //     // sprintf(filename,
    //     //        "%s/contact_angle%i.csv",
    //     //        Doc_info_pt->directory().c_str(),
    //     //        Doc_info_pt->number());
    //     // output_stream.open(filename);
    //     // output_stream << "x,y,imposed_angle,computed_angle,lambda" <<
    //     endl;
    //     // Contact_angle_mesh_pt->output(output_stream);
    //     // output_stream.close();
    //
    //     sprintf(filename,
    //             "%s/far_field%i.dat",
    //             Doc_info_pt->directory().c_str(),
    //             Doc_info_pt->number());
    //     output_stream.open(filename);
    //     output_stream
    //       <<
    //       "x,y,n_x,n_y,u_x,u_y,du_xdr,du_ydr,du_xdtheta,du_ydtheta,l_x,l_y"
    //       << endl;
    //     Far_field_mesh_pt->output(output_stream);
    //     output_stream.close();

    Doc_info_pt->number()++;
  }
} // namespace oomph
#endif
