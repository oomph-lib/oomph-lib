#ifndef AXISYM_SECTOR_PROBLEM_HEADER
#define AXISYM_SECTOR_PROBLEM_HEADER

#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"
#include "meshes/triangle_mesh.h"

/// Local headers
#include "fluid_slip_elements.h"
#include "debug_impose_impenetratibility_elements.h"
#include "domain_boundaries.h"
#include "eigensolution_elements.h"
#include "eigensolution_functions.h"
#include "far_field_element.h"
#include "free_surface_element.h"
#include "my_error_estimator.h"
#include "parameter_functions.h"
#include "parameter_values.h"
#include "parameters.h"
#include "pressure_evaluation_elements.h"
//#include "projectable_Ttaylor_hood_elements.h"
#include "refined_sector_tri_mesh.template.h"
#include "refined_sector_tri_mesh.template.cc"
#include "singular_fluid_traction_elements.h"
#include "utility_functions.h"

namespace oomph
{
  // Problem class
  template<class ELEMENT>
  class SectorProblem : public Problem
  {
  private:
    /// Private variables
    RefinedSectorTriMesh<ELEMENT>* Bulk_mesh_pt;
    Node* Contact_line_node_pt;

    Mesh* No_penetration_boundary_mesh_pt;
    Mesh* Far_field_mesh_pt;
    Mesh* Slip_boundary_mesh_pt;

    DocInfo Doc_info;
    Z2ErrorEstimator* Z2_error_estimator_pt;

  public:
    // Constructor
    SectorProblem(const unsigned& n_radial = 10,
                  const unsigned& n_azimuthal = 10)
      : Contact_line_node_pt(0), Z2_error_estimator_pt(new Z2ErrorEstimator)
    {
      // Create and add the timestepper
      add_time_stepper_pt(new BDF<2>);

      // Assign doc info pointer
      Doc_info.set_directory("RESLT");
      Doc_info.number() = 0;

      // Create an empty mesh
      add_bulk_mesh(n_radial, n_azimuthal);
      add_non_adaptive_sub_meshes();
      build_global_mesh();

      // From here should be actions_after_adapt
      actions_after_adapt();
    }

    // Actions before adapt
    void actions_after_adapt()
    {
      setup_bulk_elements();
      create_nonrefineable_elements();
      set_boundary_conditions();

      /// Is this needed?
      rebuild_global_mesh();
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void create_nonrefineable_elements()
    {
      create_slip_elements();
      create_no_penetration_elements(Slip_boundary_id);
      create_no_penetration_elements(Free_surface_boundary_id);
      create_far_field_elements();
    }

    void actions_before_adapt()
    {
      delete_nonrefineable_elements();
    }

    void delete_nonrefineable_elements()
    {
      delete_elements(No_penetration_boundary_mesh_pt);
      delete_elements(Slip_boundary_mesh_pt);
      delete_elements(Far_field_mesh_pt);
    }

    void add_non_adaptive_sub_meshes()
    {
      No_penetration_boundary_mesh_pt = new Mesh;
      add_sub_mesh(No_penetration_boundary_mesh_pt);
      Far_field_mesh_pt = new Mesh;
      add_sub_mesh(Far_field_mesh_pt);
      Slip_boundary_mesh_pt = new Mesh;
      add_sub_mesh(Slip_boundary_mesh_pt);
    }

    void setup_bulk_elements()
    {
      set_contact_line_node_pt();

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
    ~SectorProblem()
    {
      const unsigned n_sub_mesh = nsub_mesh();
      for (unsigned i = 0; i < n_sub_mesh; i++)
      {
        delete_elements(this->mesh_pt(i));
        delete this->mesh_pt(i);
        this->mesh_pt(i) = 0;
      }
      if (Z2_error_estimator_pt)
      {
        delete Z2_error_estimator_pt;
      }
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
      oomph_info << "debug_jacobian" << std::endl;

      DoubleVector residuals;
      DenseDoubleMatrix jacobian;
      DoubleVector residualsFD;
      DenseDoubleMatrix jacobianFD(ndof());

      get_jacobian(residuals, jacobian);
      jacobian.sparse_indexed_output(
        this->Doc_info.directory() + "/jacJ.dat", 16, true);
      get_fd_jacobian(residualsFD, jacobianFD);
      jacobianFD.sparse_indexed_output(
        this->Doc_info.directory() + "/jacfdJ.dat", 16, true);

      bool jacobians_are_equal = compare_matrices(jacobian, jacobianFD);

      if (jacobians_are_equal)
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
        output_stream.open(this->Doc_info.directory() + "/dofs.txt");
        this->describe_dofs(output_stream);
        output_stream.close();
      }
    }

  private:
    void add_bulk_mesh(const unsigned& n_radial, const unsigned& n_azimuthal);
    void create_slip_elements();
    void create_no_penetration_elements(const unsigned& boundary_id);
    void create_far_field_elements();

  public:
    void pin_far_field_elements();

  private:
    void refine_mesh_for_weak_contact_angle_constraint();
    void set_contact_line_node_pt();
    void compute_error_estimate(double& max_err, double& min_err);
    void find_corner_bulk_element(const unsigned& boundary_1_id,
                                  const unsigned& boundary_2_id,
                                  unsigned& element_index);
    void find_corner_bulk_node(const unsigned& boundary_1_id,
                               const unsigned& boundary_2_id,
                               unsigned& element_index,
                               unsigned& node_index);

    void setup_pressure_contribution_elements();

  private:
    void set_boundary_conditions()
    {
      oomph_info << "set_boundary_conditions" << endl;

      // Pin the pressure at one point, the top right
      unsigned element_index = 0;
      unsigned node_index = 0;
      find_corner_bulk_node(
        Slip_boundary_id, Far_field_boundary_id, element_index, node_index);
      const unsigned pressure_index = 2;
      this->Bulk_mesh_pt->boundary_element_pt(Slip_boundary_id, element_index)
        ->node_pt(node_index)
        ->pin(pressure_index);
      oomph_info << this->Bulk_mesh_pt
                      ->boundary_element_pt(Slip_boundary_id, element_index)
                      ->node_pt(node_index)
                      ->x(0)
                 << ", "
                 << this->Bulk_mesh_pt
                      ->boundary_element_pt(Slip_boundary_id, element_index)
                      ->node_pt(node_index)
                      ->x(1)
                 << " " << endl;
      this->Bulk_mesh_pt->boundary_element_pt(Slip_boundary_id, element_index)
        ->node_pt(node_index)
        ->set_value(pressure_index, 0.0);

      // Fix end points far field boundary condition lagrange_multipliers
      pin_far_field_lagrange_multiplier_end_points();
    }

    void pin_far_field_lagrange_multiplier_end_points()
    {
      const unsigned n_el = Far_field_mesh_pt->nelement();
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
            el_pt->pin_lagrange_multiplier(i_nod, 0);
            el_pt->pin_lagrange_multiplier(i_nod, 1);
          }
        }
      }
    }
  };

  template<class ELEMENT>
  void SectorProblem<ELEMENT>::add_bulk_mesh(const unsigned& n_radial,
                                             const unsigned& n_azimuthal)
  {
    oomph_info << "add_bulk_mesh" << endl;

    // Generate the mesh using the template ELEMENT
    Bulk_mesh_pt = new RefinedSectorTriMesh<ELEMENT>(
      n_radial,
      1.5,
      n_azimuthal,
      0.9,
      135.0 * MathematicalConstants::Pi / 180.0,
      this->time_stepper_pt());

    // Add mesh to problem
    add_sub_mesh(Bulk_mesh_pt);

    // refine_mesh_for_weak_contact_angle_constraint();
  }

  template<class ELEMENT>
  void SectorProblem<ELEMENT>::create_slip_elements()
  {
    oomph_info << "create_slip_elements" << endl;

    unsigned b = Slip_boundary_id;

    // How many bulk elements are adjacent to boundary b?
    unsigned n_element = Bulk_mesh_pt->nboundary_element(b);

    // Loop over the bulk elements adjacent to boundary b?
    for (unsigned e = 0; e < n_element; e++)
    {
      // Get pointer to the bulk element that is adjacent to boundary
      // b
      ELEMENT* bulk_elem_pt =
        dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b, e));

      // What is the index of the face of element e along boundary b
      int face_index = Bulk_mesh_pt->face_index_at_boundary(b, e);
      NavierStokesSlipElement<ELEMENT>* slip_element_pt = 0;

      // Build the corresponding slip element
      slip_element_pt =
        new NavierStokesSlipElement<ELEMENT>(bulk_elem_pt, face_index);

      // Set the pointer to the prescribed slip function
      slip_element_pt->slip_fct_pt() = &Slip_Parameters::prescribed_slip_fct;

      slip_element_pt->wall_velocity_fct_pt() =
        &Slip_Parameters::prescribed_wall_velocity_fct;

      // Add the prescribed-flux element to the surface mesh
      Slip_boundary_mesh_pt->add_element_pt(slip_element_pt);
    }
  }

  template<class ELEMENT>
  void SectorProblem<ELEMENT>::create_no_penetration_elements(
    const unsigned& boundary_id)
  {
    oomph_info << "create_no_penetration_elements" << endl;

    // How many bulk elements are adjacent to boundary b?
    unsigned n_element = Bulk_mesh_pt->nboundary_element(boundary_id);

    // Loop over the bulk elements adjacent to boundary b?
    for (unsigned e = 0; e < n_element; e++)
    {
      // Get pointer to the bulk element that is adjacent to boundary
      // b
      ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Bulk_mesh_pt->boundary_element_pt(boundary_id, e));

      // What is the index of the face of element e along boundary b
      int face_index = Bulk_mesh_pt->face_index_at_boundary(boundary_id, e);

      // Build the corresponding slip element
      ImposeImpenetrabilityElement<ELEMENT>* no_penetration_element_pt =
        new ImposeImpenetrabilityElement<ELEMENT>(
          bulk_elem_pt, face_index, boundary_id);

      // Add the prescribed-flux element to the surface mesh
      No_penetration_boundary_mesh_pt->add_element_pt(
        no_penetration_element_pt);
    }
  }

  template<class ELEMENT>
  void SectorProblem<ELEMENT>::refine_mesh_for_weak_contact_angle_constraint()
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
      set_contact_line_node_pt();

      // Get error estimator
      Bulk_mesh_pt->spatial_error_estimator_pt() =
        new ContactlineErrorEstimator(
          Contact_line_node_pt,
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
  void SectorProblem<ELEMENT>::set_contact_line_node_pt()
  {
    oomph_info << "set_contact_line_node_pt" << endl;

    unsigned element_index;
    unsigned node_index;

    find_corner_bulk_node(
      Slip_boundary_id, Free_surface_boundary_id, element_index, node_index);

    Contact_line_node_pt =
      Bulk_mesh_pt->boundary_element_pt(Slip_boundary_id, element_index)
        ->node_pt(node_index);
  }

  template<class ELEMENT>
  void SectorProblem<ELEMENT>::compute_error_estimate(double& max_err,
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
  void SectorProblem<ELEMENT>::create_far_field_elements()
  {
    oomph_info << "create_far_field_elements" << endl;

    // Loop over the free surface boundary and create the "interface elements
    unsigned b = Far_field_boundary_id;

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
      FarFieldElement<ELEMENT>* el_pt = new FarFieldElement<ELEMENT>(
        bulk_elem_pt, face_index, Far_field_boundary_id);

      // Add it to the mesh
      Far_field_mesh_pt->add_element_pt(el_pt);
    }
  }

  template<class ELEMENT>
  void SectorProblem<ELEMENT>::pin_far_field_elements()
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

    oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
  }

  template<class ELEMENT>
  void SectorProblem<ELEMENT>::find_corner_bulk_element(
    const unsigned& boundary_1_id,
    const unsigned& boundary_2_id,
    unsigned& element_index)
  {
    unsigned dummy_node_index;
    find_corner_bulk_node(
      boundary_1_id, boundary_2_id, element_index, dummy_node_index);
  }

  template<class ELEMENT>
  void SectorProblem<ELEMENT>::find_corner_bulk_node(
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
  void SectorProblem<ELEMENT>::doc_solution()
  {
    unsigned doc_number = Doc_info.number();

    oomph_info << "Doc Number: " << doc_number << endl;

    std::ofstream output_stream;
    char filename[100];
    unsigned npts = 3;

    if (doc_number == 0)
    {
      sprintf(filename, "%s/dofs.txt", Doc_info.directory().c_str());
      output_stream.open(filename);
      describe_dofs(output_stream);
      output_stream.close();
    }

    // Get/output error estimates
    unsigned n_elements = Bulk_mesh_pt->nelement();
    Vector<double> elemental_error(n_elements);
    // We need a dynamic cast, get_element_errors from the Bulk_mesh_pt
    // Dynamic cast is used because get_element_errors require a Mesh* ans
    // not a SolidMesh*
    Mesh* fluid_mesh_pt = dynamic_cast<Mesh*>(Bulk_mesh_pt);
    Z2_error_estimator_pt->get_element_errors(fluid_mesh_pt, elemental_error);
    // Set errors for post-processing and find extrema
    double max_err = 0.0;
    double min_err = 1e8;
    for (unsigned e = 0; e < n_elements; e++)
    {
      dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))
        ->set_error(elemental_error[e]);

      max_err = std::max(max_err, elemental_error[e]);
      min_err = std::min(min_err, elemental_error[e]);
    }

    oomph_info << "Max error is " << max_err << std::endl;
    oomph_info << "Min error is " << min_err << std::endl;

    sprintf(filename,
            "%s/soln%i.dat",
            Doc_info.directory().c_str(),
            Doc_info.number());
    output_stream.open(filename);
    output_stream.precision(15);
    Bulk_mesh_pt->output(output_stream, npts);
    output_stream.close();

    //    sprintf(filename,
    //            "%s/pressure_1_%i.csv",
    //            Doc_info.directory().c_str(),
    //            Doc_info.number());
    //    output_stream.open(filename);
    //    output_stream << "x,y,p," << endl;
    //    Pressure_contribution_mesh_1_pt->output(output_stream, npts);
    //    output_stream.close();
    //
    //    sprintf(filename,
    //            "%s/pressure_2_%i.csv",
    //            Doc_info.directory().c_str(),
    //            Doc_info.number());
    //    output_stream.open(filename);
    //    output_stream << "x,y,p," << endl;
    //    Pressure_contribution_mesh_2_pt->output(output_stream, npts);
    //    output_stream.close();

    sprintf(filename,
            "%s/slip_surface%i.csv",
            Doc_info.directory().c_str(),
            Doc_info.number());
    output_stream.open(filename);
    output_stream << "x,y,l_x,l_y,n_x,n_y,u,v,p" << endl;
    Slip_boundary_mesh_pt->output(output_stream, npts);
    output_stream.close();

    sprintf(filename,
            "%s/no_penetration_surface%i.csv",
            Doc_info.directory().c_str(),
            Doc_info.number());
    output_stream.open(filename);
    output_stream << "x,y,u,v,p,lagrange_multiplier" << endl;
    No_penetration_boundary_mesh_pt->output(output_stream, 3);
    output_stream.close();

    sprintf(filename,
            "%s/far_field%i.dat",
            Doc_info.directory().c_str(),
            Doc_info.number());
    output_stream.open(filename);
    output_stream
      << "x,y,n_x,n_y,u_x,u_y,du_xdr,du_ydr,du_xdtheta,du_ydtheta,l_x,l_y"
      << endl;
    Far_field_mesh_pt->output(output_stream);
    output_stream.close();

    Doc_info.number()++;
  }
} // namespace oomph
#endif
