#ifndef REGION_SECTOR_PROBLEM_HEADER
#define REGION_SECTOR_PROBLEM_HEADER

#include <cmath>
#include <array>

#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "meshes/triangle_mesh.h"

/// Local headers
#include "base_problem.h"
#include "parameter_functions.h"
#include "parameter_struct.h"
#include "two_region_refined_sector_tri_mesh.template.h"
#include "two_region_refined_sector_tri_mesh.template.cc"
#include "info_elements.h"

namespace oomph
{
  // Problem class
  template<class ELEMENT>
  class RegionAxisymSectorProblem : public BaseProblem
  {
  private:
    /// Private variables
    TwoRegionRefinedSectorTriMesh<ELEMENT>* Bulk_mesh_pt;

    Mesh* No_penetration_boundary_mesh1_pt;
    Mesh* No_penetration_boundary_mesh2_pt;
    Mesh* Slip_boundary_mesh_pt;
    unsigned AxisymFluxComputeMeshIndex;
    unsigned AxisymFluxOutputMeshIndex;

    Z2ErrorEstimator* Z2_error_estimator_pt;

    std::function<void(const double&,
                       const Vector<double>&,
                       const Vector<double>&,
                       Vector<double>&)>
      Slip_function;

  public:
    /// Copy of two region refined sector boundary ids
    enum
    {
      Slip_boundary_id,
      Far_field_boundary_id,
      Free_surface_boundary_id,
      Inner_free_surface_boundary_id,
      Inner_slip_boundary_id,
      Inner_boundary_id,
    };

    // Constructor
    RegionAxisymSectorProblem() : Z2_error_estimator_pt(new Z2ErrorEstimator)
    {
      // Create and add the timestepper
      add_time_stepper_pt(new BDF<2>);

      // Assign doc info pointer
      doc_info_pt()->set_directory("RESLT_axi_no_fix_region");

      /// Create parameters from parameters file.
      Slip_function = slip_function_factory(parameters().slip_length);

      // Create an empty mesh
      add_bulk_mesh();
      add_non_adaptive_sub_meshes();
      build_global_mesh();
    }

    RegionAxisymSectorProblem(Params& params)
      : BaseProblem(params), Z2_error_estimator_pt(new Z2ErrorEstimator)
    {
      // Create and add the timestepper
      add_time_stepper_pt(new BDF<2>);

      // Assign doc info pointer
      doc_info_pt()->set_directory("RESLT_axi_no_fix_region");

      /// Create parameters from parameters file.
      Slip_function = slip_function_factory(parameters().slip_length);

      // Create an empty mesh
      add_bulk_mesh();
      add_non_adaptive_sub_meshes();
      build_global_mesh();
    }

    virtual void setup()
    {
      this->setup_bulk_elements();
      create_nonrefineable_elements();
      set_boundary_conditions();

      this->rebuild_global_mesh();
      oomph_info << "Number of unknowns: " << this->assign_eqn_numbers()
                 << std::endl;
    }

    void create_nonrefineable_elements()
    {
      create_slip_elements();
      create_no_penetration1_elements();
      create_no_penetration2_elements();

      // create_axisym_flux_elements();
    }

    void delete_nonrefineable_elements()
    {
      delete_elements(Slip_boundary_mesh_pt);
      delete_elements(No_penetration_boundary_mesh1_pt);
      delete_elements(No_penetration_boundary_mesh2_pt);
    }

    void add_non_adaptive_sub_meshes()
    {
      No_penetration_boundary_mesh1_pt = new Mesh;
      add_sub_mesh(No_penetration_boundary_mesh1_pt);
      No_penetration_boundary_mesh2_pt = new Mesh;
      add_sub_mesh(No_penetration_boundary_mesh2_pt);
      Slip_boundary_mesh_pt = new Mesh;
      add_sub_mesh(Slip_boundary_mesh_pt);
    }

    void setup_bulk_elements()
    {
      unsigned n_bulk = Bulk_mesh_pt->nelement();
      for (unsigned e = 0; e < n_bulk; e++)
      {
        // Upcast from GeneralisedElement to the present element
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

        // Set the Reynolds number
        el_pt->re_pt() = &parameters().reynolds_number;

        // Set the Reynolds Strouhal number
        el_pt->re_st_pt() = &parameters().strouhal_reynolds_number;
      }
    }

    // Destructor
    ~RegionAxisymSectorProblem()
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

    TwoRegionRefinedSectorTriMesh<ELEMENT>*& bulk_mesh_pt()
    {
      return Bulk_mesh_pt;
    }


    void add_bulk_mesh();

  public:
    void pin_no_penetration_conditions()
    {
      unsigned n_el = No_penetration_boundary_mesh1_pt->nelement();
      for (unsigned i_el = 0; i_el < n_el; i_el++)
      {
        ImposeImpenetrabilityElement<ELEMENT>* el_pt =
          dynamic_cast<ImposeImpenetrabilityElement<ELEMENT>*>(
            No_penetration_boundary_mesh1_pt->element_pt(i_el));
        const unsigned n_nod = el_pt->nnode();
        for (unsigned i_nod = 0; i_nod < n_nod; i_nod++)
        {
          el_pt->pin_lagrange_multiplier(i_nod);
        }
      }

      n_el = No_penetration_boundary_mesh2_pt->nelement();
      for (unsigned i_el = 0; i_el < n_el; i_el++)
      {
        ImposeImpenetrabilityElement<ELEMENT>* el_pt =
          dynamic_cast<ImposeImpenetrabilityElement<ELEMENT>*>(
            No_penetration_boundary_mesh1_pt->element_pt(i_el));
        const unsigned n_nod = el_pt->nnode();
        for (unsigned i_nod = 0; i_nod < n_nod; i_nod++)
        {
          el_pt->pin_lagrange_multiplier(i_nod);
        }
      }
    }

    void set_boundary_conditions()
    {
      oomph_info << "set_boundary_conditions" << std::endl;

      // Find a node in the interior of the domain and pin it.
      bool has_node_been_found = false;
      unsigned element_index =
        std::floor(Bulk_mesh_pt->nboundary_element(Far_field_boundary_id) / 2);
      unsigned node_index = 0;
      while (!has_node_been_found)
      {
        // Select a candidate
        Node* node_pt =
          this->Bulk_mesh_pt
            ->boundary_element_pt(Far_field_boundary_id, element_index)
            ->node_pt(node_index);

        // Test if it is internal
        if (!node_pt->is_on_boundary(Far_field_boundary_id))
        {
          // Pin the pressure at one point
          const unsigned pressure_index = 3;
          node_pt->pin(pressure_index);
          node_pt->set_value(pressure_index, 0.0);
          oomph_info << "Pressure pinned at: " << node_pt->x(0) << ", "
                     << node_pt->x(1) << " " << std::endl;
          has_node_been_found = true;
        }
        else
        {
          // Seek a new candidate
          if (node_index < 6)
          {
            node_index++;
          }
          else
          {
            node_index = 0;
            element_index++;
          }
        }
      }

      // Pin the lagrange multiplier for the no penetration condition on the
      // vertical surface at the centre
      // pin_lagrange_multiplier_end_point(Free_surface_boundary_id,
      //                                  Slip_boundary_id);
    }

  private:
    void create_slip_elements();
    void create_no_penetration1_elements();
    void create_no_penetration2_elements();
    void create_axisym_flux_elements()
    {
      std::cout << "create_axisym_flux_elements" << std::endl;
      AxisymFluxOutputMeshIndex = add_sub_mesh(new Mesh);
      InfoElement* info_element_pt = new InfoElement;
      Data* net_flux_data_pt = info_element_pt->new_internal_data_pt();
      mesh_pt(AxisymFluxOutputMeshIndex)->add_element_pt(info_element_pt);

      AxisymFluxComputeMeshIndex = add_sub_mesh(new Mesh);
      const unsigned n_element =
        Bulk_mesh_pt->nboundary_element(Inner_boundary_id);
      std::cout << "n_element: " << n_element << std::endl;
      for (unsigned e = 0; e < n_element; e++)
      {
        ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
          Bulk_mesh_pt->boundary_element_pt(Inner_boundary_id, e));
        int face_index =
          Bulk_mesh_pt->face_index_at_boundary(Inner_boundary_id, e);
        AxisymmetricFluidFluxElement<ELEMENT>* flux_element_pt =
          new AxisymmetricFluidFluxElement<ELEMENT>(
            bulk_elem_pt, face_index, net_flux_data_pt);
        mesh_pt(AxisymFluxComputeMeshIndex)->add_element_pt(flux_element_pt);
        std::cout << "Added flux element" << std::endl;
      }
    }

    void find_corner_bulk_node(const unsigned& boundary_1_id,
                               const unsigned& boundary_2_id,
                               unsigned& element_index,
                               unsigned& node_index);

    void pin_lagrange_multiplier_end_point(const unsigned boundary_1_id,
                                           const unsigned boundary_2_id)
    {
      const unsigned n_el = No_penetration_boundary_mesh1_pt->nelement();
      for (unsigned i_el = 0; i_el < n_el; i_el++)
      {
        ImposeImpenetrabilityElement<ELEMENT>* el_pt =
          dynamic_cast<ImposeImpenetrabilityElement<ELEMENT>*>(
            No_penetration_boundary_mesh1_pt->element_pt(i_el));
        const unsigned n_nod = el_pt->nnode();
        for (unsigned i_nod = 0; i_nod < n_nod; i_nod++)
        {
          // Get boundary node
          const Node* const node_pt = el_pt->node_pt(i_nod);
          // If node is on either of the other boundaries
          if (node_pt->is_on_boundary(boundary_1_id) &&
              node_pt->is_on_boundary(boundary_2_id))
          {
            // then fix the lagrange multiplier to zero
            oomph_info << "Fix lagrange_multiplier" << std::endl;
            el_pt->pin_lagrange_multiplier(i_nod);
          }
        }
      }
    }
  };

  template<class ELEMENT>
  void RegionAxisymSectorProblem<ELEMENT>::add_bulk_mesh()
  {
    oomph_info << "add_bulk_mesh" << std::endl;

    // Generate the mesh using the template ELEMENT
    Bulk_mesh_pt = new TwoRegionRefinedSectorTriMesh<ELEMENT>(
      parameters().n_radial,
      4,
      parameters().geometric_base,
      parameters().n_azimuthal,
      parameters().sector_radius,
      parameters().sector_angle * MathematicalConstants::Pi / 180.0,
      this->time_stepper_pt());

    // Add mesh to problem
    add_sub_mesh(Bulk_mesh_pt);
  }

  template<class ELEMENT>
  void RegionAxisymSectorProblem<ELEMENT>::create_slip_elements()
  {
    oomph_info << "create_slip_elements" << std::endl;

    std::array boundary_ids = {Inner_slip_boundary_id, Slip_boundary_id};
    for (unsigned b : boundary_ids)
    {
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
        AxisymmetricNavierStokesSlipElement<ELEMENT>* slip_element_pt = 0;

        // Build the corresponding slip element
        slip_element_pt = new AxisymmetricNavierStokesSlipElement<ELEMENT>(
          bulk_elem_pt, face_index);

        // Set the pointer to the prescribed slip function
        slip_element_pt->set_slip_function(Slip_function);

        slip_element_pt->set_wall_velocity_function(wall_velocity_function);

        // Add the prescribed-flux element to the surface mesh
        Slip_boundary_mesh_pt->add_element_pt(slip_element_pt);
      }
    }
  }

  template<class ELEMENT>
  void RegionAxisymSectorProblem<ELEMENT>::create_no_penetration1_elements()
  {
    oomph_info << "create_no_penetration1_elements" << std::endl;

    std::array boundary_ids = {Inner_slip_boundary_id, Slip_boundary_id};
    for (unsigned b : boundary_ids)
    {
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

        // Build the corresponding slip element
        ImposeImpenetrabilityElement<ELEMENT>* no_penetration_element_pt =
          new ImposeImpenetrabilityElement<ELEMENT>(
            bulk_elem_pt, face_index, Slip_boundary_id);

        // Add the prescribed-flux element to the surface mesh
        No_penetration_boundary_mesh1_pt->add_element_pt(
          no_penetration_element_pt);
      }
    }
  }

  template<class ELEMENT>
  void RegionAxisymSectorProblem<ELEMENT>::create_no_penetration2_elements()
  {
    oomph_info << "create_no_penetration2_elements" << std::endl;

    std::array boundary_ids = {Inner_free_surface_boundary_id,
                               Free_surface_boundary_id};
    for (unsigned b : boundary_ids)
    {
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

        // Build the corresponding slip element
        ImposeImpenetrabilityElement<ELEMENT>* no_penetration_element_pt =
          new ImposeImpenetrabilityElement<ELEMENT>(
            bulk_elem_pt, face_index, Free_surface_boundary_id);

        // Add the prescribed-flux element to the surface mesh
        No_penetration_boundary_mesh2_pt->add_element_pt(
          no_penetration_element_pt);
      }
    }
  }

  template<class ELEMENT>
  void RegionAxisymSectorProblem<ELEMENT>::find_corner_bulk_node(
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
    oomph_info << "Warning: No corner node found!" << std::endl;

    return;
  }

  template<class ELEMENT>
  void RegionAxisymSectorProblem<ELEMENT>::doc_solution()
  {
    unsigned doc_number = doc_info_pt()->number();

    oomph_info << "Doc Number: " << doc_number << std::endl;

    std::ofstream output_stream;
    char filename[100];
    unsigned npts = 3;

    double max_err = 0.0;
    double min_err = 1e8;
    /*
    // Get/output error estimates
    unsigned n_elements = Bulk_mesh_pt->nelement();
    Vector<double> elemental_error(n_elements);
    // We need a dynamic cast, get_element_errors from the Bulk_mesh_pt
    // Dynamic cast is used because get_element_errors require a Mesh* ans
    // not a SolidMesh*
    Mesh* fluid_mesh_pt = dynamic_cast<Mesh*>(Bulk_mesh_pt);
    Z2_error_estimator_pt->get_element_errors(fluid_mesh_pt, elemental_error);
    // Set errors for post-processing and find extrema
    for (unsigned e = 0; e < n_elements; e++)
    {
      dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))
        ->set_error(elemental_error[e]);

      max_err = std::max(max_err, elemental_error[e]);
      min_err = std::min(min_err, elemental_error[e]);
    }
    */

    oomph_info << "Max error is " << max_err << std::endl;
    oomph_info << "Min error is " << min_err << std::endl;

    sprintf(filename,
            "%s/soln%i.dat",
            doc_info_pt()->directory().c_str(),
            doc_info_pt()->number());
    output_stream.open(filename);
    output_stream.precision(15);
    Bulk_mesh_pt->output(output_stream, npts);
    output_stream.close();

    //    sprintf(filename,
    //            "%s/pressure_1_%i.csv",
    //            doc_info_pt()->directory().c_str(),
    //            doc_info_pt()->number());
    //    output_stream.open(filename);
    //    output_stream << "x,y,p," << std::endl;
    //    Pressure_contribution_mesh_1_pt->output(output_stream, npts);
    //    output_stream.close();
    //
    //    sprintf(filename,
    //            "%s/pressure_2_%i.csv",
    //            doc_info_pt()->directory().c_str(),
    //            doc_info_pt()->number());
    //    output_stream.open(filename);
    //    output_stream << "x,y,p," << std::endl;
    //    Pressure_contribution_mesh_2_pt->output(output_stream, npts);
    //    output_stream.close();

    sprintf(filename,
            "%s/slip_surface%i.csv",
            doc_info_pt()->directory().c_str(),
            doc_info_pt()->number());
    output_stream.open(filename);
    output_stream << "x y l_x l_y n_x n_y u_w v_w z_w u v w p" << std::endl;
    Slip_boundary_mesh_pt->output(output_stream, npts);
    output_stream.close();

    sprintf(filename,
            "%s/no_penetration_surface%i.csv",
            doc_info_pt()->directory().c_str(),
            doc_info_pt()->number());
    output_stream.open(filename);
    output_stream << "x y u v p lagrange_multiplier nx ny " << std::endl;
    No_penetration_boundary_mesh1_pt->output(output_stream, 3);
    output_stream.close();

    sprintf(filename,
            "%s/no_penetration_surface2_%i.csv",
            doc_info_pt()->directory().c_str(),
            doc_info_pt()->number());
    output_stream.open(filename);
    output_stream << "x y u v p lagrange_multiplier nx ny " << std::endl;
    No_penetration_boundary_mesh2_pt->output(output_stream, 3);
    output_stream.close();

    // sprintf(filename,
    //         "%s/axisym_flux%i.csv",
    //         doc_info_pt()->directory().c_str(),
    //         doc_info_pt()->number());
    // output_stream.open(filename);
    // output_stream << "x y u v" << std::endl;
    // this->mesh_pt(AxisymFluxComputeMeshIndex)->output(output_stream);
    // output_stream.close();

    // sprintf(filename,
    //         "%s/total_flux%i.csv",
    //         doc_info_pt()->directory().c_str(),
    //         doc_info_pt()->number());
    // output_stream.open(filename);
    // output_stream << "Q" << std::endl;
    // dynamic_cast<InfoElement*>(
    //   this->mesh_pt(AxisymFluxOutputMeshIndex)->element_pt(0))
    //   ->output(output_stream);
    // output_stream.close();

    doc_info_pt()->number()++;
  }
} // namespace oomph
#endif
