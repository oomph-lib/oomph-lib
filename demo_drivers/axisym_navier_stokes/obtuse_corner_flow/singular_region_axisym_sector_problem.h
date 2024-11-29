#ifndef SINGULAR_REGION_SECTOR_PROBLEM_HEADER
#define SINGULAR_REGION_SECTOR_PROBLEM_HEADER

#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"
#include "meshes/triangle_mesh.h"

/// Local headers
#include "region_axisym_sector_problem.h"
#include "two_region_refined_sector_tri_mesh.template.h"

#include "eigensolution_functions.h"
#include "pressure_evaluation_elements.h"
#include "point_pressure_evaluation_elements.h"
#include "singular_fluid_traction_elements.h"

namespace oomph
{
  // Problem class
  template<class ELEMENT>
  class SingularRegionAxisymSectorProblem
    : public RegionAxisymSectorProblem<ELEMENT>
  {
  private:
    Node* Contact_line_node_pt;

    Vector<unsigned> Augmented_bulk_element_number;
    Mesh* Singularity_scaling_mesh_pt;
    Mesh* Pressure_contribution_mesh_1_pt;
    Mesh* Pressure_contribution_mesh_2_pt;
    Mesh* Eigensolution_slip_mesh_pt;

    std::function<Vector<double>(const Vector<double>&)>
      Velocity_singular_function;

    std::function<Vector<Vector<double>>(const Vector<double>&)>
      Grad_velocity_singular_function;

    std::function<void(const double&,
                       const Vector<double>&,
                       const Vector<double>&,
                       Vector<double>&)>
      Eigensolution_slip_function;


  public:
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
    SingularRegionAxisymSectorProblem()
      : RegionAxisymSectorProblem<ELEMENT>(), Contact_line_node_pt(0)
    {
      // Re-assign doc info pointer
      this->doc_info_pt()->set_directory("RESLT_axi_fix_region");

      add_singular_sub_meshes();

      this->rebuild_global_mesh();
    }

    void setup()
    {
      // Augment the bulk elements
      setup_and_augment_bulk_elements();

      RegionAxisymSectorProblem<ELEMENT>::setup();

      set_contact_line_node_pt();
      Velocity_singular_function = velocity_singular_function_factory(
        this->my_parameters().sector_angle * MathematicalConstants::Pi / 180.0,
        Contact_line_node_pt);
      Grad_velocity_singular_function = grad_velocity_singular_function_factory(
        this->my_parameters().sector_angle * MathematicalConstants::Pi / 180.0,
        Contact_line_node_pt);
      Eigensolution_slip_function = eigensolution_slip_function_factory(
        this->my_parameters().slip_length, Velocity_singular_function);

      create_singular_elements();

      this->rebuild_global_mesh();
      oomph_info << "Number of unknowns: " << this->assign_eqn_numbers()
                 << std::endl;
    }

    void create_singular_elements()
    {
      // Create the other meshes
      if (!Augmented_bulk_element_number.empty())
      {
        cout << "Make augmented elements" << endl;
        create_singularity_scaling_elements();
        create_pressure_contribution_1_elements();
        create_pressure_contribution_2_elements();

        create_slip_eigen_elements();

        // Setup the mesh interaction between the bulk and singularity meshes
        setup_mesh_interaction();
      }
    }

    void add_singular_sub_meshes()
    {
      Eigensolution_slip_mesh_pt = new Mesh;
      this->add_sub_mesh(Eigensolution_slip_mesh_pt);
      Singularity_scaling_mesh_pt = new Mesh;
      this->add_sub_mesh(Singularity_scaling_mesh_pt);
      Pressure_contribution_mesh_1_pt = new Mesh;
      this->add_sub_mesh(Pressure_contribution_mesh_1_pt);
      Pressure_contribution_mesh_2_pt = new Mesh;
      this->add_sub_mesh(Pressure_contribution_mesh_2_pt);
    }

    void setup_and_augment_bulk_elements()
    {
      set_contact_line_node_pt();

      unsigned n_bulk = this->bulk_mesh_pt()->nelement();
      for (unsigned e = 0; e < n_bulk; e++)
      {
        // Upcast from GeneralisedElement to the present element
        ELEMENT* el_pt =
          dynamic_cast<ELEMENT*>(this->bulk_mesh_pt()->element_pt(e));

        // Make augmented elements
        // Check distance from
        // s centre is centre of mass of a uniform triangle, so (1/3,1/3) for
        // Triangle[(1,0),(0,1),(0,0)]
        Vector<double> s_centre(2, 1.0 / 3.0);
        Vector<double> element_centre_x(2, 0.0);
        el_pt->get_x(s_centre, element_centre_x);
        double dist = 0;
        for (unsigned i = 0; i < 2; i++)
        {
          dist += pow(element_centre_x[i] - Contact_line_node_pt->x(i), 2.0);
        }
        dist = pow(dist, 0.5);


        // If the distance to the corner is within the "inner" region, ...
        if (el_pt->get_region_id() ==
            TwoRegionRefinedSectorTriMesh<ELEMENT>::Inner_region_id)
        {
          // ... augment element
          el_pt->augment();

          el_pt->add_additional_terms();

          el_pt->swap_unknowns();

          Augmented_bulk_element_number.push_back(e);
          // for (unsigned n = 0; n < 6; n++)
          //{
          //   for (unsigned d = 0; d < 2; d++)
          //   {
          //     el_pt->pin_total_velocity_eqn(n, d);
          //   }
          // }
        }
      }
      oomph_info << Augmented_bulk_element_number.size()
                 << " augmented elements" << std::endl;
    }

    void set_contact_line_node_pt()
    {
      oomph_info << "set_contact_line_node_pt" << endl;

      unsigned element_index;
      unsigned node_index;

      find_corner_bulk_node(Inner_slip_boundary_id,
                            Inner_free_surface_boundary_id,
                            element_index,
                            node_index);

      Contact_line_node_pt =
        this->bulk_mesh_pt()
          ->boundary_element_pt(Inner_slip_boundary_id, element_index)
          ->node_pt(node_index);
      oomph_info << Contact_line_node_pt->x(0) << ", "
                 << Contact_line_node_pt->x(1) << std::endl;
    }

    void doc_solution()
    {
      char filename[100];
      sprintf(filename,
              "%s/eigenslip_surface%i.csv",
              this->doc_info_pt()->directory().c_str(),
              this->doc_info_pt()->number());
      std::ofstream output_stream;
      output_stream.open(filename);
      output_stream << "x,y,n_x,n_y,t_x,t_y" << endl;
      const unsigned npts = 3;
      Eigensolution_slip_mesh_pt->output(output_stream, npts);
      output_stream.close();

      sprintf(filename,
              "%s/scaling%i.csv",
              this->doc_info_pt()->directory().c_str(),
              this->doc_info_pt()->number());
      output_stream.open(filename);
      output_stream << "scaling" << endl;
      dynamic_cast<SingularNavierStokesSolutionElement<ELEMENT>*>(
        Singularity_scaling_mesh_pt->element_pt(0))
        ->output(output_stream);
      output_stream.close();

      RegionAxisymSectorProblem<ELEMENT>::doc_solution();
    }

  private:
    void create_slip_eigen_elements();
    void create_singularity_scaling_elements();
    void create_pressure_contribution_1_elements();
    void create_pressure_contribution_2_elements();
    void find_corner_bulk_element_and_face_index(const unsigned& boundary_1_id,
                                                 const unsigned& boundary_2_id,
                                                 ELEMENT*& element_pt,
                                                 int& face_index);

    void find_corner_bulk_element(const unsigned& boundary_1_id,
                                  const unsigned& boundary_2_id,
                                  unsigned& element_index);
    void find_corner_bulk_node(const unsigned& boundary_1_id,
                               const unsigned& boundary_2_id,
                               unsigned& element_index,
                               unsigned& node_index);

    void setup_mesh_interaction();

    void fix_c(const double& value)
    {
      SingularNavierStokesSolutionElement<ELEMENT>* el_pt =
        dynamic_cast<SingularNavierStokesSolutionElement<ELEMENT>*>(
          Singularity_scaling_mesh_pt->element_pt(0));

      el_pt->pin_c();
      el_pt->set_c(value);
    }
  };

  template<class ELEMENT>
  void SingularRegionAxisymSectorProblem<ELEMENT>::create_slip_eigen_elements()
  {
    oomph_info << "create_slip_eigen_elements" << endl;

    // Loop over the free surface boundary and create the "interface elements
    unsigned b = Inner_slip_boundary_id;

    // How many bulk fluid elements are adjacent to boundary b?
    unsigned n_element = this->bulk_mesh_pt()->nboundary_element(b);

    // Loop over the bulk fluid elements adjacent to boundary b?
    for (unsigned e = 0; e < n_element; e++)
    {
      // Get pointer to the bulk fluid element that is
      // adjacent to boundary b
      ELEMENT* bulk_elem_pt =
        dynamic_cast<ELEMENT*>(this->bulk_mesh_pt()->boundary_element_pt(b, e));

      if (bulk_elem_pt->is_augmented())
      {
        // Find the index of the face of element e along boundary b
        int face_index = this->bulk_mesh_pt()->face_index_at_boundary(b, e);

        // Create new element
        SingularNavierStokesTractionElement<ELEMENT>* el_pt =
          new SingularNavierStokesTractionElement<ELEMENT>(
            bulk_elem_pt,
            face_index,
            Singularity_scaling_mesh_pt->element_pt(0)->internal_data_pt(0));

        el_pt->set_traction_fct(Eigensolution_slip_function);

        // Add it to the mesh
        Eigensolution_slip_mesh_pt->add_element_pt(el_pt);
      }
    }
  }

  template<class ELEMENT>
  void SingularRegionAxisymSectorProblem<
    ELEMENT>::create_singularity_scaling_elements()
  {
    oomph_info << "create_singularity_scaling_elements" << endl;
    SingularNavierStokesSolutionElement<ELEMENT>* el_pt =
      new SingularNavierStokesSolutionElement<ELEMENT>;

    // Set the pointer to the velocity singular function for this
    // element, defined in parameters namespace
    el_pt->velocity_singular_fct() = Velocity_singular_function;

    // Set the pointer to the gradient of the velocity singular
    // function for this element, defined in parameters namespace
    el_pt->grad_velocity_singular_fct() = Grad_velocity_singular_function;

    // Set the pointer to the first pressure singular function for this
    // element, defined in parameters namespace
    el_pt->pressure_singular_fct_pt() = &pressure_singular_fct;

    // The singular function satisfies the Stokes equation
    el_pt->singular_function_satisfies_stokes_equation() = false;

    // el_pt->pin_c();
    el_pt->set_c(0.0);

    // Add element to the mesh
    Singularity_scaling_mesh_pt->add_element_pt(el_pt);
  }

  template<class ELEMENT>
  void SingularRegionAxisymSectorProblem<ELEMENT>::find_corner_bulk_element(
    const unsigned& boundary_1_id,
    const unsigned& boundary_2_id,
    unsigned& element_index)
  {
    unsigned dummy_node_index;
    find_corner_bulk_node(
      boundary_1_id, boundary_2_id, element_index, dummy_node_index);
  }

  template<class ELEMENT>
  void SingularRegionAxisymSectorProblem<ELEMENT>::find_corner_bulk_node(
    const unsigned& boundary_1_id,
    const unsigned& boundary_2_id,
    unsigned& element_index,
    unsigned& node_index)
  {
    unsigned n_boundary_element =
      this->bulk_mesh_pt()->nboundary_element(boundary_1_id);
    for (unsigned e = 0; e < n_boundary_element; e++)
    {
      // Locally cache the element pointer
      FiniteElement* bulk_el_pt =
        this->bulk_mesh_pt()->boundary_element_pt(boundary_1_id, e);

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
  void SingularRegionAxisymSectorProblem<ELEMENT>::setup_mesh_interaction()
  {
    oomph_info << "setup_mesh_interaction" << endl;

    SingularNavierStokesSolutionElement<ELEMENT>* singular_el_pt =
      dynamic_cast<SingularNavierStokesSolutionElement<ELEMENT>*>(
        Singularity_scaling_mesh_pt->element_pt(0));

    // Loop over the augmented bulk elements
    unsigned n_aug_bulk = Augmented_bulk_element_number.size();
    for (unsigned e = 0; e < n_aug_bulk; e++)
    {
      // Augment elements
      // Upcast from GeneralisedElement to the present element
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
        this->bulk_mesh_pt()->element_pt(Augmented_bulk_element_number[e]));

      // Set the pointer to the element that determines the amplitude
      // of the singular fct
      el_pt->add_c_equation_element_pt(singular_el_pt);
    }
  }

  template<class ELEMENT>
  void SingularRegionAxisymSectorProblem<
    ELEMENT>::create_pressure_contribution_1_elements()
  {
    oomph_info << "create_pressure_contribution_1_elements" << std::endl;

    ELEMENT* element_pt = 0;
    int face_index = 0;
    find_corner_bulk_element_and_face_index(Inner_slip_boundary_id,
                                            Inner_free_surface_boundary_id,
                                            element_pt,
                                            face_index);

    // PressureEvaluationElement<ELEMENT>* el_pt =
    //   new PressureEvaluationElement<ELEMENT>(
    //     element_pt, face_index, dynamic_cast<Node*>(Contact_line_node_pt));
    // el_pt->set_boundary_number_in_bulk_mesh(Slip_boundary_id);

    Node* node_pt = 0;
    for (unsigned n = 0; n < 3; n++)
    {
      node_pt = element_pt->node_pt(n);
      if (node_pt->is_on_boundary(Inner_slip_boundary_id) &&
          !node_pt->is_on_boundary(Inner_free_surface_boundary_id))
      {
        break;
      }
    }

    oomph_info << node_pt->x(0) << ", " << node_pt->x(1) << std::endl;

    PointPressureEvaluationElement* el_pt =
      new PointPressureEvaluationElement(node_pt);

    el_pt->set_pressure_data_pt(
      Singularity_scaling_mesh_pt->element_pt(0)->internal_data_pt(0));

    Pressure_contribution_mesh_1_pt->add_element_pt(el_pt);
  }

  template<class ELEMENT>
  void SingularRegionAxisymSectorProblem<
    ELEMENT>::create_pressure_contribution_2_elements()
  {
    oomph_info << "create_pressure_contribution_2_elements" << std::endl;

    ELEMENT* element_pt = 0;
    int face_index = 0;
    find_corner_bulk_element_and_face_index(Inner_free_surface_boundary_id,
                                            Inner_slip_boundary_id,
                                            element_pt,
                                            face_index);

    // PressureEvaluationElement<ELEMENT>* el_pt =
    //   new PressureEvaluationElement<ELEMENT>(
    //     element_pt, face_index, dynamic_cast<Node*>(Contact_line_node_pt));
    // el_pt->set_boundary_number_in_bulk_mesh(Free_surface_boundary_id);

    Node* node_pt = 0;
    for (unsigned n = 0; n < 3; n++)
    {
      node_pt = element_pt->node_pt(n);
      if (!node_pt->is_on_boundary(Inner_slip_boundary_id) &&
          node_pt->is_on_boundary(Inner_free_surface_boundary_id))
      {
        break;
      }
    }

    oomph_info << node_pt->x(0) << ", " << node_pt->x(1) << std::endl;

    PointPressureEvaluationElement* el_pt =
      new PointPressureEvaluationElement(node_pt);

    el_pt->set_pressure_data_pt(
      Singularity_scaling_mesh_pt->element_pt(0)->internal_data_pt(0));
    el_pt->set_subtract_from_residuals();

    Pressure_contribution_mesh_2_pt->add_element_pt(el_pt);
  }

  template<class ELEMENT>
  void SingularRegionAxisymSectorProblem<ELEMENT>::
    find_corner_bulk_element_and_face_index(const unsigned& boundary_1_id,
                                            const unsigned& boundary_2_id,
                                            ELEMENT*& element_pt,
                                            int& face_index)
  {
    unsigned n_boundary_element =
      this->bulk_mesh_pt()->nboundary_element(boundary_1_id);
    for (unsigned e = 0; e < n_boundary_element; e++)
    {
      // Locally cache the element pointer
      FiniteElement* bulk_el_pt =
        this->bulk_mesh_pt()->boundary_element_pt(boundary_1_id, e);

      // Read out number of nodes in the element
      unsigned n_node = bulk_el_pt->nnode();
      for (unsigned i_node = 0; i_node < n_node; i_node++)
      {
        // If the node is on the free surface boundary as well then ...
        if (bulk_el_pt->node_pt(i_node)->is_on_boundary(boundary_2_id) &&
            bulk_el_pt->node_pt(i_node)->is_on_boundary(boundary_1_id))
        {
          // set the output arguments,
          element_pt = dynamic_cast<ELEMENT*>(bulk_el_pt);
          face_index =
            this->bulk_mesh_pt()->face_index_at_boundary(boundary_1_id, e);

          // Return to exit both loops and end function
          return;
        }
      }
    }
    // If not found, issue warning and return anyway
    oomph_info << "Warning: No corner node found!" << std::endl;
  }
} // namespace oomph
#endif
