#ifndef PERTURBED_LINEAR_STABILITY_CAP_PROBLEM_HEADER
#define PERTURBED_LINEAR_STABILITY_CAP_PROBLEM_HEADER

//#include "axisym_linear_stability_cap_problem.h"
#include "overlaying_elastic_linearised_axisym_fluid_slip_elements.h"
#include "overlaying_linearised_elastic_axisym_fluid_interface_element.h"
#include "parameters.h"
#include "volume_constraint_elements_with_output.h"
#include "decomposed_integral_elements.h"
#include "decomposed_flux_elements.h"
#include "linearised_contact_angle_elements.h"
#include "my_eigenproblem.h"
#include "my_triangle_mesh.h"
#include "symmetry_velocity_condition_element.h"
#include "utility_functions.h"
#include "net_flux_elements.h"

namespace oomph
{
  //===========start_of_perturbed_linear_stability_cap_problem_class============
  // A class that solves the linearised Axisymmetric Navier--Stokes equations
  // to compute the stability of a dynamic interface in a cylindrical container
  // with imposed contact angle at the boundary.
  //============================================================================
  template<class BASE_ELEMENT, class PERTURBED_ELEMENT, class TIMESTEPPER>
  class PerturbedLinearStabilityCapProblem : public MyEigenproblem
  {
  private:
    // Typedef all the elements so that they are in one place and it is
    // clearer.
    typedef OverlayingLinearisedElasticAxisymmetricFluidInterfaceElement<
      BASE_ELEMENT,
      PERTURBED_ELEMENT>
      FREE_SURFACE_ELEMENT;
    typedef OverlayingLinearisedElasticAxisymmetricFluidSlipElement<
      BASE_ELEMENT,
      PERTURBED_ELEMENT>
      SLIP_ELEMENT;
    typedef IntegralElement VOLUME_ELEMENT;
    typedef NetFluxElement NET_FLUX_ELEMENT;
    typedef DataTradingElement VOLUME_CONSTRAINT_ELEMENT;
    typedef DecomposedIntegralElements<PERTURBED_ELEMENT>
      VOLUME_COMPUTATION_ELEMENT;
    typedef DecomposedFluxElements<PERTURBED_ELEMENT> FLUX_COMPUTATION_ELEMENT;
    typedef SymmetryVelocityConditionElement SYMMETRY_CONDITION_ELEMENT;
    typedef LinearisedContactAngleElement<FREE_SURFACE_ELEMENT>
      CONTACT_LINE_ELEMENT;

    MyTriangleMesh<PERTURBED_ELEMENT>* Fluid_mesh_pt;
    Mesh* Free_surface_mesh_pt;
    Mesh* Slip_mesh_pt;
    Mesh* Integral_mesh_pt;
    Mesh* Volume_mesh_pt;
    Mesh* Flux_mesh_pt;
    Mesh* Net_flux_mesh_pt;
    Mesh* Centre_mesh_pt;
    Mesh* Contact_line_mesh_pt;

    Mesh* External_base_mesh_pt;
    Mesh* External_free_surface_mesh_pt;
    Mesh* External_slip_surface_mesh_pt;

    Data* Tilde_external_pressure_data_pt;
    Data* Volume_data_pt;
    Data* Flux_lagrange_multiplier_data_pt;

    int Azimuthal_mode_number;

    double Volume;
    double Bo;
    double Ca;
    double Contact_angle;
    double Re;
    double St;
    double ReSt;
    double ReInvFr;
    ConstitutiveLaw* Constitutive_law_pt;
    IsotropicElasticityTensor* ElasticityTensor_pt;

    enum Lagrange_id
    {
      Kinematic,
      Centre,
    };

    enum
    {
      Upper_boundary_id,
      Outer_boundary_with_slip_id,
      Free_surface_boundary_id,
      Inner_boundary_id,
    };

    // Radial, Vertical, Azimuthal
    enum
    {
      rc_index,
      rs_index,
      zc_index,
      zs_index,
      uc_index,
      us_index,
      wc_index,
      ws_index,
      vc_index,
      vs_index,
      pc_index,
      ps_index,
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

  public:
    // Constructor
    // Uses the external base and surface mesh of the base state to create a
    // new problem solving the linearised system on the same mesh geometry.
    // This was easier than using the multi domain elements due to the
    // complexities on the free surface boundary and the contact line.
    // The azimuthal_mode_number is passed in as a parameter.
    PerturbedLinearStabilityCapProblem(Mesh* external_base_mesh_pt,
                                       Mesh* external_free_surface_mesh_pt,
                                       Mesh* external_slip_surface_mesh_pt,
                                       const int& azimuthal_mode_number)
      : Free_surface_mesh_pt(0),
        Slip_mesh_pt(0),
        Integral_mesh_pt(0),
        Volume_mesh_pt(0),
        Flux_mesh_pt(0),
        Net_flux_mesh_pt(0),
        Centre_mesh_pt(0),
        Contact_line_mesh_pt(0),
        External_base_mesh_pt(external_base_mesh_pt),
        External_free_surface_mesh_pt(external_free_surface_mesh_pt),
        External_slip_surface_mesh_pt(external_slip_surface_mesh_pt),
        Volume_data_pt(0),
        Flux_lagrange_multiplier_data_pt(0),
        Azimuthal_mode_number(azimuthal_mode_number),
        Volume((Global_Physical_Parameters::Initial_fluid_height) / 2.0),
        Bo(1.0),
        Ca(1.0),
        Contact_angle(1.0),
        Re(0.0),
        St(1.0),
        ReSt(0.0),
        ReInvFr(1.0),
        problem_type(Bulk_only_problem),
        Is_steady(true)
    {
      // Problem::Always_take_one_newton_step = true;

      // Add the time stepper
      add_time_stepper_pt(new TIMESTEPPER);

      // Set up the other variables.
      this->Constitutive_law_pt =
        new GeneralisedHookean(&Mesh_Control_Parameters::Nu);
      ElasticityTensor_pt = new IsotropicElasticityTensor(0.1);
      Tilde_external_pressure_data_pt = new Data(1);
      Tilde_external_pressure_data_pt->set_value(0, 0.0);
      this->doc_info().number() = 0;

      // Create the fluid mesh
      Fluid_mesh_pt =
        new MyTriangleMesh<PERTURBED_ELEMENT>(this->time_stepper_pt());
      Fluid_mesh_pt->deep_copy(External_base_mesh_pt);

      // Setup up the bulk mesh elements with the parameters etc.
      setup_bulk_elements();

      // Create the free surface elements
      Free_surface_mesh_pt = new Mesh;
      create_free_surface_elements();

      // Create the slip surface elements
      Slip_mesh_pt = new Mesh;
      create_slip_elements();

      // Create the integral elements if we are solving the axisymmetric
      // problem
      if (this->Azimuthal_mode_number == 0)
      {
        Integral_mesh_pt = new Mesh;
        Volume_mesh_pt = new Mesh;
        create_integral_elements();
      }

      // Create the flux elements if we are solving the axisymmetric
      // problem
      if (this->Azimuthal_mode_number == 0)
      {
        Flux_mesh_pt = new Mesh;
        Net_flux_mesh_pt = new Mesh;
        create_flux_elements();
      }

      // Create the velocity symmetry condition elements if we are solving the
      // mode 1 problem
      if (this->Azimuthal_mode_number == 1)
      {
        Centre_mesh_pt = new Mesh;
        create_centre_elements();
      }

      // Create the contact line elements
      // Contact_line_mesh_pt = new Mesh;

      // Add the sub meshes and data to the problem
      add_global_data(Tilde_external_pressure_data_pt);
      this->add_sub_mesh(Fluid_mesh_pt);
      this->add_sub_mesh(Free_surface_mesh_pt);
      if (Slip_mesh_pt)
      {
        this->add_sub_mesh(Slip_mesh_pt);
      }
      if (Flux_mesh_pt)
      {
        add_sub_mesh(Flux_mesh_pt);
      }
      if (Contact_line_mesh_pt)
      {
        create_contact_line_elements();
        add_sub_mesh(Contact_line_mesh_pt);
      }
      if (this->Azimuthal_mode_number == 0)
      {
        add_sub_mesh(Volume_mesh_pt);
        add_sub_mesh(Integral_mesh_pt);
      }
      if (this->Azimuthal_mode_number == 1)
      {
        add_sub_mesh(Centre_mesh_pt);
      }
      if (Net_flux_mesh_pt)
      {
        add_sub_mesh(Net_flux_mesh_pt);
      }
      // add_global_data(Flux_lagrange_multiplier_data_pt);

      // Set the boundary conditions
      set_boundary_conditions();

      // Set up the connections to the base state
      set_up_overlapping_domain_functions();

      // Build the global mesh
      build_global_mesh();

      // Set up the equation numbering so we are ready to solve the problem.
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;

      // Pin the horizontal mesh deformation.
      if (this->Azimuthal_mode_number > 0)
      {
        pin_horizontal_mesh_deformation();
      }
    }

    // Destructor
    ~PerturbedLinearStabilityCapProblem()
    {
      // delete time_stepper_pt;
      if (Constitutive_law_pt) delete Constitutive_law_pt;
      if (ElasticityTensor_pt) delete ElasticityTensor_pt;
      if (Free_surface_mesh_pt) delete Free_surface_mesh_pt;
      if (Slip_mesh_pt) delete Slip_mesh_pt;
      if (Flux_mesh_pt) delete Flux_mesh_pt;
      if (Integral_mesh_pt) delete Integral_mesh_pt;
      if (Volume_mesh_pt) delete Volume_mesh_pt;
      if (Centre_mesh_pt) delete Centre_mesh_pt;
    }

    void setup_bulk_elements()
    {
      // Loop over the elements to set the consitutive law and jacobian
      unsigned n_bulk = Fluid_mesh_pt->nelement();
      for (unsigned e = 0; e < n_bulk; e++)
      {
        // Upcast from GeneralisedElement to the present element
        PERTURBED_ELEMENT* el_pt =
          dynamic_cast<PERTURBED_ELEMENT*>(Fluid_mesh_pt->element_pt(e));


        // Set the Reynolds number
        el_pt->re_pt() = &Re;

        // Set the Womersley number
        el_pt->re_st_pt() = &ReSt;

        el_pt->azimuthal_mode_number_pt() = &Azimuthal_mode_number;

        // el_pt->constitutive_law_pt() = Constitutive_law_pt;
        // el_pt->enable_evaluate_jacobian_by_fd();
        el_pt->elasticity_tensor_pt() = ElasticityTensor_pt;
        el_pt->disable_inertia();

        // // Set the product of the Reynolds number and the inverse of the
        // // Froude number
        el_pt->re_invfr_pt() = &ReInvFr;

        // // Set the direction of gravity
        el_pt->g_pt() = &Global_Physical_Parameters::G;

        el_pt->set_time_stepper_pt(time_stepper_pt());
        // el_pt->disable_ALE();
      }
    }

    void create_free_surface_elements()
    {
      // Loop over the free surface boundary and create the "interface
      // elements
      unsigned b = Free_surface_boundary_id;

      // How many bulk fluid elements are adjacent to boundary b?
      unsigned n_element = Fluid_mesh_pt->nboundary_element(b);

      // Loop over the bulk fluid elements adjacent to boundary b?
      for (unsigned e = 0; e < n_element; e++)
      {
        // Get pointer to the bulk fluid element that is
        // adjacent to boundary b
        PERTURBED_ELEMENT* bulk_elem_pt = dynamic_cast<PERTURBED_ELEMENT*>(
          Fluid_mesh_pt->boundary_element_pt(b, e));

        // Find the index of the face of element e along boundary b
        int face_index = Fluid_mesh_pt->face_index_at_boundary(b, e);

        // Create new element
        FREE_SURFACE_ELEMENT* el_pt = new FREE_SURFACE_ELEMENT(
          bulk_elem_pt, face_index, time_stepper_pt(), Lagrange_id::Kinematic);

        // Add the appropriate boundary number
        el_pt->set_boundary_number_in_bulk_mesh(b);

        // Add the capillary number
        el_pt->ca_pt() = &Ca;
        el_pt->st_pt() = &St;

        el_pt->azimuthal_mode_number_pt() = &Azimuthal_mode_number;

        const unsigned value_index = 0;
        el_pt->set_external_pressure_data_pt(Tilde_external_pressure_data_pt,
                                             value_index);
        // Add it to the mesh
        Free_surface_mesh_pt->add_element_pt(el_pt);
      }
    }

    void create_slip_elements()
    {
      // Loop over the free surface boundary and create the "interface
      // elements
      unsigned b = Outer_boundary_with_slip_id;

      // How many bulk fluid elements are adjacent to boundary b?
      unsigned n_element = Fluid_mesh_pt->nboundary_element(b);

      // Loop over the bulk fluid elements adjacent to boundary b?
      for (unsigned e = 0; e < n_element; e++)
      {
        // Get pointer to the bulk fluid element that is
        // adjacent to boundary b
        PERTURBED_ELEMENT* bulk_elem_pt = dynamic_cast<PERTURBED_ELEMENT*>(
          Fluid_mesh_pt->boundary_element_pt(b, e));

        // Find the index of the face of element e along boundary b
        int face_index = Fluid_mesh_pt->face_index_at_boundary(b, e);

        // Create new element
        SLIP_ELEMENT* el_pt = new SLIP_ELEMENT(bulk_elem_pt, face_index);

        // Add the appropriate boundary number
        el_pt->set_boundary_number_in_bulk_mesh(b);

        el_pt->slip_fct_pt() = &Slip_Parameters::prescribed_slip_fct;

        // Add it to the mesh
        Slip_mesh_pt->add_element_pt(el_pt);
      }
    }

    void create_integral_elements()
    {
      // VOLUME_ELEMENT* int_el_pt = new VOLUME_ELEMENT;
      // Volume_mesh_pt->add_element_pt(int_el_pt);
      // Volume_data_pt = int_el_pt->get_data_pt();


      // if (this->Azimuthal_mode_number == 0)
      // {
      //   VOLUME_CONSTRAINT_ELEMENT* data_el_pt = new
      //   VOLUME_CONSTRAINT_ELEMENT(
      //     Volume_data_pt, Tilde_external_pressure_data_pt);
      //   Integral_mesh_pt->add_element_pt(data_el_pt);
      // }

      // How many bulk fluid elements are adjacent to boundary b?
      unsigned n_element =
        Fluid_mesh_pt->nboundary_element(Free_surface_boundary_id);

      // Loop over the bulk fluid elements adjacent to boundary b?
      for (unsigned e = 0; e < n_element; e++)
      {
        // Get pointer to the bulk fluid element that is
        // adjacent to boundary b
        PERTURBED_ELEMENT* bulk_elem_pt = dynamic_cast<PERTURBED_ELEMENT*>(
          Fluid_mesh_pt->boundary_element_pt(Free_surface_boundary_id, e));

        // Find the index of the face of element e along boundary b
        int face_index =
          Fluid_mesh_pt->face_index_at_boundary(Free_surface_boundary_id, e);

        // Create new element
        VOLUME_COMPUTATION_ELEMENT* el_pt = new VOLUME_COMPUTATION_ELEMENT(
          bulk_elem_pt, face_index, Tilde_external_pressure_data_pt, 0);
        // VOLUME_COMPUTATION_ELEMENT* el_pt = new VOLUME_COMPUTATION_ELEMENT(
        //   bulk_elem_pt, face_index, Volume_data_pt, 0);

        // Add it to the mesh
        Integral_mesh_pt->add_element_pt(el_pt);
      }
    }

    void create_flux_elements()
    {
      NET_FLUX_ELEMENT* net_flux_el_pt =
        new NET_FLUX_ELEMENT(&Flux_Parameters::flux_fct, time_pt());
      Net_flux_mesh_pt->add_element_pt(net_flux_el_pt);
      Flux_lagrange_multiplier_data_pt = net_flux_el_pt->internal_data_pt(0);

      // How many bulk fluid elements are adjacent to boundary b?
      unsigned n_element = Fluid_mesh_pt->nboundary_element(Upper_boundary_id);

      // Loop over the bulk fluid elements adjacent to boundary b?
      for (unsigned e = 0; e < n_element; e++)
      {
        // Get pointer to the bulk fluid element that is
        // adjacent to boundary b
        PERTURBED_ELEMENT* bulk_elem_pt = dynamic_cast<PERTURBED_ELEMENT*>(
          Fluid_mesh_pt->boundary_element_pt(Upper_boundary_id, e));

        // Find the index of the face of element e along boundary b
        int face_index =
          Fluid_mesh_pt->face_index_at_boundary(Upper_boundary_id, e);

        // Create new element
        FLUX_COMPUTATION_ELEMENT* el_pt = new FLUX_COMPUTATION_ELEMENT(
          bulk_elem_pt, face_index, Flux_lagrange_multiplier_data_pt, 0);

        // Add it to the mesh
        Flux_mesh_pt->add_element_pt(el_pt);
      }
    }

    void create_centre_elements()
    {
      oomph_info << "create_centre_elements" << std::endl;

      // Loop over all the inner boundary elements and hijack the cosine
      // horizontal dofs at the inner nodes.
      // Visit each node only once.

      // Loop over elements
      // const unsigned n_element =
      //  Fluid_mesh_pt->nboundary_element(Inner_boundary_id);
      const unsigned n_element = Fluid_mesh_pt->nelement();
      Vector<Node*> visited_nodes;
      for (unsigned n = 0; n < n_element; n++)
      {
        // PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>(
        //   Fluid_mesh_pt->boundary_element_pt(Inner_boundary_id, n));
        PERTURBED_ELEMENT* el_pt =
          dynamic_cast<PERTURBED_ELEMENT*>(Fluid_mesh_pt->element_pt(n));

        // Loop over nodes
        unsigned n_node = el_pt->nnode();
        for (unsigned i_node = 0; i_node < n_node; i_node++)
        {
          Node* current_node_pt = el_pt->node_pt(i_node);
          if (current_node_pt->is_on_boundary(Inner_boundary_id))
          {
            Vector<Node*>::iterator found_node_pt =
              find(visited_nodes.begin(), visited_nodes.end(), current_node_pt);

            // Hijack the velocity values
            Data* hijacked_data_pt =
              el_pt->hijack_nodal_value(i_node, uc_index);
            add_global_data(hijacked_data_pt);
            hijacked_data_pt = el_pt->hijack_nodal_value(i_node, us_index);
            add_global_data(hijacked_data_pt);

            if (found_node_pt == visited_nodes.end())
            {
              visited_nodes.push_back(current_node_pt);

              // add to the new elements
              // UC
              SYMMETRY_CONDITION_ELEMENT* new_el_pt =
                new SYMMETRY_CONDITION_ELEMENT(current_node_pt);
              Centre_mesh_pt->add_element_pt(new_el_pt);
            }
          }
        }
      }

      // Hijack the free surface element at the centre
      unsigned corner_element_node_index = 0;
      FREE_SURFACE_ELEMENT* fs_el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
        find_corner_free_surface_element(corner_element_node_index));
      fs_el_pt->hijack_nodal_value(corner_element_node_index, uc_index);
      fs_el_pt->hijack_nodal_value(corner_element_node_index, us_index);

      std::cout << visited_nodes.size() << " visited of "
           << Fluid_mesh_pt->nboundary_node(Inner_boundary_id) << std::endl;
    }

    void create_contact_line_elements()
    {
      oomph_info << "create_contact_line_elements" << std::endl;
      const unsigned n_element = Free_surface_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(n));

        bool is_contact_element = false;
        int node_index = -1;
        unsigned n_node = el_pt->nnode();
        for (unsigned i_node = 0; i_node < n_node; i_node++)
        {
          // If the node is on the free surface boundary as well then ...
          if (el_pt->node_pt(i_node)->is_on_boundary(
                Outer_boundary_with_slip_id))
          {
            is_contact_element = true;
            node_index = i_node;
            break;
          }
        }

        if (is_contact_element)
        {
          Vector<double> gradient = el_pt->displacement_gradient(1, 0);

          if (node_index == 0)
          {
            CONTACT_LINE_ELEMENT* new_el_pt =
              new CONTACT_LINE_ELEMENT(el_pt, -1);

            new_el_pt->ca_pt() = &Ca;
            new_el_pt->contact_angle_pt() = &Contact_angle;
            new_el_pt->wall_unit_normal_fct_pt() =
              &Global_Physical_Parameters::wall_unit_normal_fct;
            new_el_pt->azimuthal_mode_number_pt() = &Azimuthal_mode_number;
            // new_el_pt->external_element_pt(0, 0) = el_pt;
            new_el_pt->add_external_data(el_pt->node_pt(1));
            new_el_pt->add_external_data(el_pt->node_pt(2));
            Contact_line_mesh_pt->add_element_pt(new_el_pt);
          }
          else
          {
            CONTACT_LINE_ELEMENT* new_el_pt =
              new CONTACT_LINE_ELEMENT(el_pt, 1);

            new_el_pt->ca_pt() = &Ca;
            new_el_pt->contact_angle_pt() = &Contact_angle;
            new_el_pt->wall_unit_normal_fct_pt() =
              &Global_Physical_Parameters::wall_unit_normal_fct;
            new_el_pt->azimuthal_mode_number_pt() = &Azimuthal_mode_number;
            // new_el_pt->external_element_pt(0, 0) = el_pt;
            new_el_pt->add_external_data(el_pt->node_pt(0));
            new_el_pt->add_external_data(el_pt->node_pt(1));
            Contact_line_mesh_pt->add_element_pt(new_el_pt);
          }
        }
      }
    }

    // Make the problem steady.
    // This will change what boundary conditions we are imposing
    void make_steady()
    {
      oomph_info << "make_steady" << std::endl;

      this->Is_steady = true;

      set_boundary_conditions();

      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    // Make the problem unsteady.
    // This will change what boundary conditions we are imposing
    void make_unsteady()
    {
      oomph_info << "make_unsteady" << std::endl;

      this->Is_steady = false;

      set_boundary_conditions();

      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    const bool is_steady()
    {
      return this->Is_steady;
    }

    // Set the problem boundary conditions
    void set_boundary_conditions()
    {
      // Simplify by pinning the horizontal "deformation"/deviation
      // pin_horizontal_mesh_deformation();

      // If m = 0
      if (this->Azimuthal_mode_number == 0)
      {
        pin_sine_components();
        pin_azimuthal_velocity();
      }

      // Set each of the boundary conditions
      set_outer_boundary_condition();
      set_upper_boundary_condition();
      set_inner_boundary_condition();
      set_free_surface_boundary_condition();

      if (this->Azimuthal_mode_number == 0)
      {
        // If the problem is steady
        if (this->Is_steady)
        {
          Tilde_external_pressure_data_pt->unpin(0);
          pin_flux_constraint();
          pin_interior_pressure();
        }
        else
        {
          Flux_lagrange_multiplier_data_pt->unpin(0);
          pin_volume_constraint();
          unpin_interior_pressure();
        }
      }
      else
      {
        pin_volume_constraint();
        // pin_interior_pressure();
        //  pin_flux_constraint();
      }

      // If we have no slip
      if (Slip_Parameters::slip_length == 0)
      {
        // Impose a fixed contact line
        pin_velocity_on_boundary(Outer_boundary_with_slip_id, wc_index);
        pin_velocity_on_boundary(Outer_boundary_with_slip_id, ws_index);
        pin_contact_line();
      }

      // pin_fluid();
      // pin_all_pressure();
      // pin_vertical_mesh_deformation();

      // set_constant_lagrange_free_surface_boundary_condition(0.001, 0.0);
      // Tilde_external_pressure_data_pt->pin(0);

      // pin_fluid_boundary(Free_surface_boundary_id);
      // pin_fluid_boundary(Upper_boundary_id);
      // pin_fluid_boundary(Outer_boundary_with_slip_id);
      // pin_fluid_boundary(Inner_boundary_id);
    }

    void isolate_volume()
    {
      pin_fluid();
      pin_horizontal_mesh_deformation();
      pin_vertical_mesh_deformation();
      set_constant_lagrange_free_surface_boundary_condition(0.0, 0.0);
      pin_flux_constraint();
      Tilde_external_pressure_data_pt->pin(0);

      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void pin_flux_constraint()
    {
      Flux_lagrange_multiplier_data_pt->pin(0);
      Flux_lagrange_multiplier_data_pt->set_value(0, 0.0);
    }

  public:
    void pin_fluid()
    {
      unsigned n_element = Fluid_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)

      {
        dynamic_cast<PERTURBED_ELEMENT*>(Fluid_mesh_pt->element_pt(n))->pin();
      }

      set_constant_lagrange_free_surface_boundary_condition(0.0, 0.0);

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void pin_solid()
    {
      pin_horizontal_mesh_deformation();
      pin_vertical_mesh_deformation();
      pin_volume_constraint();

      set_constant_lagrange_free_surface_boundary_condition(0.0, 0.0);

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void pin_solid_boundaries()
    {
      for (unsigned j = 0; j < 2; j++)
      {
        for (unsigned b = 0; b < 4; b++)
        {
          for (unsigned cs = 0; cs < 2; cs++)
          {
            pin_position_on_boundary(j, b, cs);
          }
        }
      }
      pin_volume_constraint();

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void pin_fs_solid_boundary()
    {
      for (unsigned j = 0; j < 2; j++)
      {
        for (unsigned b = 0; b < 4; b++)
        {
          for (unsigned cs = 0; cs < 2; cs++)
          {
            if (b == Free_surface_boundary_id)
            {
              pin_position_on_boundary(j, b, cs);
            }
          }
        }
      }
      pin_volume_constraint();

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void pin_solid_boundaries_except_upper()
    {
      for (unsigned j = 0; j < 2; j++)
      {
        for (unsigned b = 0; b < 4; b++)
        {
          for (unsigned cs = 0; cs < 2; cs++)
          {
            if (b != Upper_boundary_id)
            {
              pin_position_on_boundary(j, b, cs);
            }
          }
        }
      }
      pin_volume_constraint();

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void pin_solid_boundaries(const unsigned& b)
    {
      for (unsigned j = 0; j < 2; j++)
      {
        for (unsigned cs = 0; cs < 2; cs++)
        {
          pin_position_on_boundary(j, b, cs);
        }
      }
      if (b == Free_surface_boundary_id)
      {
        pin_volume_constraint();
      }

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

  private:
    void pin_volume_constraint()
    {
      Tilde_external_pressure_data_pt->pin(0);
      Tilde_external_pressure_data_pt->set_value(0, 0.0);
    }

  public:
    void pin_position_on_boundary(const unsigned& direction,
                                  const unsigned& b,
                                  const unsigned& cs)
    {
      unsigned n_element;
      n_element = Fluid_mesh_pt->nboundary_element(b);
      for (unsigned n = 0; n < n_element; n++)
      {
        PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>(
          Fluid_mesh_pt->boundary_element_pt(b, n));
        for (unsigned m = 0; m < 6; m++)
        {
          if (el_pt->node_pt(m)->is_on_boundary(b))
          {
            el_pt->pin_Xhat(m, direction, cs);
          }
        }
      }
    }

    void pin_horizontal_mesh_deformation()
    {
      unsigned n_element = Fluid_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)

      {
        PERTURBED_ELEMENT* el_pt =
          dynamic_cast<PERTURBED_ELEMENT*>(Fluid_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 6; m++)
        {
          el_pt->pin_Xhat(m, 0, 0);
          el_pt->pin_Xhat(m, 0, 1);
        }
      }
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void pin_vertical_mesh_deformation()
    {
      unsigned n_element = Fluid_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)

      {
        PERTURBED_ELEMENT* el_pt =
          dynamic_cast<PERTURBED_ELEMENT*>(Fluid_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 6; m++)
        {
          el_pt->pin_Xhat(m, 1, 0);
          el_pt->pin_Xhat(m, 1, 1);
        }
      }
    }

    void pin_azimuthal_velocity()
    {
      unsigned n_node = Fluid_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        Fluid_mesh_pt->node_pt(n)->pin(vc_index);
        Fluid_mesh_pt->node_pt(n)->pin(vs_index);
      }
    }

    void pin_sine_components()
    {
      unsigned n_node = Fluid_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        Fluid_mesh_pt->node_pt(n)->pin(us_index);
        Fluid_mesh_pt->node_pt(n)->pin(vs_index);
        Fluid_mesh_pt->node_pt(n)->pin(ws_index);
      }

      unsigned n_element = Fluid_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        PERTURBED_ELEMENT* el_pt =
          dynamic_cast<PERTURBED_ELEMENT*>(Fluid_mesh_pt->element_pt(n));
        el_pt->pin_pressure(1);
        for (unsigned m = 0; m < 6; m++)
        {
          el_pt->pin_Xhat(m, 0, 1);
          el_pt->pin_Xhat(m, 1, 1);
        }
      }

      n_element = Free_surface_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 3; m++)
        {
          el_pt->pin_Xhat(m, 0, 1);
          el_pt->pin_Xhat(m, 1, 1);
          el_pt->pin_lagrange_multiplier(m, 1, 0.0);
        }
      }
    }

    void pin_all_pressure()
    {
      unsigned n_element = Fluid_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        PERTURBED_ELEMENT* el_pt =
          dynamic_cast<PERTURBED_ELEMENT*>(Fluid_mesh_pt->element_pt(n));
        el_pt->pin_pressure(0);
        el_pt->pin_pressure(1);
      }
    }

    void pin_interior_pressure()
    {
      unsigned element_index = 0;
      unsigned node_index = 0;
      unsigned cs_index = 0;
      find_corner_bulk_element_and_node(Outer_boundary_with_slip_id,
                                        Upper_boundary_id,
                                        element_index,
                                        node_index);

      PERTURBED_ELEMENT* el_pt =
        dynamic_cast<PERTURBED_ELEMENT*>(Fluid_mesh_pt->boundary_element_pt(
          Outer_boundary_with_slip_id, element_index));
      el_pt->pin_pressure(node_index, cs_index);
    }

    void unpin_interior_pressure()
    {
      unsigned element_index = 0;
      unsigned node_index = 0;
      unsigned cs_index = 0;
      find_corner_bulk_element_and_node(Outer_boundary_with_slip_id,
                                        Upper_boundary_id,
                                        element_index,
                                        node_index);

      PERTURBED_ELEMENT* el_pt =
        dynamic_cast<PERTURBED_ELEMENT*>(Fluid_mesh_pt->boundary_element_pt(
          Outer_boundary_with_slip_id, element_index));
      el_pt->unpin_pressure(node_index, cs_index);
    }

    void find_corner_bulk_element_and_node(const unsigned& boundary_1_id,
                                           const unsigned& boundary_2_id,
                                           unsigned& element_index,
                                           unsigned& node_index)
    {
      unsigned n_boundary_element =
        Fluid_mesh_pt->nboundary_element(boundary_1_id);
      for (unsigned e = 0; e < n_boundary_element; e++)
      {
        // Locally cache the element pointer
        FiniteElement* bulk_el_pt =
          Fluid_mesh_pt->boundary_element_pt(boundary_1_id, e);

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

    // Find corner node and return whether it has been found
    bool find_corner_node(const unsigned& first_boundary_id,
                          const unsigned& second_boundary_id,
                          SolidNode*& node_pt)
    {
      unsigned n_nod = Fluid_mesh_pt->nboundary_node(first_boundary_id);
      for (unsigned inod = 0; inod < n_nod; inod++)
      {
        // Get boundary node
        Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(first_boundary_id, inod);

        if (nod_pt->is_on_boundary(second_boundary_id))
        {
          node_pt = dynamic_cast<SolidNode*>(nod_pt);

          return true;
        }
      }
      return false;
    }

    void set_outer_boundary_condition()
    {
      oomph_info << "set_outer_boundary_condition" << std::endl;
      // Loop over the nodes on the boundary
      unsigned n_boundary_node;
      n_boundary_node =
        Fluid_mesh_pt->nboundary_node(Outer_boundary_with_slip_id);
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        Fluid_mesh_pt->boundary_node_pt(Outer_boundary_with_slip_id, n)
          ->pin(uc_index);
        Fluid_mesh_pt->boundary_node_pt(Outer_boundary_with_slip_id, n)
          ->pin(vc_index);
        // Fluid_mesh_pt->boundary_node_pt(Outer_boundary_with_slip_id, n)
        //  ->pin(wc_index);
        Fluid_mesh_pt->boundary_node_pt(Outer_boundary_with_slip_id, n)
          ->pin(us_index);
        Fluid_mesh_pt->boundary_node_pt(Outer_boundary_with_slip_id, n)
          ->pin(vs_index);
        // Fluid_mesh_pt->boundary_node_pt(Outer_boundary_with_slip_id, n)
        //  ->pin(ws_index);
      }

      unsigned n_element;
      n_element = Fluid_mesh_pt->nboundary_element(Outer_boundary_with_slip_id);
      for (unsigned n = 0; n < n_element; n++)
      {
        PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>(
          Fluid_mesh_pt->boundary_element_pt(Outer_boundary_with_slip_id, n));
        for (unsigned m = 0; m < 6; m++)
        {
          if (el_pt->node_pt(m)->is_on_boundary(Outer_boundary_with_slip_id))
          {
            for (unsigned j = 0; j < 2; j++)
            {
              el_pt->pin_Xhat(m, 0, j);
            }
          }
        }
      }
    }

    void pin_velocity_on_boundary(const unsigned boundary_id,
                                  const unsigned velocity_id)
    {
      // Loop over the nodes on the boundary
      unsigned n_boundary_node = Fluid_mesh_pt->nboundary_node(boundary_id);
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        Fluid_mesh_pt->boundary_node_pt(boundary_id, n)->pin(velocity_id);
      }
    }

    void set_velocity_on_boundary(const unsigned boundary_id,
                                  const unsigned velocity_id,
                                  const double value)
    {
      // Loop over the nodes on the boundary
      unsigned n_boundary_node = Fluid_mesh_pt->nboundary_node(boundary_id);
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        Fluid_mesh_pt->boundary_node_pt(boundary_id, n)
          ->set_value(velocity_id, value);
      }
    }

    void set_upper_boundary_condition()
    {
      oomph_info << "set_upper_boundary_condition" << std::endl;
      // Velocity constraints
      unsigned n_boundary_node =
        Fluid_mesh_pt->nboundary_node(Upper_boundary_id);
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        Fluid_mesh_pt->boundary_node_pt(Upper_boundary_id, n)->pin(uc_index);
        Fluid_mesh_pt->boundary_node_pt(Upper_boundary_id, n)->pin(vc_index);
        // Fluid_mesh_pt->boundary_node_pt(Upper_boundary_id,
        // n)->pin(wc_index);
        Fluid_mesh_pt->boundary_node_pt(Upper_boundary_id, n)->pin(us_index);
        Fluid_mesh_pt->boundary_node_pt(Upper_boundary_id, n)->pin(vs_index);
        // Fluid_mesh_pt->boundary_node_pt(Upper_boundary_id,
        // n)->pin(ws_index);
      }

      // Displacement constraints
      unsigned n_element;
      n_element = Fluid_mesh_pt->nboundary_element(Upper_boundary_id);
      for (unsigned n = 0; n < n_element; n++)
      {
        PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>(
          Fluid_mesh_pt->boundary_element_pt(Upper_boundary_id, n));
        for (unsigned m = 0; m < 6; m++)
        {
          if (el_pt->node_pt(m)->is_on_boundary(Upper_boundary_id))
          {
            for (unsigned j = 0; j < 2; j++)
            {
              el_pt->pin_Xhat(m, 1, j);
            }
          }
        }
      }
    }

    void pin_fluid_boundary(const unsigned& boundary_id)
    {
      unsigned n_boundary_node = Fluid_mesh_pt->nboundary_node(boundary_id);
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        Fluid_mesh_pt->boundary_node_pt(boundary_id, n)->pin(uc_index);
        Fluid_mesh_pt->boundary_node_pt(boundary_id, n)->pin(vc_index);
        Fluid_mesh_pt->boundary_node_pt(boundary_id, n)->pin(wc_index);
        Fluid_mesh_pt->boundary_node_pt(boundary_id, n)->pin(us_index);
        Fluid_mesh_pt->boundary_node_pt(boundary_id, n)->pin(vs_index);
        Fluid_mesh_pt->boundary_node_pt(boundary_id, n)->pin(ws_index);
      }
      if (boundary_id == Free_surface_boundary_id)
      {
        set_constant_lagrange_free_surface_boundary_condition(0.0, 0.0);
      }

      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void set_inner_boundary_condition()
    {
      oomph_info << "set_inner_boundary_condition" << std::endl;
      unsigned n_boundary_node =
        Fluid_mesh_pt->nboundary_node(Inner_boundary_id);
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        if (Azimuthal_mode_number == 0)
        {
          Fluid_mesh_pt->boundary_node_pt(Inner_boundary_id, n)->pin(uc_index);
          Fluid_mesh_pt->boundary_node_pt(Inner_boundary_id, n)->pin(us_index);
          Fluid_mesh_pt->boundary_node_pt(Inner_boundary_id, n)->pin(vc_index);
          Fluid_mesh_pt->boundary_node_pt(Inner_boundary_id, n)->pin(vs_index);
        }
        else if (Azimuthal_mode_number == 1)
        {
          Fluid_mesh_pt->boundary_node_pt(Inner_boundary_id, n)->pin(wc_index);
          Fluid_mesh_pt->boundary_node_pt(Inner_boundary_id, n)->pin(ws_index);
        }
        else if (Azimuthal_mode_number > 1)
        {
          Fluid_mesh_pt->boundary_node_pt(Inner_boundary_id, n)->pin(uc_index);
          Fluid_mesh_pt->boundary_node_pt(Inner_boundary_id, n)->pin(us_index);
          Fluid_mesh_pt->boundary_node_pt(Inner_boundary_id, n)->pin(vc_index);
          Fluid_mesh_pt->boundary_node_pt(Inner_boundary_id, n)->pin(vs_index);
          Fluid_mesh_pt->boundary_node_pt(Inner_boundary_id, n)->pin(wc_index);
          Fluid_mesh_pt->boundary_node_pt(Inner_boundary_id, n)->pin(ws_index);
        }
      }

      unsigned n_element;
      n_element = Fluid_mesh_pt->nboundary_element(Inner_boundary_id);
      for (unsigned n = 0; n < n_element; n++)
      {
        PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>(
          Fluid_mesh_pt->boundary_element_pt(Inner_boundary_id, n));
        for (unsigned m = 0; m < 6; m++)
        {
          if (el_pt->node_pt(m)->is_on_boundary(Inner_boundary_id))
          {
            for (unsigned j = 0; j < 2; j++)
            {
              el_pt->pin_Xhat(m, 0, j);
              if (Azimuthal_mode_number > 0)
              {
                el_pt->pin_Xhat(m, 1, j);
              }
            }
          }
        }
      }

      if (Azimuthal_mode_number == 1)
      {
        // pin_pressure_along_inner_boundary();
      }
    }

    void pin_pressure_along_inner_boundary()
    {
      unsigned n_element = Fluid_mesh_pt->nboundary_element(Inner_boundary_id);
      for (unsigned n = 0; n < n_element; n++)
      {
        PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>(
          Fluid_mesh_pt->boundary_element_pt(Inner_boundary_id, n));
        for (unsigned m = 0; m < 3; m++)
        {
          if (el_pt->node_pt(m)->is_on_boundary(Inner_boundary_id))

          {
            el_pt->node_pt(m)->pin(4 + 6);
            el_pt->node_pt(m)->pin(4 + 6 + 1);
          }
        }
      }
    }

    void set_free_surface_boundary_condition()
    {
      oomph_info << "set_free_surface_boundary_condition" << std::endl;

      if (this->Azimuthal_mode_number == 0)
      {
        // pin_centre_corner_lagrange_multipler();
      }
      else if (this->Azimuthal_mode_number == 1)
      {
        // pin_centre_corner_displacement();
        // pin_centre_corner_horizontal_displacement();
        // pin_centre_corner_velocities();
        pin_centre_corner_lagrange_multipler();
      }
      else if (this->Azimuthal_mode_number > 1)
      {
        pin_centre_corner_lagrange_multipler();
        pin_centre_corner_velocities();
      }
    }

    void pin_centre_pressure()
    {
      oomph_info << "pin_centre_pressure" << std::endl;
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
            el_pt->node_pt(m)->pin(4 + 6);
            el_pt->node_pt(m)->pin(4 + 6 + 1);
          }
        }
      }
    }

    void pin_centre_corner_lagrange_multipler()
    {
      oomph_info << "pin_centre_corner_lagrange_multipler" << std::endl;
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
            el_pt->pin_lagrange_multiplier(m, 0, 0.0);
            el_pt->pin_lagrange_multiplier(m, 1, 0.0);
          }
        }
      }
    }

    void pin_centre_corner_displacement()
    {
      oomph_info << "pin_centre_corner_displacement" << std::endl;
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
            // Pin RC and RS
            el_pt->pin_Xhat(m, 0, 0);
            el_pt->pin_Xhat(m, 0, 1);
            el_pt->pin_Xhat(m, 1, 0);
            el_pt->pin_Xhat(m, 1, 1);
          }
        }
      }
    }

    void pin_centre_corner_horizontal_displacement()
    {
      oomph_info << "pin_centre_corner_displacement" << std::endl;
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
            // Pin RC and RS
            el_pt->pin_Xhat(m, 0, 0);
            el_pt->pin_Xhat(m, 0, 1);
          }
        }
      }
    }

    void pin_centre_corner_velocities()
    {
      oomph_info << "pin_centre_corner_velocities" << std::endl;
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
            // Pin RC and RS
            el_pt->node_pt(m)->pin(uc_index);
            el_pt->node_pt(m)->pin(us_index);
            el_pt->node_pt(m)->pin(vc_index);
            el_pt->node_pt(m)->pin(vs_index);
            el_pt->node_pt(m)->pin(wc_index);
            el_pt->node_pt(m)->pin(ws_index);
          }
        }
      }
    }

    void pin_centre()
    {
      oomph_info << "pin_centre" << std::endl;
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
            el_pt->pin_lagrange_multiplier(m, 0, 0.0);
            el_pt->pin_lagrange_multiplier(m, 1, 0.0);
            // Pin RC and RS
            el_pt->pin_Xhat(m, 0, 0);
            el_pt->pin_Xhat(m, 0, 1);
            el_pt->pin_Xhat(m, 1, 0);
            el_pt->pin_Xhat(m, 1, 1);

            el_pt->node_pt(m)->pin(0);
            el_pt->node_pt(m)->pin(1);
            el_pt->node_pt(m)->pin(2);
            el_pt->node_pt(m)->pin(3);
            el_pt->node_pt(m)->pin(4);
            el_pt->node_pt(m)->pin(5);
          }
        }
      }

      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void pin_centre_element()
    {
      oomph_info << "pin_centre_element" << std::endl;
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
              el_pt->pin_lagrange_multiplier(p, 0, 0.0);
              el_pt->pin_lagrange_multiplier(p, 1, 0.0);
              // Pin RC and RS
              el_pt->pin_Xhat(p, 0, 0);
              el_pt->pin_Xhat(p, 0, 1);
              el_pt->pin_Xhat(p, 1, 0);
              el_pt->pin_Xhat(p, 1, 1);

              el_pt->node_pt(p)->pin(0);
              el_pt->node_pt(p)->pin(1);
              el_pt->node_pt(p)->pin(2);
              el_pt->node_pt(p)->pin(3);
              el_pt->node_pt(p)->pin(4);
              el_pt->node_pt(p)->pin(5);
            }
          }
        }
      }
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
            el_pt->pin_lagrange_multiplier(m, 0, 0.0);
            el_pt->pin_lagrange_multiplier(m, 1, 0.0);
            // Pin RC and RS
            el_pt->pin_Xhat(m, 0, 0);
            el_pt->pin_Xhat(m, 0, 1);
            el_pt->pin_Xhat(m, 1, 0);
            el_pt->pin_Xhat(m, 1, 1);

            el_pt->node_pt(m)->pin(0);
            el_pt->node_pt(m)->pin(1);
            el_pt->node_pt(m)->pin(2);
            el_pt->node_pt(m)->pin(3);
            el_pt->node_pt(m)->pin(4);
            el_pt->node_pt(m)->pin(5);
          }
        }
      }

      // Rebuild the global mesh
      this->rebuild_global_mesh();

      // Setup all the equation numbering and look-up schemes
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }


    void pin_contact_line_and_lagrange_multiplier()
    {
      oomph_info << "pin_contact_line_and_lagrange_multiplier" << std::endl;
      unsigned n_element;
      n_element = Free_surface_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 3; m++)
        {
          if (el_pt->node_pt(m)->is_on_boundary(Outer_boundary_with_slip_id))
          {
            el_pt->pin_lagrange_multiplier(m, 0, 0.0);
            el_pt->pin_lagrange_multiplier(m, 1, 0.0);
            // Pin RC and RS
            el_pt->pin_Xhat(m, 0, 0);
            el_pt->pin_Xhat(m, 0, 1);
            el_pt->pin_Xhat(m, 1, 0);
            el_pt->pin_Xhat(m, 1, 1);
          }
        }
      }
    }

    void pin_contact_line()
    {
      oomph_info << "pin_contact_line" << std::endl;
      unsigned n_element;
      n_element = Free_surface_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 3; m++)
        {
          if (el_pt->node_pt(m)->is_on_boundary(Outer_boundary_with_slip_id))
          {
            el_pt->pin_lagrange_multiplier(m, 0, 0.0);
            el_pt->pin_lagrange_multiplier(m, 1, 0.0);
            // Pin RC and RS
            el_pt->pin_Xhat(m, 0, 0);
            el_pt->pin_Xhat(m, 0, 1);
            el_pt->pin_Xhat(m, 1, 0);
            el_pt->pin_Xhat(m, 1, 1);

            el_pt->node_pt(m)->pin(0);
            el_pt->node_pt(m)->pin(1);
            el_pt->node_pt(m)->pin(2);
            el_pt->node_pt(m)->pin(3);
            el_pt->node_pt(m)->pin(4);
            el_pt->node_pt(m)->pin(5);
          }
        }
      }
    }

    void pin_sin_free_surface_outer_point()
    {
      unsigned n_element = Free_surface_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 3; m++)
        {
          if (el_pt->node_pt(m)->is_on_boundary(Outer_boundary_with_slip_id))
          {
            el_pt->pin_Xhat(m, 1, 1);
          }
        }
      }
    }

    void pin_horizontal_free_surface_displacements()
    {
      unsigned n_element = Free_surface_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 3; m++)
        {
          // // Pin RC and RS
          el_pt->pin_Xhat(m, 0, 0);
          el_pt->pin_Xhat(m, 0, 1);
        }
      }
    }

    void pin_fluid_on_free_surface()
    {
      unsigned n_boundary_node =
        Fluid_mesh_pt->nboundary_node(Free_surface_boundary_id);
      for (unsigned n = 0; n < n_boundary_node; n++)
      {
        Fluid_mesh_pt->boundary_node_pt(Free_surface_boundary_id, n)
          ->pin(uc_index);
        Fluid_mesh_pt->boundary_node_pt(Free_surface_boundary_id, n)
          ->pin(vc_index);
        Fluid_mesh_pt->boundary_node_pt(Free_surface_boundary_id, n)
          ->pin(wc_index);
        Fluid_mesh_pt->boundary_node_pt(Free_surface_boundary_id, n)
          ->pin(us_index);
        Fluid_mesh_pt->boundary_node_pt(Free_surface_boundary_id, n)
          ->pin(vs_index);
        Fluid_mesh_pt->boundary_node_pt(Free_surface_boundary_id, n)
          ->pin(ws_index);
      }
    }

    void pin_displacement_on_free_surface()
    {
      unsigned n_element = Free_surface_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 3; m++)
        {
          double r = el_pt->node_pt(m)->x(0);
          el_pt->pin_Xhat(m, 0, 0);
          el_pt->pin_Xhat(m, 1, 0);
          el_pt->pin_Xhat(m, 0, 1);
          el_pt->pin_Xhat(m, 1, 1);
        }
      }
    }

    void set_constant_lagrange_free_surface_boundary_condition(const double& a,
                                                               const double& b)
    {
      unsigned n_element = Free_surface_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 3; m++)
        {
          double r = el_pt->node_pt(m)->x(0);
          el_pt->pin_lagrange_multiplier(m, 0, a * r);
          el_pt->pin_lagrange_multiplier(m, 1, b * r);
        }
      }

      // Set up the equation numbering so we are ready to solve the problem.
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;
    }

    void set_initial_condition()
    {
      oomph_info << "set_initial_condition" << std::endl;

      perturb_z_displacement();
    }

    void perturb_free_surface()
    {
      const unsigned n_element = Free_surface_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 3; m++)
        {
          double r = el_pt->node_pt(m)->x(0);
          const double amplitude = 0.01;

          // Smoother perturbation at the origin
          if (this->Azimuthal_mode_number == 1)
          {
            if (r < 0.5)
            {
              el_pt->set_value_Xhat(m, 1, 0, 0.0);
            }
            else
            {
              el_pt->set_value_Xhat(
                m,
                1,
                0,
                0.5 * amplitude *
                  (1 - cos(4 * MathematicalConstants::Pi * (r - 0.5))));
            }
          }
          else if (this->Azimuthal_mode_number > 1)
          {
            el_pt->set_value_Xhat(
              m,
              1,
              0,
              0.5 * amplitude * (1 - cos(2 * MathematicalConstants::Pi * r)));
          }
          else
          {
            // el_pt->set_value_Xhat(
            //  m,
            //  1,
            //  0,
            //  amplitude * (1.0 / 60.0 - pow(r, 2.0) * pow(1 - r, 2.0)));
            el_pt->set_value_Xhat(m,
                                  1,
                                  0,
                                  amplitude * (-7.0 / 60.0 + pow(r, 2.0) / 2.0 -
                                               pow(r, 3.0) / 3.0));
          }
        }
      }
    }

    void perturb_free_surface2()
    {
      const unsigned n_element = Free_surface_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        FREE_SURFACE_ELEMENT* el_pt = dynamic_cast<FREE_SURFACE_ELEMENT*>(
          Free_surface_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 3; m++)
        {
          double r = el_pt->node_pt(m)->x(0);
          if (r > 0.5)
          {
            const double amplitude = 0.01;

            // el_pt->set_value_Xhat(
            //  m,
            //  1,
            //  0,
            //  amplitude * (1.0 / 60.0 - pow(r, 2.0) * pow(1 - r, 2.0)));
            el_pt->set_value_Xhat(
              m, 1, 0, amplitude * (r - 0.5) * (r - 0.5) * (1.0 - r));
          }
        }
      }
    }

    void perturb_z_displacement()
    {
      const unsigned n_element = Fluid_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        PERTURBED_ELEMENT* el_pt =
          dynamic_cast<PERTURBED_ELEMENT*>(Fluid_mesh_pt->element_pt(n));
        for (unsigned m = 0; m < 6; m++)
        {
          double r = el_pt->node_pt(m)->x(0);
          double z = el_pt->node_pt(m)->x(1);
          const double amplitude = 0.01;
          el_pt->set_value_Xhat(m,
                                1,
                                0,
                                0.5 * amplitude *
                                  (1 - cos(2 * MathematicalConstants::Pi * r)) *
                                  (3.5 - z));
        }
      }
    }

    void perturb_vertical_velocity()
    {
      const unsigned n_node = Fluid_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        double r = Fluid_mesh_pt->node_pt(n)->x(0);
        double z = Fluid_mesh_pt->node_pt(n)->x(1);
        Fluid_mesh_pt->node_pt(n)->set_value(wc_index,
                                             r * z * (1.0 - r) * (3.5 - z));
        Fluid_mesh_pt->node_pt(n)->set_value(ws_index,
                                             r * z * (1.0 - r) * (3.5 - z));
      }
    }

    void perturb_azimuthal_velocity()
    {
      const unsigned n_node = Fluid_mesh_pt->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        double r = Fluid_mesh_pt->node_pt(n)->x(0);
        double z = Fluid_mesh_pt->node_pt(n)->x(1);
        Fluid_mesh_pt->node_pt(n)->set_value(vc_index,
                                             r * z * (1.0 - r) * (1.0 - z));
        Fluid_mesh_pt->node_pt(n)->set_value(vs_index,
                                             r * z * (1.0 - r) * (1.0 - z));
      }
    }

    void set_up_overlapping_domain_functions()
    {
      unsigned n_element = Fluid_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        dynamic_cast<PERTURBED_ELEMENT*>(Fluid_mesh_pt->element_pt(n))
          ->set_base_element_pt(this->External_base_mesh_pt->element_pt(n));
      }

      n_element = Free_surface_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        GeneralisedElement* gen_el_pt =
          this->External_free_surface_mesh_pt->element_pt(n);
        dynamic_cast<FREE_SURFACE_ELEMENT*>(Free_surface_mesh_pt->element_pt(n))
          ->set_base_element_pt(gen_el_pt);
      }

      n_element = Slip_mesh_pt->nelement();
      for (unsigned n = 0; n < n_element; n++)
      {
        dynamic_cast<SLIP_ELEMENT*>(Slip_mesh_pt->element_pt(n))
          ->set_base_element_pt(
            this->External_slip_surface_mesh_pt->element_pt(n));
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
          local_max_adapt = 1; // this->Max_adapt;
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
          // this->create_restart_file();
        }
        it++;
      }
    }


    void set_directory(const std::string& dir_name)
    {
      this->doc_info().set_directory(dir_name);
    }

    void set_doc_number(const unsigned& doc_number)
    {
      this->doc_info().number() = doc_number;
    }

    // Set the bond number represented by ReInvFr
    void set_bond_number(const double& bond_number)
    {
      this->Bo = bond_number;

      update_ReInvFr();
    }

    // Set the bond number represented by ReInvFr
    double* reynolds_number_inverse_froude_number_pt()
    {
      return &this->ReInvFr;
    }

    // Set the Capillary without changing the Bond number
    void set_capillary_number(const double& capillary_number)
    {
      // Set the new capillary number
      this->Ca = capillary_number;

      update_ReInvFr();
    }

  private:
    void update_ReInvFr()
    {
      this->ReInvFr = this->Bo / this->Ca;
    }

  public:
    // Set the contact angle, not this is overridden by the contact angle
    // function.
    void set_contact_angle(const double& contact_angle)
    {
      this->Contact_angle = contact_angle;
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

    Vector<double> estimate_gradient_at_origin()
    {
      // Find the corner node
      unsigned node_index = 0;
      FiniteElement* el_pt = find_corner_free_surface_element(node_index);

      // Find another node of the boundary
      unsigned second_node_index = 1;

      const double x0 = el_pt->node_pt(node_index)->x(0);
      const double x1 = el_pt->node_pt(second_node_index)->x(0);
      const double dx = x0 - x1;
      // Compute differences
      Vector<double> gradu(6);
      for (unsigned n = 0; n < 6; n++)
      {
        const double u0 = el_pt->node_pt(node_index)->value(4 + n);
        const double u1 = el_pt->node_pt(second_node_index)->value(4 + n);
        const double du = u0 - u1;
        gradu[n] = du / dx;
      }
      return gradu;
    }

    FiniteElement* find_corner_free_surface_element(unsigned& node_index)
    {
      FiniteElement* bulk_el_pt = 0;
      unsigned n_element = Free_surface_mesh_pt->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        // Locally cache the element pointer
        bulk_el_pt =
          dynamic_cast<FiniteElement*>(Free_surface_mesh_pt->element_pt(e));

        // Read out number of nodes in the element
        unsigned n_node = bulk_el_pt->nnode();
        for (unsigned i_node = 0; i_node < n_node; i_node++)
        {
          // If the node is on the free surface boundary as well then ...
          if (bulk_el_pt->node_pt(i_node)->is_on_boundary(Inner_boundary_id))
          {
            // set the output arguments,
            node_index = i_node;

            // Return to exit both loops and end function
            return bulk_el_pt;
          }
        }
      }
      // If not found, issue warning and return anyway
      oomph_info << "Warning: No corner node found!" << std::endl;

      return bulk_el_pt;
    }


    void doc_solution()
    {
      oomph_info << "doc_solution: " << this->doc_info().number() << std::endl;
      std::ofstream output_stream;
      output_stream.precision(16);
      std::string filename;
      unsigned npts = 3;
      filename = this->doc_info().directory() + "/perturbed" +
                 to_string(this->doc_info().number()) + ".dat";
      output_stream.open(filename);
      Fluid_mesh_pt->output(output_stream, npts);
      output_stream.close();

      filename = this->doc_info().directory() + "/perturbed_free_surface" +
                 to_string(this->doc_info().number()) + ".dat";
      output_stream.open(filename);
      Free_surface_mesh_pt->output(output_stream);
      output_stream.close();

      filename = this->doc_info().directory() + "/perturbed_slip_surface" +
                 to_string(this->doc_info().number()) + ".dat";
      output_stream.open(filename);
      if (Slip_mesh_pt)
      {
        Slip_mesh_pt->output(output_stream);
      }
      output_stream.close();

      if (Contact_line_mesh_pt)
      {
        filename = this->doc_info().directory() + "/contact" +
                   to_string(this->doc_info().number()) + ".dat";
        output_stream.open(filename);
        Contact_line_mesh_pt->output(output_stream);
        output_stream.close();
      }

      if (Volume_data_pt)
      {
        filename = this->doc_info().directory() + "/perturbed_volume_trace.dat";
        output_stream.open(filename, std::ios_base::app);
        output_stream << this->doc_info().number() << " "
                      << Volume_data_pt->value(0) << std::endl;
        output_stream.close();
      }

      if (Flux_mesh_pt)
      {
        filename = this->doc_info().directory() + "/perturbed_flux_trace.dat";
        output_stream.open(filename, std::ios_base::app);
        output_stream << this->doc_info().number() << " "
                      << Flux_lagrange_multiplier_data_pt->value(0) << std::endl;
        output_stream.close();

        filename = this->doc_info().directory() + "/flux_surface" +
                   to_string(this->doc_info().number()) + ".dat";
        output_stream.open(filename);
        const unsigned n_element = Flux_mesh_pt->nelement();
        for (unsigned e = 2; e < n_element; e++)
        {
          dynamic_cast<FLUX_COMPUTATION_ELEMENT*>(Flux_mesh_pt->element_pt(e))
            ->output(output_stream);
        }
        output_stream.close();
      }

      bool Has_doc_trace = true;
      if (Has_doc_trace)
      {
        Node* left_hand_node_pt = 0;
        find_interface_end_node(false, left_hand_node_pt);
        Node* right_hand_node_pt = 0;
        find_interface_end_node(true, right_hand_node_pt);

        filename = this->doc_info().directory() + "/perturbed_trace.dat";
        output_stream.open(filename, std::ios_base::app);
        output_stream << this->doc_info().number() << " ";
        output_stream << this->time_pt()->time() << " ";
        output_stream << left_hand_node_pt->value(10) << " ";
        output_stream << right_hand_node_pt->value(10) << " ";
        output_stream << Tilde_external_pressure_data_pt->value(0) << " ";
        // output_stream << Tilde_external_pressure_data_pt->value(1) << " ";
        output_stream << std::endl;
        output_stream.close();
      }

      // Bump up counter
      this->doc_info().number()++;
    }

    double get_centre_point_ZC()
    {
      unsigned element_index = 0;
      unsigned node_index = 0;
      find_corner_bulk_element_and_node(
        Inner_boundary_id, Free_surface_boundary_id, element_index, node_index);
      const unsigned cs_index = 0;
      return dynamic_cast<FiniteElement*>(
               this->Fluid_mesh_pt->element_pt(element_index))
        ->node_pt(node_index)
        ->value(2 + cs_index);
    }

    double get_volume()
    {
      return Volume_data_pt->value(0);
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

      unsigned n_nod = Fluid_mesh_pt->nboundary_node(Free_surface_boundary_id);
      for (unsigned inod = 0; inod < n_nod; inod++)
      {
        // Get boundary node
        Node* nod_pt =
          Fluid_mesh_pt->boundary_node_pt(Free_surface_boundary_id, inod);

        if (nod_pt->is_on_boundary(boundary_id))
        {
          node_pt = nod_pt;

          return;
        }
      }
    }

    void doc_full_eigenmatrices()
    {
      oomph_info << "doc_eigenmatrices" << std::endl;
      DoubleVector residuals;
      CRDoubleMatrix jacobianOnly;
      get_jacobian(residuals, jacobianOnly);
      std::ofstream output_stream(this->doc_info().directory() +
                                  "/jacobianOnly.dat");
      output_stream.precision(16);
      output_matrices(output_stream, jacobianOnly);
      output_stream.close();

      CRDoubleMatrix jacobian;
      CRDoubleMatrix mass_matrix;
      get_eigenproblem_matrices(mass_matrix, jacobian);
      output_stream.open(this->doc_info().directory() + "/mass_matrix.dat");
      output_matrices(output_stream, mass_matrix);
      output_stream.close();
      output_stream.open(this->doc_info().directory() + "/jacobian.dat");
      output_matrices(output_stream, jacobian);
      output_stream.close();

      output_stream.open(this->doc_info().directory() + "/dofs.txt");
      this->describe_dofs(output_stream);
      output_stream.close();
    }

    void doc_sparse_eigenmatrices()
    {
      oomph_info << "doc_sparse_eigenmatrices" << std::endl;
      DoubleVector residuals;
      CRDoubleMatrix jacobianOnly;
      get_jacobian(residuals, jacobianOnly);
      jacobianOnly.sparse_indexed_output(
        this->doc_info().directory() + "/jacobianOnly.dat", 16, true);

      CRDoubleMatrix jacobian;
      CRDoubleMatrix mass_matrix;
      get_eigenproblem_matrices(mass_matrix, jacobian);
      mass_matrix.sparse_indexed_output(
        this->doc_info().directory() + "/mass_matrix.dat", 16, true);
      jacobian.sparse_indexed_output(
        this->doc_info().directory() + "/jacobian.dat", 16, true);

      std::ofstream output_stream;
      output_stream.open(this->doc_info().directory() + "/dofs.txt");
      this->describe_dofs(output_stream);
      output_stream.close();
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


        // unsigned n_element = Bulk_mesh_pt->nelement();
        // for (unsigned n = 0; n < n_element; n++)
        // {
        //   Vector<double> element_residuals(ndof());
        //   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(n))
        //     ->debug_jacobian(ndof(), element_residuals, jacobian,
        //     jacobianFD);
        // }

        // unsigned n_element = Volume_computation_mesh_pt->nelement();
        // Vector<double> element_residuals(ndof());
        // for (unsigned n = 0; n < n_element; n++)
        //{
        //  dynamic_cast<
        //    DebugElasticAxisymmetricVolumeConstraintBoundingElement<ELEMENT>*>(
        //    Volume_computation_mesh_pt->element_pt(n))
        //    ->debug_jacobian(ndof(), element_residuals, jacobian,
        //    jacobianFD);
        //}

        // unsigned n_element = No_penetration_boundary_mesh_pt->nelement();
        // Vector<double> element_residuals(ndof());
        // for (unsigned n = 0; n < n_element; n++)
        //{
        //  dynamic_cast<DebugImposeImpenetrabilityElement<ELEMENT>*>(
        //    No_penetration_boundary_mesh_pt->element_pt(n))
        //    ->debug_jacobian(ndof(), element_residuals, jacobian,
        //    jacobianFD);
        //}
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
  };
} // namespace oomph

#endif
