#ifndef ZERO_STRESS_SECTOR_PROBLEM_HEADER
#define ZERO_STRESS_SECTOR_PROBLEM_HEADER

#include <cmath>

#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"
#include "meshes/triangle_mesh.h"

/// Local headers
#include "parameter_functions.h"
#include "parameter_struct.h"
#include "refined_sector_tri_mesh.template.h"
#include "refined_sector_tri_mesh.template.cc"

namespace oomph
{
  // Problem class
  template<class ELEMENT>
  class ZeroStressSectorProblem : public Problem
  {
  private:
    /// Private variables
    RefinedSectorTriMesh<ELEMENT>* Bulk_mesh_pt;

    DocInfo Doc_info;
    Z2ErrorEstimator* Z2_error_estimator_pt;

    Params My_params;

  public:
    enum
    {
      Slip_boundary_id,
      Far_field_boundary_id,
      Free_surface_boundary_id,
    };

    // Constructor
    ZeroStressSectorProblem() : Z2_error_estimator_pt(new Z2ErrorEstimator)
    {
      // Create and add the timestepper
      add_time_stepper_pt(new BDF<2>);

      // Assign doc info pointer
      Doc_info.set_directory("RESLT_zero");
      Doc_info.number() = 0;

      /// Create parameters from parameters file.
      My_params = create_parameters_from_file("parameters.dat");

      // Create an empty mesh
      add_bulk_mesh();

      build_global_mesh();
    }

    virtual void setup()
    {
      this->setup_bulk_elements();

      set_boundary_conditions();

      this->rebuild_global_mesh();
      oomph_info << "Number of unknowns: " << this->assign_eqn_numbers()
                 << std::endl;
    }


    void setup_bulk_elements()
    {
      unsigned n_bulk = Bulk_mesh_pt->nelement();
      for (unsigned e = 0; e < n_bulk; e++)
      {
        // Upcast from GeneralisedElement to the present element
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

        // Set the Reynolds number
        el_pt->re_pt() = &My_params.reynolds_number;

        // Set the Reynolds Strouhal number
        el_pt->re_st_pt() = &My_params.strouhal_reynolds_number;
      }
    }

    Params& my_parameters()
    {
      return My_params;
    }

    // Destructor
    ~ZeroStressSectorProblem()
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

    DocInfo* doc_info_pt()
    {
      return &Doc_info;
    }

    RefinedSectorTriMesh<ELEMENT>*& bulk_mesh_pt()
    {
      return Bulk_mesh_pt;
    }


    void add_bulk_mesh();

  public:
    void set_boundary_conditions()
    {
      oomph_info << "set_boundary_conditions" << endl;

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
          const unsigned pressure_index = 2;
          node_pt->pin(pressure_index);
          node_pt->set_value(pressure_index, 0.0);
          oomph_info << node_pt->x(0) << ", " << node_pt->x(1) << " " << endl;
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
    }

  private:

    void find_corner_bulk_node(const unsigned& boundary_1_id,
                               const unsigned& boundary_2_id,
                               unsigned& element_index,
                               unsigned& node_index);
  };

  template<class ELEMENT>
  void ZeroStressSectorProblem<ELEMENT>::add_bulk_mesh()
  {
    oomph_info << "add_bulk_mesh" << endl;

    // Generate the mesh using the template ELEMENT
    Bulk_mesh_pt = new RefinedSectorTriMesh<ELEMENT>(
      My_params.n_radial,
      My_params.geometric_base,
      My_params.n_azimuthal,
      My_params.sector_radius,
      My_params.sector_angle * MathematicalConstants::Pi / 180.0,
      this->time_stepper_pt());

    // Add mesh to problem
    add_sub_mesh(Bulk_mesh_pt);
  }

  template<class ELEMENT>
  void ZeroStressSectorProblem<ELEMENT>::find_corner_bulk_node(
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
  void ZeroStressSectorProblem<ELEMENT>::doc_solution()
  {
    unsigned doc_number = Doc_info.number();

    oomph_info << "Doc Number: " << doc_number << endl;

    std::ofstream output_stream;
    char filename[100];
    unsigned npts = 3;

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

    Doc_info.number()++;
  }
} // namespace oomph
#endif
