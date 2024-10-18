#ifndef PRESSURE_EVALUATION_ELEMENTS_HEADER
#define PRESSURE_EVALUATION_ELEMENTS_HEADER

#include "generic.h"
#include "debug_jacobian_elements.h"

namespace oomph
{
  //==================CLASS FOR THE PRESSURE CONTRIBUTION================
  /// This class adds the finite element pressure at the evaluation point
  /// to the residual for the singular eigensolution function.
  ///
  /// R_C += +- p_FE (Evaluation_point)
  ///
  /// and thus regularises the FE solution by matching the pressure at
  /// two locations. If the amplitude of the  singular solution is known,
  /// pin C.
  //=====================================================================
  template<class ELEMENT>
  class PressureEvaluationElement
    : public virtual FaceGeometry<ELEMENT>,
      public virtual SolidFaceElement,
      public virtual DebugJacobianSolidFiniteElement

  {
  private:
    // Storage for the bulk element
    ELEMENT* Cast_bulk_element_pt;
    int Pressure_index;
    bool Is_adding_to_residuals;
    Vector<int> Node_index;
    double evalution_point_r;
    Vector<double> evalution_point_s;
    Node* corner_node_pt;


  public:
    // Constructor
    PressureEvaluationElement(FiniteElement* const& element_pt,
                              const int& face_index,
                              Node* const& node_pt)
      : FaceGeometry<ELEMENT>(),
        SolidFaceElement(),
        Cast_bulk_element_pt(dynamic_cast<ELEMENT*>(element_pt)),
        Pressure_index(-1),
        Is_adding_to_residuals(true),
        evalution_point_r(1e-3),
        evalution_point_s(1, 0.0),
        corner_node_pt(node_pt)
    {
      // Attach the geometrical information to the element. N.B. This function
      // also assigns nbulk_value from the required_nvalue of the bulk element
      element_pt->build_face_element(face_index, this);

      // Add the nodes (which are data) where the pressure is stored in the bulk
      // element as external data.
      // add_pressure_nodes_not_on_face_as_external_data();
      this->add_other_bulk_node_positions_as_external_data();
    }

    void compute_s(Node* const& node_pt)
    {
      // Find which node is the corner node
      unsigned corner_node_index = 0;
      for (unsigned n = 0; n < this->nnode(); n++)
      {
        if (node_pt == this->node_pt(n))
        {
          corner_node_index = n;
        }
      }

      // Set the evaluation point s
      if (corner_node_index == 0)
      {
        evalution_point_s[0] = evalution_point_r / this->size();
      }
      else
      {
        evalution_point_s[0] = 1 - evalution_point_r / this->size();
      }
    }

  private:
    void add_pressure_nodes_not_on_face_as_external_data()
    {
      // Find number of nodes where the pressure is stored in the element
      unsigned nnod = Cast_bulk_element_pt->npres_nst();

      // Loop over the nodes of the bulk element
      for (unsigned j = 0; j < nnod; j++)
      {
        // Add the node as external data in the
        // SingularAxisymmetricNavierStokesSolutionElement class. Note that this
        // assumes that the pressure is stored at the nodes (Taylor Hood type
        // NSt elements, which is assumed elsewhere too...)
        // If the node is a pressure node and not in the element, ...
        if (Cast_bulk_element_pt->is_pressure_node(j) &&
            this->get_node_number(Cast_bulk_element_pt->node_pt(j)) == -1)
        {
          // ... add it as external data.
          Node_index.push_back(
            add_external_data(Cast_bulk_element_pt->node_pt(j)));
        }
      }
    }

  public:
    void set_add_to_residuals()
    {
      Is_adding_to_residuals = true;
    }
    void set_subtract_from_residuals()
    {
      Is_adding_to_residuals = false;
    }

    // Set and add the pressure data as external data
    void set_pressure_data_pt(Data* const& pressure_data_pt)
    {
      Pressure_index = add_external_data(pressure_data_pt);
    }

    // Calculate the element's residual vector
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_pressure_contribution(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    // Calculate the element's residual vector and Jacobian
    // void fill_in_contribution_to_jacobian(Vector<double>& residuals,
    //                                       DenseMatrix<double>& jacobian)
    // {
    //   // Call the generic routine with the flag set to 1
    //   fill_in_generic_residual_contribution_pressure_contribution(
    //     residuals, jacobian, 1);

    //   // Call the generic finite difference routine to handle the solid
    //   // variables
    //   this->fill_in_jacobian_from_solid_position_by_fd(jacobian);
    // }

    void fill_in_contribution_to_dresiduals_dparameter(
      double* const& parameter_pt, Vector<double>& dres_dparam)
    {
    }

    /// Specify the value of nodal zeta from the face geometry
    /// The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default (needed to break
    /// indeterminacy if bulk element is SolidElement)
    double zeta_nodal(const unsigned& n,
                      const unsigned& k,
                      const unsigned& i) const
    {
      return FaceElement::zeta_nodal(n, k, i);
    }

  protected:
    // Generic residual and Jacobian routine
    void fill_in_generic_residual_contribution_pressure_contribution(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
    {
      // Find out the number of surface coordinates
      const unsigned el_dim = this->dim();

      // Throw error if we have more than one dimension
      if (el_dim != 1)
      {
        std::stringstream error_stream;
        error_stream
          << "Currently, the PressureEvaluationElement elements\n"
          << "only work for 1D elements, but we have detected\n"
          << "more than 1 dimension. "
          << "If your problem has more dimensions, you're welcome to \n"
          << "volunteer and implement the required functionality in \n\n"
          << "   "
             "PressureEvaluationElement::fill_in_generic_residual_"
             "contribution_pressure_contribution()"
             "\n\n"
          << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Initialise the index at which pressure value is stored
      const unsigned pressure_value_index = 0;

      // Set the residual multiplier, dependent on whether we are adding or
      // subtracting to the residuals
      double multiplier = 1.0;
      if (!Is_adding_to_residuals)
      {
        multiplier = -1.0;
      }

      // Find out how many pressure dofs there are in the bulk element
      unsigned n_pres = Cast_bulk_element_pt->npres_nst();

      // Set up memory for pressure shape and test functions
      Shape psip(n_pres);

      // Find the coordinate in the bulk element
      compute_s(corner_node_pt);
      Vector<double> s_bulk(dim() + 1);
      s_bulk = local_coordinate_in_bulk(evalution_point_s);

      // Evaluate the pressure shape functions at the coordinate in the bulk
      // element
      Cast_bulk_element_pt->pshape_axi_nst(s_bulk, psip);

      // Set the local equation
      int local_eqn = 0;

      // Add to singular function scaling residual
      local_eqn =
        this->external_local_eqn(Pressure_index, pressure_value_index);

      // If the equation is not pinned
      if (local_eqn >= 0)
      {
        // Add (or subtract) the pressure at the evaluation point
        residuals[local_eqn] +=
          Cast_bulk_element_pt->interpolated_p_nst_fe_only(s_bulk) * multiplier;

        // If the Jacobian flag is on, add to the Jacobian
        if (flag)
        {
          // Initialise a variable for the local_unknown
          int local_unknown = 0;

          // Loop over shape functions
          const unsigned n_local_pres = Node_index.size();
          for (unsigned j = 0; j < n_local_pres; j++)
          {
            // The residual depends on the pressure at each of the bulk
            // elements nodes, which are stored here as external data.
            local_unknown = this->external_local_eqn(
              Node_index[j], Cast_bulk_element_pt->p_nodal_index_nst());

            // If not pinned
            if (local_unknown > 0)
            {
              // Add the contribution of the node to the local jacobian
              jacobian(local_eqn, local_unknown) += psip(j) * multiplier;
            }
          }
        }
      }
    }

    // Return the value of the FE pressure at the local coordinate s.
    double interpolated_p(const Vector<double>& s) const
    {
      // Initialised value of the pressure
      // Local coordinates in bulk element
      Vector<double> s_bulk(dim() + 1);
      s_bulk = local_coordinate_in_bulk(s);

      return Cast_bulk_element_pt->interpolated_p_nst_fe_only(s_bulk);
    }

    // Overwrite the output function
    void output(std::ostream& outfile)
    {
      // Vector of local coordinates
      unsigned n_dim = dim();

      // Spatial coordinates are one higher
      for (unsigned i = 0; i < n_dim + 1; i++)
      {
        outfile << interpolated_x(evalution_point_s, i) << ",";
      }

      // Output the pressure
      outfile << interpolated_p(evalution_point_s) << ",";

      // End of line
      outfile << std::endl;
    }
  };
} // namespace oomph
#endif
