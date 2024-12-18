#ifndef SINGULAR_FAR_FIELD_ELEMENT_HEADER
#define SINGULAR_FAR_FIELD_ELEMENT_HEADER


namespace oomph
{
  // My free surface element, overloads the sigma function to remove the
  // surface tension term
  template<class ELEMENT>
  class SingularFarFieldElement : public virtual NavierStokesFaceElement,
                                  public virtual FaceGeometry<ELEMENT>
  {
  private:
    std::function<Vector<Vector<double>>(const Vector<double>&)>
      grad_velocity_singular_fct;

    /// Lagrange multiplier id
    const unsigned Lagrange_id;

  public:
    SingularFarFieldElement(FiniteElement* const& element_pt,
                            const int& face_index,
                            Data* const& singular_scaling_data_pt,
                            const unsigned& id = 0)
      : NavierStokesFaceElement(), FaceGeometry<ELEMENT>(), Lagrange_id(id)
    {
      // Attach the geometrical information to the element
      // This function also assigned nbulk_value from required_nvalue of the
      // bulk element
      element_pt->build_face_element(face_index, this);

      this->add_external_data(singular_scaling_data_pt);
    }

    void set_grad_velocity_fct(
      std::function<Vector<Vector<double>>(const Vector<double>&)>
        Grad_velocity_singular_function)
    {
      grad_velocity_singular_fct = Grad_velocity_singular_function;
    }

    /// Fill in the residuals
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic routine with the flag set to 0
      fill_in_generic_contribution_to_residuals(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

  private:
    /// Return the index at which the lagrange multiplier is
    /// stored at the n-th node
    inline unsigned lagrange_index(const unsigned& n, const unsigned& i)
    {
      return dynamic_cast<BoundaryNodeBase*>(this->node_pt(n))
               ->index_of_first_value_assigned_by_face_element(Lagrange_id) +
             i;
    }

    /// Helper function to compute the residuals and, if flag==1, the
    /// Jacobian
    void fill_in_generic_contribution_to_residuals(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      // Get the pointer to the bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

      // Find out how many nodes there are
      unsigned n_node = this->nnode();
      unsigned n_bulk_node = bulk_el_pt->nnode();

      // Number of dimensions of the nodes
      unsigned nodal_dim = this->nodal_dimension();

      // Set up memory for the shape and test functions
      // Lagrange multiplier shape functions
      Shape psil(n_node);
      // Fluid shape functions
      Shape psif(n_bulk_node);
      DShape dpsifdx(n_bulk_node, nodal_dim);

      // to store normal vector
      Vector<double> norm_vec(nodal_dim);

      double scaling = this->external_data_pt(0)->value(0);

      // Loop over the all the integration points
      unsigned n_intpt = this->integral_pt()->nweight();
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Set up the shape elements
        // Lagrange multiplier shape elements
        this->shape_at_knot(ipt, psil);

        // Get the local coordinates of the integral knot
        Vector<double> s(nodal_dim - 1);
        for (unsigned i = 0; i < nodal_dim - 1; i++)
        {
          s[i] = this->integral_pt()->knot(ipt, i);
        }

        // Get the coordinates in the bulk element
        Vector<double> s_bulk(nodal_dim);
        this->get_local_coordinate_in_bulk(s, s_bulk);

        // Calculate the fluid shape functions
        double J = bulk_el_pt->dshape_eulerian(s_bulk, psif, dpsifdx);

        // Get the integral weight
        double w = this->integral_pt()->weight(ipt);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        Vector<double> x(nodal_dim);
        this->interpolated_x(s, x);
        Vector<Vector<double>> interpolated_dudx =
          grad_velocity_singular_fct(x);

        // Get the unit normal
        this->outer_unit_normal(ipt, norm_vec);

        // Assemble residuals and jacobian

        // loop over the lagrange multiplier components
        for (unsigned l = 0; l < nodal_dim; l++)
        {
          // Loop over the nodes
          for (unsigned j = 0; j < n_node; j++)
          {
            // Local eqn number for the l-th component of lamdba
            // in the j-th element
            int local_eqn = this->nodal_local_eqn(j, lagrange_index(j, l));

            if (local_eqn >= 0)
            {
              // Loop over coordinates
              for (unsigned i = 0; i < nodal_dim; i++)
              {
                // Assemble residual for lagrange multiplier
                residuals[local_eqn] +=
                  scaling * psil(j) * norm_vec[i] * interpolated_dudx[l][i] * W;
              }
            }
          }
        }
      }
    }
  };

} // namespace oomph
#endif
