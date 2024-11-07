#ifndef FAR_FIELD_ELEMENT_HEADER
#define FAR_FIELD_ELEMENT_HEADER

#include "generic.h"
#include "fluid_interface.h"
#include <algorithm>

namespace oomph
{
  // My free surface element, overloads the sigma function to remove the
  // surface tension term
  template<class ELEMENT>
  class FarFieldElement : public virtual NavierStokesFaceElement,
                          public virtual FaceGeometry<ELEMENT>
  {
  private:
    /// Lagrange multiplier id
    const unsigned Lagrange_id;

    virtual double nst_u(const unsigned& j, const unsigned& i)
    {
      return nodal_value(j, u_index_nst(j, i));
    }

  public:
    FarFieldElement(FiniteElement* const& element_pt,
                    const int& face_index,
                    const unsigned& id = 0)
      : NavierStokesFaceElement(), FaceGeometry<ELEMENT>(), Lagrange_id(id)
    {
      // Attach the geometrical information to the element
      // This function also assigned nbulk_value from required_nvalue of the
      // bulk element
      element_pt->build_face_element(face_index, this);

      // Add storage for the Lagrange multipliers
      Vector<unsigned> n_additional_values(nnode(), bulk_element_pt()->dim());

      // Now add storage for Lagrange multipliers and set the map containing
      // the position of the first entry of this face element's
      // additional values.
      this->add_additional_values(n_additional_values, Lagrange_id);

      // Add the other nodes as external data as we are imposing a constraint
      // on dudx
      this->add_other_bulk_nodes_as_external_data();
    }

    /// Return the index at which the lagrange multiplier is
    /// stored at the n-th node
    inline unsigned lagrange_index(const unsigned& n, const unsigned& i)
    {
      return dynamic_cast<BoundaryNodeBase*>(this->node_pt(n))
               ->index_of_first_value_assigned_by_face_element(Lagrange_id) +
             i;
    }

    // Fix the lagrange multiplier at node i_node in the given direction
    void pin_lagrange_multiplier(const unsigned& i_node,
                                 const unsigned& direction)
    {
      const unsigned value_index = lagrange_index(i_node, direction);
      this->node_pt(i_node)->pin(value_index);
    }

    // Free the lagrange multiplier at node i_node in the given direction
    void unpin_lagrange_multiplier(const unsigned& i_node,
                                   const unsigned& direction)
    {
      const unsigned value_index = lagrange_index(i_node, direction);
      this->node_pt(i_node)->unpin(value_index);
    }

    // Return the lagrange multiplier at node i_node in the given direction
    double lagrange_multiplier(const unsigned& i_node,
                               const unsigned& direction)
    {
      return this->node_pt(i_node)->value(lagrange_index(i_node, direction));
    }

    // Equation number of the Lagrange multiplier associated with node i_node
    // in the given direction.
    inline int lagrange_local_eqn(const unsigned& i_node,
                                  const unsigned& direction)
    {
      // Get the index of the nodal value associated with Lagrange multiplier
      return this->nodal_local_eqn(i_node, lagrange_index(i_node, direction));
    }

    /// Fill in the residuals
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic routine with the flag set to 0
      fill_in_generic_contribution_to_residuals(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    // Fill in contribution from Jacobian
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
      // Call the generic routine with the flag set to 1
      // fill_in_generic_contribution_to_residuals(residuals, jacobian, 1);

      // fill_in_jacobian_from_internal_by_fd(residuals, jacobian, true);
      // fill_in_jacobian_from_external_by_fd(residuals, jacobian, true);
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


    /// Overload the output function
    void output(std::ostream& outfile)
    {
      const unsigned n_plot = 5;
      output(outfile, n_plot);
    }

    /// Output function: x,y,[z],u,v,[w],p in csv format
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      // Get the pointer to the bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      // Number of dimensions
      const unsigned n_dim = this->nodal_dimension();

      // Find out how many nodes there are
      const unsigned n_node = nnode();
      const unsigned n_bulk_node = bulk_el_pt->nnode();

      // Set up memory for the shape functions
      Shape psil(n_node);
      Shape psif(n_bulk_node);
      DShape dpsifdx(n_bulk_node, n_dim);

      // Local and bulk coordinates
      Vector<double> s(n_dim - 1);
      Vector<double> s_bulk(n_dim);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(n_plot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, n_plot, s);
        this->get_local_coordinate_in_bulk(s, s_bulk);

        // Outer unit normal
        Vector<double> norm_vec(n_dim);
        outer_unit_normal(s, norm_vec);

        // Find the shape functions
        shape(s, psil);
        bulk_el_pt->dshape_eulerian(s_bulk, psif, dpsifdx);

        // Initialise to zero
        Vector<double> interpolated_x(n_dim);
        Vector<double> interpolated_u(n_dim);
        // Vector<Vector<double>> interpolated_dudx(n_dim,
        //                                          Vector<double>(n_dim, 0.0));
        Vector<double> interpolated_lambda(n_dim);

        // Loop over all nodes
        // for (unsigned j = 0; j < n_bulk_node; j++)
        //{
        //  // Loop over directions
        //  for (unsigned i = 0; i < n_dim; i++)
        //  {
        //    for (unsigned k = 0; k < n_dim; k++)
        //    {
        //      interpolated_dudx[i][k] += this->nst_u(j, i) * dpsifdx(j, k);
        //    }
        //  }
        //}

        // Loop over local nodes
        for (unsigned j = 0; j < n_node; j++)
        {
          // Loop over directions
          for (unsigned i = 0; i < n_dim; i++)
          {
            interpolated_x[i] += this->nodal_position(j, i) * psil[j];
            interpolated_u[i] += this->nst_u(j, i) * psif(bulk_node_number(j));
            interpolated_lambda[i] += lagrange_multiplier(j, i) * psil[j];
          }
        }


        // Output the x,y,..
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << interpolated_x[i] << ",";
        }

        // Output normal
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << norm_vec[i] << ",";
        }

        // Output the velocity
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << interpolated_u[i] << ",";
        }

        // Output the velocity derivatives
        // for (unsigned i = 0; i < n_dim; i++)
        //{
        //  for (unsigned k = 0; k < n_dim; k++)
        //  {
        //    outfile << interpolated_dudx[i][k] << ",";
        //  }
        //}

        // Output the lagrange multiplier, don't include the comma on the
        // final value
        for (unsigned i = 0; i < n_dim - 1; i++)
        {
          outfile << interpolated_lambda[i] << ",";
        }
        outfile << interpolated_lambda[n_dim - 1] << std::endl;
      }
    }

  protected:
    /// Return FE interpolated derivatives of velocity component u[i]
    /// w.r.t spatial global coordinate direction x[j] at local coordinate s
    virtual double interpolated_dudx_bulk(const Vector<double>& s,
                                          const unsigned& i,
                                          const unsigned& j) const
    {
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      // Vector for local coordinates in bulk element
      Vector<double> s_bulk(this->nodal_dimension());

      // Get the bulk coordinates
      this->get_local_coordinate_in_bulk(s, s_bulk);

      return bulk_el_pt->interpolated_dudx_nst(s_bulk, i, j);
    }

    /// Helper function to compute the residuals and, if flag==1, the
    /// Jacobian
    void fill_in_generic_contribution_to_residuals(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      // Get the pointer to the bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      // Find out how many nodes there are
      unsigned n_node = nnode();
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

      // Loop over the all the integration points
      unsigned n_intpt = integral_pt()->nweight();
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Set up the shape elements
        // Lagrange multiplier shape elements
        shape_at_knot(ipt, psil);

        // Get the local coordinates of the integral knot
        Vector<double> s(nodal_dim - 1);
        for (unsigned i = 0; i < nodal_dim - 1; i++)
        {
          s[i] = integral_pt()->knot(ipt, i);
        }

        // Get the coordinates in the bulk element
        Vector<double> s_bulk(nodal_dim);
        this->get_local_coordinate_in_bulk(s, s_bulk);

        // Calculate the fluid shape functions
        double J = bulk_el_pt->dshape_eulerian(s_bulk, psif, dpsifdx);

        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // compute  the velocity and the Lagrange multiplier
        Vector<double> interpolated_u(nodal_dim, 0.0);
        Vector<Vector<double>> interpolated_dudx(
          nodal_dim, Vector<double>(nodal_dim, 0.0));
        Vector<double> interpolated_lambda(nodal_dim, 0.0);

        // Loop over nodes
        for (unsigned j = 0; j < n_node; j++)
        {
          // Assemble the velocity component
          for (unsigned i = 0; i < nodal_dim; i++)
          {
            interpolated_u[i] +=
              this->nst_u(j, i) * psif(this->bulk_node_number(j));
            interpolated_lambda[i] += lagrange_multiplier(j, i) * psil(j);
          }
        }

        // Get velocity gradient at the knot from the bulk element
        for (unsigned i = 0; i < nodal_dim; i++)
        {
          for (unsigned j = 0; j < nodal_dim; j++)
          {
            interpolated_dudx[i][j] = interpolated_dudx_bulk(s, i, j);
          }
        }

        // Get the unit normal
        outer_unit_normal(ipt, norm_vec);

        // Assemble residuals and jacobian

        // loop over the lagrange multiplier components
        for (unsigned l = 0; l < nodal_dim; l++)
        {
          // Loop over the nodes
          for (unsigned j = 0; j < n_node; j++)
          {
            // Local eqn number for the l-th component of lamdba
            // in the j-th element
            int local_eqn = lagrange_local_eqn(j, l);

            if (local_eqn >= 0)
            {
              // Loop over coordinates
              for (unsigned i = 0; i < nodal_dim; i++)
              {
                // Assemble residual for lagrange multiplier
                residuals[local_eqn] +=
                  psil(j) * norm_vec[i] * interpolated_dudx[l][i] * W;

                // Assemble Jacobian for Lagrange multiplier:
                if (flag == 1)
                {
                  // Loop over the nodes again for unknowns
                  for (unsigned jj = 0; jj < n_node; jj++)
                  {
                    // Local eqn number for the i-th component
                    // of the velocity in the jj-th element
                    int local_unknown = this->nst_u_local_unknown(jj, i);
                    if (local_unknown >= 0)
                    {
                      // jacobian(local_eqn, local_unknown) +=
                      //   norm_vec[i] * psi(jj) * dpsidx(j, i) * W;
                    }
                  }
                }
              }
              // Local eqn number for the i-th component of the
              // velocity in the j-th element
              local_eqn = this->nst_momentum_local_eqn(j, l);

              if (local_eqn >= 0)
              {
                // Loop over the directions
                for (unsigned i = 0; i < nodal_dim; i++)
                {
                  // Add to residual
                  residuals[local_eqn] += interpolated_lambda[l] * norm_vec[i] *
                                          dpsifdx(bulk_node_number(j), i) * W;
                }

                // // Do Jacobian too?
                // if (flag == 1)
                // {
                //   // Loop over the nodes again for unknowns
                //   for (unsigned jj = 0; jj < n_node; jj++)
                //   {
                //     // Local eqn number for the l-th component of lamdba
                //     // in the jj-th element
                //     int local_unknown = lagrange_local_eqn(jj, l);
                //     if (local_unknown >= 0)
                //     {
                //       //  jacobian(local_eqn, local_unknown) +=
                //       //    psi(jj) * norm_vec[i] * dpsidx(j, i) * W;
                //     }
                //   }
                // }
              }
            }
          }
        }
      }
    }
  };

} // namespace oomph
#endif
