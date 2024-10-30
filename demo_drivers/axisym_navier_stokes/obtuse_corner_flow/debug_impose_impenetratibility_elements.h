#ifndef DEBUG_IMPOSE_IMPENETRATIBILITY_ELEMENTS_HEADER
#define DEBUG_IMPOSE_IMPENETRATIBILITY_ELEMENTS_HEADER

#include "solid/solid_elements.h"
#include "navier_stokes.h"
#include "debug_jacobian_elements.h"

namespace oomph
{
  template<class ELEMENT>
  class DebugImposeImpenetrabilityElement
    : public virtual ImposeImpenetrabilityElement<ELEMENT>,
      public virtual DebugJacobianSolidFiniteElement,
      public virtual SolidFaceElement
  {
  public:
    DebugImposeImpenetrabilityElement(FiniteElement* const& element_pt,
                                      const int& face_index,
                                      const unsigned& id = 0)
      : ImposeImpenetrabilityElement<ELEMENT>(element_pt, face_index, id),
        DebugJacobianSolidFiniteElement(),
        SolidFaceElement()
    {
      this->add_other_bulk_node_positions_as_external_data();
    }

    /// The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default
    /// This final over-ride is required because both SolidFiniteElements
    /// and FaceElements overload zeta_nodal
    virtual double zeta_nodal(const unsigned& n,
                              const unsigned& k,
                              const unsigned& i) const
    {
      return FaceElement::zeta_nodal(n, k, i);
    }

    /// Set pointer to MacroElement -- overloads generic version
    /// and uses the MacroElement
    /// also as the default for the "undeformed" configuration.
    /// This assignment must be overwritten with
    /// set_undeformed_macro_elem_pt(...) if the deformation of
    /// the solid body is driven by a deformation of the
    /// "current" Domain/MacroElement representation of it's boundary.
    /// Can be overloaded in derived classes to perform additional
    /// tasks
    virtual void set_macro_elem_pt(MacroElement* macro_elem_pt)
    {
      Macro_elem_pt = macro_elem_pt;
      Undeformed_macro_elem_pt = macro_elem_pt;
    }

    /// Fill in contribution from Jacobian
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // // Call the generic routine with the flag set to 1
      // this->fill_in_generic_contribution_to_residuals_parall_lagr_multiplier(
      //   residuals, jacobian, 1);

      // // Get the solid entries in the jacobian using finite differences
      // fill_in_jacobian_from_solid_position_by_fd(residuals, jacobian);

      // FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
      SolidFiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }


    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    /// Add the elemental contribution to the derivatives of
    /// the residuals with respect to a parameter.
    void fill_in_contribution_to_dresiduals_dparameter(
      double* const& parameter_pt, Vector<double>& dres_dparam)
    {
    }

    /// Compute the element's residual Vector and the jacobian matrix
    /// Virtual function can be overloaded by hanging-node version
    void fill_in_contribution_to_djacobian_dparameter(
      double* const& parameter_pt,
      Vector<double>& dres_dparam,
      DenseMatrix<double>& djac_dparam)
    {
    }

    void output(std::ostream& outfile)
    {
      unsigned nplot = 5;
      output(outfile, nplot);
    }

    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      // Number of dimensions
      unsigned n_dim = this->nodal_dimension();

      // Find out how many nodes there are
      const unsigned n_node = this->nnode();

      // Set up memory for the shape functions
      Shape psi(n_node);

      // Local and global coordinates
      Vector<double> s(n_dim - 1);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(n_plot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, n_plot, s);

        // Find the shape functions
        this->shape(s, psi);

        // Initialise to zero
        Vector<double> interpolated_x(n_dim);
        Vector<double> interpolated_u(n_dim + 1);
        double interpolated_lambda = 0.0;

        // Calculate stuff
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned i = 0; i < n_dim; i++)
          {
            interpolated_x[i] += this->nodal_position(l, i) * psi[l];
          }

          for (unsigned i = 0; i < n_dim + 1; i++)
          {
            interpolated_u[i] +=
              this->nodal_value(l, this->u_index_nst(l, i)) * psi[l];
          }
          interpolated_lambda +=
            this->nodal_value(l, this->lagrange_multiplier_index(l)) * psi[l];
        }

        // Output x,y,..
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << interpolated_x[i] << " ";
        }

        // Output u,v,w
        for (unsigned i = 0; i < n_dim + 1; i++)
        {
          outfile << interpolated_u[i] << " ";
        }

        // Output the lagrange multiplier
        outfile << interpolated_lambda << std::endl;
      }
    }
  };
} // namespace oomph

#endif
