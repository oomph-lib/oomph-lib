#ifndef FREE_SURFACE_ELEMENTS_HEADER
#define FREE_SURFACE_ELEMENTS_HEADER

#include "fluid_interface.h"
#include "debug_jacobian_elements.h"

namespace oomph
{
  template<class ELEMENT>
  class FreeSurfaceElement
    : public virtual ElasticAxisymmetricFluidInterfaceElement<ELEMENT>,
      public virtual DebugJacobianSolidFiniteElement,
      public virtual SolidFaceElement

  {
  private:
    TimeStepper* Time_stepper_pt;
    double Error;

  public:
    FreeSurfaceElement(FiniteElement* const& element_pt,
                       const int& face_index,
                       TimeStepper* const& time_stepper_pt,
                       const unsigned& id = 0)
      : ElasticAxisymmetricFluidInterfaceElement<ELEMENT>(
          element_pt, face_index, id),
        SolidFaceElement(),
        Time_stepper_pt(time_stepper_pt),
        Error(0.0)
    {
      this->add_other_bulk_node_positions_as_external_data();
    }

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

    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      SolidFiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    void my_fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                             DenseMatrix<double>& jacobian)
    {
      // Add the contribution to the residuals
      this->fill_in_contribution_to_residuals(residuals);

      // Allocate storage for the full residuals (residuals of entire element)
      unsigned n_dof = this->ndof();
      Vector<double> full_residuals(n_dof);
      // Get the residuals for the entire element
      this->get_residuals(full_residuals);
      // Get the solid entries in the jacobian using finite differences
      this->fill_in_jacobian_from_solid_position_by_fd(full_residuals,
                                                       jacobian);
      // There could be internal data
      //(finite-difference the lot by default)
      this->fill_in_jacobian_from_internal_by_fd(
        full_residuals, jacobian, true);
      // There could also be external data
      //(finite-difference the lot by default)
      this->fill_in_jacobian_from_external_by_fd(
        full_residuals, jacobian, true);
      // There could also be nodal data
      this->fill_in_jacobian_from_nodal_by_fd(full_residuals, jacobian);
    }


    /// Add the elemental contribution to the jacobian and mass matrices
    /// and the residuals vector. Note that
    /// this function will NOT initialise the residuals vector or the jacobian
    /// matrix. It must be called after the residuals vector and
    /// jacobian matrix have been initialised to zero. The default
    /// is to use finite differences to calculate the jacobian
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // Add the contribution to the residuals
      this->fill_in_contribution_to_residuals(residuals);

      // Allocate storage for the full residuals (residuals of entire
      // element)
      unsigned n_dof = this->ndof();
      Vector<double> full_residuals(n_dof);
      // Get the residuals for the entire element
      this->get_residuals(full_residuals);

      // Make the timestepper steady
      Time_stepper_pt->make_steady();

      // Use the fill in contribution to jacobian
      my_fill_in_contribution_to_jacobian(full_residuals, jacobian);

      Time_stepper_pt->undo_make_steady();

      DenseMatrix<double> unsteady_jacobian(n_dof, n_dof);
      my_fill_in_contribution_to_jacobian(full_residuals, unsteady_jacobian);

      const double dt = Time_stepper_pt->time_pt()->dt();
      for (unsigned i = 0; i < n_dof; i++)
      {
        for (unsigned j = 0; j < n_dof; j++)
        {
          /// The 2/3 is due to the BDF<2> scheme.
          mass_matrix(j, i) +=
            (2.0 / 3.0) * dt * (-unsteady_jacobian(j, i) + jacobian(j, i));
        }
      }
    }

    void compute_error()
    {
      // Maths from http://www.cgafaq.info/wiki/Circle_Through_Three_Points
      double a_x = this->node_pt(0)->x(0);
      double a_y = this->node_pt(0)->x(1);
      double b_x = this->node_pt(1)->x(0);
      double b_y = this->node_pt(1)->x(1);
      double c_x = this->node_pt(2)->x(0);
      double c_y = this->node_pt(2)->x(1);

      double a = b_x - a_x;
      double b = b_y - a_y;
      double c = c_x - a_x;
      double d = c_y - a_y;

      double e = a * (a_x + b_x) + b * (a_y + b_y);
      double f = c * (a_x + c_x) + d * (a_y + c_y);

      double g = 2.0 * (a * (c_y - b_y) - b * (c_x - b_x));

      double error = 0.0;
      if (std::fabs(g) >= 1.0e-14)
      {
        double p_x = (d * e - b * f) / g;
        double p_y = (a * f - c * e) / g;

        double r = sqrt(pow((a_x - p_x), 2) + pow((a_y - p_y), 2));

        double rhalfca_x = 0.5 * (a_x - c_x);
        double rhalfca_y = 0.5 * (a_y - c_y);

        double halfca_squared = pow(rhalfca_x, 2) + pow(rhalfca_y, 2);

        double sticky_out_bit = r - sqrt(std::fabs((r * r) - halfca_squared));

        // If sticky out bit divided by distance between end nodes
        // is less than tolerance the boundary is so flat that we
        // can safely kill the node
        error = sticky_out_bit / (2.0 * sqrt(halfca_squared));
      }
      Error = error;
    }

    double get_error()
    {
      return Error;
    }

    void output(std::ostream& outfile, const unsigned& npts)
    {
      ElasticAxisymmetricFluidInterfaceElement<ELEMENT>::output(outfile, npts);

      outfile << get_error() << endl;
    }

    void output(std::ostream& outfile)
    {
      ElasticAxisymmetricFluidInterfaceElement<ELEMENT>::output(outfile);

      outfile << get_error() << endl;
    }
  };

} // namespace oomph

#endif // FREE_SURFACE_ELEMENTS_HEADER
