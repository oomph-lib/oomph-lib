#ifndef FREE_SURFACE_ELEMENTS_HEADER
#define FREE_SURFACE_ELEMENTS_HEADER

#include "fluid_interface.h"

namespace oomph
{
  template<class ELEMENT>
  class FreeSurfaceElement
    : public ElasticAxisymmetricFluidInterfaceElement<ELEMENT>
  {
  private:
    TimeStepper* Time_stepper_pt;

  public:
    FreeSurfaceElement(FiniteElement* const& element_pt,
                       const int& face_index,
                       TimeStepper* const& time_stepper_pt,
                       const unsigned& id = 0)
      : ElasticAxisymmetricFluidInterfaceElement<ELEMENT>(
          element_pt, face_index, id),
        Time_stepper_pt(time_stepper_pt)
    {
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
  };
} // namespace oomph

#endif // FREE_SURFACE_ELEMENTS_HEADER
