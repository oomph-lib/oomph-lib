#ifndef CONSTRAINED_VOLUME_ELEMENTS_HEADER
#define CONSTRAINED_VOLUME_ELEMENTS_HEADER

#include "fluid_interface.h"

namespace oomph
{
  class VolumeConstraintElementWithOutput : public VolumeConstraintElement
  {
  protected:
    double* Upgraded_prescribed_volume_pt;

  public:
    VolumeConstraintElementWithOutput(double* prescribed_volume_pt)
      : VolumeConstraintElement(prescribed_volume_pt),
        Upgraded_prescribed_volume_pt(prescribed_volume_pt)
    {
    }

    VolumeConstraintElementWithOutput(double* prescribed_volume_pt,
                                      Data* p_traded_data_pt,
                                      const unsigned& index_of_traded_pressure)
      : VolumeConstraintElement(
          prescribed_volume_pt, p_traded_data_pt, index_of_traded_pressure),
        Upgraded_prescribed_volume_pt(prescribed_volume_pt)
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
      GeneralisedElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    // Add output function
    void output(std::ostream& outfile)
    {
      outfile << *this->Upgraded_prescribed_volume_pt << " ";
      outfile << this->p_traded() << std::endl;
    }
  };

} // namespace oomph

#endif
