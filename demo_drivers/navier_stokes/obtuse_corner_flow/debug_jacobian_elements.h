#ifndef DEBUG_JACOBIAN_ELEMENTS_HEADER
#define DEBUG_JACOBIAN_ELEMENTS_HEADER

#include "generic/elements.h"

namespace oomph
{
  class DebugJacobianFiniteElement : public virtual FiniteElement
  {
  public:
    unsigned debug_jacobian(const unsigned& ndof,
                            Vector<double>& residuals,
                            DenseMatrix<double>& jacobian,
                            DenseMatrix<double>& jacobianFD)
    {
      // Zero the residuals vector
      residuals.initialise(0.0);
      // Zero the jacobian matrix
      jacobian.initialise(0.0);
      // Zero the jacobian matrix
      jacobianFD.initialise(0.0);

      this->fill_in_contribution_to_jacobian(residuals, jacobian);

      residuals.initialise(0.0);

      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobianFD);

      double abs_err;
      double rel_err;
      bool jacobians_are_different = false;
      for (unsigned i = 0; i < ndof; i++)
      {
        for (unsigned j = 0; j < ndof; j++)
        {
          abs_err = abs(jacobian(j, i) - jacobianFD(j, i));
          rel_err = abs_err / abs(jacobianFD(j, i));
          if (rel_err > 1e-4 && abs_err > 2e-6)
          {
            jacobians_are_different = true;
          }
        }
      }

      if (jacobians_are_different)
      {
        oomph_info << "WARNING: Computed Jacobian is different to the finite "
                      "differenced Jacobian"
                   << std::endl;
      }
      return jacobians_are_different;
    }
  };

  class DebugJacobianSolidFiniteElement : public virtual SolidFiniteElement
  {
  public:
    unsigned debug_jacobian(const unsigned& ndof,
                            Vector<double>& residuals,
                            DenseMatrix<double>& jacobian,
                            DenseMatrix<double>& jacobianFD)
    {
      // Zero the residuals vector
      residuals.initialise(0.0);
      // Zero the jacobian matrix
      jacobian.initialise(0.0);
      // Zero the jacobian matrix
      jacobianFD.initialise(0.0);

      this->fill_in_contribution_to_jacobian(residuals, jacobian);

      residuals.initialise(0.0);

      SolidFiniteElement::fill_in_contribution_to_jacobian(residuals,
                                                           jacobianFD);

      double abs_err;
      double rel_err;
      bool jacobians_are_different = false;
      for (unsigned i = 0; i < ndof; i++)
      {
        for (unsigned j = 0; j < ndof; j++)
        {
          abs_err = abs(jacobian(j, i) - jacobianFD(j, i));
          rel_err = abs_err / abs(jacobianFD(j, i));
          if (rel_err > 1e-4 && abs_err > 2e-6)
          {
            jacobians_are_different = true;
          }
        }
      }

      if (jacobians_are_different)
      {
        oomph_info << "WARNING: Computed Jacobian is different to the finite "
                      "differenced Jacobian"
                   << std::endl;
      }
      return jacobians_are_different;
    }
  };
}; // namespace oomph
#endif
