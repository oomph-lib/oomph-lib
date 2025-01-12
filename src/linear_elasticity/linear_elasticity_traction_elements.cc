#include "linear_elasticity_traction_elements.h"

namespace oomph
{
  //=======================================================================
  /// Namespace containing the zero traction function for linear elasticity
  /// traction elements
  //=======================================================================
  namespace LinearElasticityTractionElementHelper
  {
    void Zero_traction_fct(const double& time,
                           const Vector<double>& x,
                           const Vector<double>& N,
                           Vector<double>& load)
    {
      unsigned n_dim = load.size();
      for (unsigned i = 0; i < n_dim; i++)
      {
        load[i] = 0.0;
      }
    }
  } // namespace LinearElasticityTractionElementHelper
} // namespace oomph
