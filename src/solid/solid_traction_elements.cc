#include "solid_traction_elements.h"

namespace oomph
{
  //=======================================================================
  /// Namespace containing the zero traction function for solid traction
  /// elements
  //=======================================================================
  namespace SolidTractionElementHelper
  {
    //=======================================================================
    /// Default load function (zero traction)
    //=======================================================================
    void Zero_traction_fct(const Vector<double>& xi,
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
  } // namespace SolidTractionElementHelper
} // namespace oomph
