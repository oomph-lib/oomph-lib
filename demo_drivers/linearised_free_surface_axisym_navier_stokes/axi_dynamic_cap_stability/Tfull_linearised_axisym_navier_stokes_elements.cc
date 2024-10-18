#include "Tfull_linearised_axisym_navier_stokes_elements.h"

namespace oomph
{
  // Set the data for the pressure conversion array
  const unsigned
    TaylorHoodFullLinearisedAxisymmetricNavierStokesEquations::Pconv[3] = {
      0, 1, 2};
  const unsigned TaylorHoodFullLinearisedAxisymmetricNavierStokesEquations::
    Initial_Nvalue[6] = {12, 12, 12, 10, 10, 10};
} // namespace oomph
