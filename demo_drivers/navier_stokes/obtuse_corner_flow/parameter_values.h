#ifndef PARAMETER_VALUES_HEADER
#define PARAMETER_VALUES_HEADER

#include "generic.h"

namespace oomph
{
  namespace parameters
  {
    double contact_angle = 2.0 / 3.0 * MathematicalConstants::Pi;
    const double slip_length = 1e0;
    double capillary_number = 1e0;
    const double element_area = 1e-2;
    double nu = 0.25;
    double reynolds_number = 0.0;
    double strouhal_reynolds_number = 1.0;
    const double small_r = 1e-4;
    const double inner_radius = 0.001;
    const unsigned nsegment = 32;
  } // namespace parameters
} // namespace oomph
#endif
