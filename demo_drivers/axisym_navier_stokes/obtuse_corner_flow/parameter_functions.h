#ifndef PARAMETER_FUNCTIONS_HEADER
#define PARAMETER_FUNCTIONS_HEADER

#include "generic.h"

namespace oomph
{
  std::function<void(const double&,
                     const Vector<double>&,
                     const Vector<double>&,
                     Vector<double>&)>
  slip_function_factory(const double& slip_length)
  {
    // Return a lambda that captures 's' and takes 't' and 'x' as arguments
    return [slip_length](const double& t,
                         const Vector<double>& x,
                         const Vector<double>& n,
                         Vector<double>& slip) -> void
    {
      slip[0] = -1;
      slip[1] = slip_length;
      slip[2] = -1;
    };
  }

  // Wall velocity
  void wall_velocity_function(const double& t,
                              const Vector<double>& x,
                              const Vector<double>& n,
                              Vector<double>& wall_velocity)
  {
    // Assign solution
    wall_velocity[0] = 0.0;
    wall_velocity[1] = 1.0;
    wall_velocity[2] = 0.0;
  }
} // namespace oomph
#endif
