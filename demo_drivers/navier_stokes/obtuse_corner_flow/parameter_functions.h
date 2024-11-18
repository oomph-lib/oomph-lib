#ifndef PARAMETER_FUNCTIONS_HEADER
#define PARAMETER_FUNCTIONS_HEADER

#include "generic.h"
#include "parameter_values.h"
#include "parameter_struct.h"

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
    wall_velocity[1] = -1.0;
  }

  namespace parameters
  {
    // Slip length near the on the outer wall near the contact line
    void prescribed_slip_fct(const double& t,
                             const Vector<double>& x,
                             const Vector<double>& n,
                             Vector<double>& slip)
    {
      const double local_slip = parameters::slip_length;

      // Assign solution
      slip[0] = -1;
      slip[1] = local_slip;
    }

    // Wall velocity
    void prescribed_wall_velocity_fct(const double& t,
                                      const Vector<double>& x,
                                      const Vector<double>& n,
                                      Vector<double>& wall_velocity)
    {
      // Assign solution
      wall_velocity[0] = 0.0;
      wall_velocity[1] = -1.0;
    }

    void contact_angle_fct_pt(const oomph::Vector<double>& local_u,
                              double& contact_angle)
    {
      contact_angle = parameters::contact_angle;
    }

    // Function that specifies the wall unit normal, into the fluid
    void wall_unit_normal_fct(const oomph::Vector<double>& x,
                              oomph::Vector<double>& normal)
    {
      normal.resize(2);
      normal[0] = -1.0;
      normal[1] = 0.0;
    }
  } // namespace parameters
} // namespace oomph
#endif
