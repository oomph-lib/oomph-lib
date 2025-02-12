#ifndef OOMPH_EIGENSOLUTION_FUNCTIONS_HEADER
#define OOMPH_EIGENSOLUTION_FUNCTIONS_HEADER

// Std includes
#include <functional>

// Generic includes
#include "../generic/Vector.h"
#include "../generic/nodes.h"

namespace oomph
{
  /// function to return radial distance from a point
  double get_radial_distance(const Vector<double>& x,
                             const Vector<double>& x_centre);

  /// function to return angle of radial line from positive y-axis
  double get_azimuthal_angle(const Vector<double>& x,
                             const Vector<double>& x_centre);

  /// function to convert velocity vectors from polar to Cartesian coords
  void polar_to_cartesian_velocity(const Vector<double>& u_polar,
                                   const double& theta,
                                   Vector<double>& u_cart);

  /// Function that computes the singular pressure solution
  double pressure_singular_fct(const Vector<double>& x);

  /// Function that computes the gradient of the singular pressure solution
  Vector<double> grad_pressure_singular_fct(const Vector<double>& coord);

  /// Function that return the singular velocity function
  std::function<Vector<double>(const Vector<double>&)> velocity_singular_function_factory(
    const double& contact_angle, Node*& x_centre_node_pt);

  /// Function that returns the gradient of the singular velocity function with
  /// the appropriate parameters
  std::function<Vector<Vector<double>>(const Vector<double>&)> grad_velocity_singular_function_factory(
    const double& contact_angle, Node*& x_centre_node_pt);

  /// Function that returns an eigensolution traction function with the
  /// appropriate parameters.
  std::function<void(const double&,
                     const Vector<double>&,
                     const Vector<double>&,
                     Vector<double>&)>
  eigensolution_traction_function_factory(
    const double& contact_angle,
    const std::function<Vector<Vector<double>>(const Vector<double>&)>&
      grad_velocity_singular_fct);

  /// Function that returns an eigensolution slip function with the
  /// appropriate parameters.
  std::function<void(const double&,
                     const Vector<double>&,
                     const Vector<double>&,
                     Vector<double>&)>
  eigensolution_slip_function_factory(
    const double& slip_length,
    const std::function<Vector<double>(const Vector<double>&)>&
      velocity_singular_fct);

} // namespace oomph
#endif
