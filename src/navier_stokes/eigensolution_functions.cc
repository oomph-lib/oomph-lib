#include "eigensolution_functions.h"

namespace oomph
{
  /// function to return radial distance from a point
  double get_radial_distance(const Vector<double>& x,
                             const Vector<double>& x_centre)
  {
    double sum = 0;
    const unsigned DIM = 2;
    for (unsigned i = 0; i < DIM; i++)
    {
      sum += pow(x[i] - x_centre[i], 2);
    }
    return sqrt(sum);
  }

  /// function to return angle of radial line from positive y-axis
  double get_azimuthal_angle(const Vector<double>& x,
                             const Vector<double>& x_centre)
  {
    return atan2(x_centre[0] - x[0], x[1] - x_centre[1]);
  }

  /// function to convert velocity vectors from polar to Cartesian coords
  void polar_to_cartesian_velocity(const Vector<double>& u_polar,
                                   const double& theta,
                                   Vector<double>& u_cart)
  {
    u_cart[0] = -u_polar[0] * sin(theta) - u_polar[1] * cos(theta);
    u_cart[1] = u_polar[0] * cos(theta) - u_polar[1] * sin(theta);
  }

  /// Function that computes the singular pressure solution
  double pressure_singular_fct(const Vector<double>& x)
  {
    double p = 0.0;

    return p;
  }

  /// Function that computes the gradient of the singular pressure solution
  Vector<double> grad_pressure_singular_fct(const Vector<double>& coord)
  {
    // Initialise the gradient vector
    Vector<double> grad_p(2, 0.0);

    return grad_p;
  }

  /// Function that return the singular velocity function
  std::function<Vector<double>(const Vector<double>&)> velocity_singular_function_factory(
    const double& contact_angle, Node*& x_centre_node_pt)
  {
    /// Function that computes the singular velocity solution near the corner at
    /// x_centre_node_pt for a given contact angle
    return [contact_angle,
            &x_centre_node_pt](const Vector<double>& x) -> Vector<double>
    {
      // Initialise velocity vector to return
      Vector<double> u(3, 0.0);

      // Find radial distance of polar coordinates centered at the point
      // x_centre
      Vector<double> x_centre = x_centre_node_pt->position();
      double r = get_radial_distance(x, x_centre);
      // Find angle of polar coordinates centered at the point x_centre
      double theta = get_azimuthal_angle(x, x_centre);

      // Compute lambda, the ratio of
      double lambda = MathematicalConstants::Pi / contact_angle;

      u[0] = lambda * pow(r, lambda - 1) * sin((lambda - 1) * theta);
      u[1] = lambda * pow(r, lambda - 1) * cos((lambda - 1) * theta);
      u[2] = 0;

      // return the velocity vector
      return u;
    }; // End of lambda function
  }

  /// Function that returns the gradient of the singular velocity function with
  /// the appropriate parameters
  std::function<Vector<Vector<double>>(const Vector<double>&)> grad_velocity_singular_function_factory(
    const double& contact_angle, Node*& x_centre_node_pt)
  {
    /// Function that computes the gradient of the singular velocity
    /// near the corner x_centre_node_pt: \f$grad[i][j] = du_i/dx_j\f$
    return [contact_angle, &x_centre_node_pt](
             const Vector<double>& x) -> Vector<Vector<double>>
    {
      // Initialise the gradient matrix to return
      Vector<Vector<double>> grad_u(2);
      for (unsigned d = 0; d < 2; d++)
      {
        grad_u[d].resize(2);
        for (unsigned k = 0; k < 2; k++)
        {
          grad_u[d][k] = 0.0;
        }
      }

      // Find radial distance of polar coordinates centered at the point
      // x_centre
      Vector<double> x_centre = x_centre_node_pt->position();
      double r = get_radial_distance(x, x_centre);
      // Find angle of polar coordinates centered at the point x_centre
      double theta = get_azimuthal_angle(x, x_centre);

      double lambda = MathematicalConstants::Pi / contact_angle;

      const double du0dr =
        lambda * (lambda - 1) * pow(r, lambda - 2) * sin((lambda - 1) * theta);
      const double du0dtheta =
        lambda * (lambda - 1) * pow(r, lambda - 1) * cos((lambda - 1) * theta);
      const double du1dr =
        lambda * (lambda - 1) * pow(r, lambda - 2) * cos((lambda - 1) * theta);
      const double du1dtheta =
        -lambda * (lambda - 1) * pow(r, lambda - 1) * sin((lambda - 1) * theta);

      // New cartesian coordinates centred at x_centre
      Vector<double> x_bar(2, 0.0);
      x_bar[0] = x[0] - x_centre[0];
      x_bar[1] = x[1] - x_centre[1];

      const double drdx = x_bar[0] / r;
      const double dthetadx = -x_bar[1] / r / r;
      const double drdy = x_bar[1] / r;
      const double dthetady = x_bar[0] / r / r;

      // Value of \frac{\partial u}{\partial x}
      grad_u[0][0] = du0dr * drdx + du0dtheta * dthetadx;
      // Value of \frac{\partial u}{\partial y}
      grad_u[0][1] = du0dr * drdy + du0dtheta * dthetady;
      // Value of \frac{\partial v}{\partial x}
      grad_u[1][0] = du1dr * drdx + du1dtheta * dthetadx;
      // Value of \frac{\partial u}{\partial y}
      grad_u[1][1] = du1dr * drdy + du1dtheta * dthetady;

      return grad_u;
    }; // End of lambda function
  }

  /// Function that returns an eigensolution traction function with the
  /// appropriate parameters.
  std::function<void(const double&,
                     const Vector<double>&,
                     const Vector<double>&,
                     Vector<double>&)>
  eigensolution_traction_function_factory(
    const double& contact_angle,
    const std::function<Vector<Vector<double>>(const Vector<double>&)>&
      grad_velocity_singular_fct)
  {
    // Return a lambda function that captures the contact angle and the gradient
    // of the singular velocity function grad_velocity_singular_fct and takes
    // 't', 'x', 'n' and 'result' as arguments The lambda function computes the
    // traction at a point x with normal n using the gradient of the singular
    // velocity function grad_velocity_singular_fct The traction is computed as
    // - (grad u + grad u ^ T) dot n where u is the singular velocity function
    return [contact_angle,
            &grad_velocity_singular_fct](const double& t,
                                         const Vector<double>& x,
                                         const Vector<double>& n,
                                         Vector<double>& result) -> void
    {
      Vector<Vector<double>> grad_u = grad_velocity_singular_fct(x);
      const unsigned dim = x.size();
      for (unsigned i = 0; i < dim; i++)
      {
        for (unsigned j = 0; j < dim; j++)
        {
          result[i] -= (grad_u[i][j] + grad_u[j][i]) * n[j];
        }
      }
    };
  }


  /// Function that returns an eigensolution slip function with the
  /// appropriate parameters.
  std::function<void(const double&,
                     const Vector<double>&,
                     const Vector<double>&,
                     Vector<double>&)>
  eigensolution_slip_function_factory(
    const double& slip_length,
    const std::function<Vector<double>(const Vector<double>&)>&
      velocity_singular_fct)
  {
    // Return a lambda that captures 's' and takes 't' and 'x' as arguments
    // The lambda function computes the eigensolution's  contribuiton to the
    // slip condition.
    return [slip_length, velocity_singular_fct](const double& t,
                                                const Vector<double>& x,
                                                const Vector<double>& n,
                                                Vector<double>& result) -> void
    {
      Vector<double> v = velocity_singular_fct(x);
      result[0] = 0.0;
      result[1] = -v[1] / slip_length;
    };
  }

} // namespace oomph
