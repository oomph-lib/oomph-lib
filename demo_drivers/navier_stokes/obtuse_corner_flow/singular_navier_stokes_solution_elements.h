// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Header file for Navier Stokes elements with singularity
#ifndef SINGULAR_NAVIER_STOKES_SOLUTION_ELEMENTS_HEADER
#define SINGULAR_NAVIER_STOKES_SOLUTION_ELEMENTS_HEADER

#include "navier_stokes.h"

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

namespace oomph
{
  //==================CLASS FOR THE ADDITIONAL UNKNOWN==================
  /// We consider a singularity in the solution at the point O.
  ///
  /// This class defines the singular functions.
  ///
  /// velocity_singularity:
  ///
  ///    u_bar = C*velocity_singular_function
  ///
  /// pressure_singularity:
  ///
  ///    p_bar = C*pressure_singular_function
  ///
  /// and their gradients
  ///
  /// The class also defines the function that computes the
  /// residual associated to C which is:
  ///
  /// R_C = \frac{\partial p_FE}{\partial x_i} (O)
  ///
  /// and thus regularises the FE solution by setting the pressure gradient
  /// in the coordinate direction x_i to zero. If the amplitude of the
  /// singular solution is known, pin C.
  //=====================================================================
  template<class WRAPPED_NAVIER_STOKES_ELEMENT>
  class SingularNavierStokesSolutionElement : public virtual GeneralisedElement
  {
  public:
    /// Function pointer to the velocity singular function:
    typedef std::function<Vector<double>(const Vector<double>&)>
      NavierStokesVelocitySingularFct;

    /// Function pointer to the gradient of the velocity singular function:
    typedef std::function<Vector<Vector<double>>(const Vector<double>&)>
      NavierStokesGradVelocitySingularFct;

    /// Function pointer to the pressure singular function:
    typedef double (*NavierStokesPressureSingularFctPt)(
      const Vector<double>& x);

    /// Function pointer to the gradient of the pressure singular function:
    typedef Vector<double> (*NavierStokesGradPressureSingularFctPt)(
      const Vector<double>& x);

  private:
    /// Pointer to velocity singular function
    NavierStokesVelocitySingularFct Velocity_singular_fct;

    /// Pointer to gradient of velocity singular function;
    /// grad[i][j] = du_i/dx_j
    NavierStokesGradVelocitySingularFct Grad_velocity_singular_fct;

    /// Pointer to pressure singular function
    NavierStokesPressureSingularFctPt Pressure_singular_fct_pt;

    /// Pointer to gradient of pressure singular function
    NavierStokesGradPressureSingularFctPt Grad_pressure_singular_fct_pt;

    /// Local coordinates of singulariity in wrapped Navier-Stokes element
    Vector<double> S_in_wrapped_navier_stokes_element;

    /// Direction of the derivative used for the residual of the element
    unsigned* Direction_pt;

    // Does singular fct satisfy Stokes eqn?
    bool Singular_function_satisfies_stokes_equation;

  public:
    /// Constructor
    SingularNavierStokesSolutionElement()
    {
      // Initialise Function pointer to pressure singular function to NULL
      Pressure_singular_fct_pt = 0;

      // Initialise Function pointer to gradient of pressure singular
      // function to NULL
      Grad_pressure_singular_fct_pt = 0;

      // Initialise the pointer to the direction of the derivative used
      // for the residual of this element
      Direction_pt = 0;

      // Safe assumption: Singular fct does not satisfy Stokes eqn
      Singular_function_satisfies_stokes_equation = true;

      // Create a single item of internal Data, storing one unknown which
      // represents the unknown C.
      add_internal_data(new Data(1));
    }

    void fill_in_contribution_to_dresiduals_dparameter(
      double* const& parameter_pt, Vector<double>& dres_dparam)
    {
    }

    /// Assert that singular function satisfies the Stokes equations by setting
    /// this to true or false.
    bool& singular_function_satisfies_stokes_equation()
    {
      return Singular_function_satisfies_stokes_equation;
    }

    /// Find the value of the unknown amplitude C
    double c() const
    {
      return internal_data_pt(0)->value(0);
    }


    /// Find the value of the unknown amplitude C
    void set_c(const double& value)
    {
      internal_data_pt(0)->set_value(0, value);
    }

    /// Pin the value of the unknown amplitude C
    void pin_c()
    {
      return internal_data_pt(0)->pin(0);
    }

    /// Unpin the value of the unknown amplitude C
    void unpin_c()
    {
      return internal_data_pt(0)->unpin(0);
    }

    /// Access function to pointer to velocity singular function
    NavierStokesVelocitySingularFct& velocity_singular_fct()
    {
      return Velocity_singular_fct;
    }

    /// Access function to pointer to gradient of velocity singular function
    NavierStokesGradVelocitySingularFct& grad_velocity_singular_fct()
    {
      return Grad_velocity_singular_fct;
    }

    /// Access function to pointer to pressure singular function
    NavierStokesPressureSingularFctPt& pressure_singular_fct_pt()
    {
      return Pressure_singular_fct_pt;
    }

    /// Access function to pointer to gradient of pressure singular function
    NavierStokesGradPressureSingularFctPt& grad_pressure_singular_fct_pt()
    {
      return Grad_pressure_singular_fct_pt;
    }

    /// Evaluate velocity singular function at Eulerian position x
    Vector<double> velocity_singular_function(const Vector<double>& x) const
    {
#ifdef PARANOID
      if (Velocity_singular_fct == 0)
      {
        std::stringstream error_stream;
        error_stream
          << "Pointer to velocity singular function hasn't been defined!"
          << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Evaluate velocity singular function
      return Velocity_singular_fct(x);
    }

    /// Evaluate gradient of velocity singular function at Eulerian
    /// position x. grad[i][j] = du_i/dx_j
    Vector<Vector<double>> grad_velocity_singular_function(
      const Vector<double>& x) const
    {
#ifdef PARANOID
      if (Grad_velocity_singular_fct == 0)
      {
        std::stringstream error_stream;
        error_stream << "Pointer to gradient of velocity singular function "
                     << "hasn't been defined!" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Evaluate gradient of velocity singular function
      return Grad_velocity_singular_fct(x);
    }

    /// Evaluate pressure singular function at Eulerian position x
    double pressure_singular_function(const Vector<double>& x) const
    {
#ifdef PARANOID
      if (Pressure_singular_fct_pt == 0)
      {
        std::stringstream error_stream;
        error_stream
          << "Pointer to pressure singular function hasn't been defined!"
          << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Evaluate pressure singular function
      return (*Pressure_singular_fct_pt)(x);
    }

    /// Evaluate gradient of pressure singular function at Eulerian position x
    Vector<double> grad_pressure_singular_function(
      const Vector<double>& x) const
    {
#ifdef PARANOID
      if (Grad_pressure_singular_fct_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Pointer to gradient of pressure singular function "
                     << "hasn't been defined!" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Evaluate gradient of pressure singular function
      return (*Grad_pressure_singular_fct_pt)(x);
    }

    /// Evaluate velocity singular function (including its amplitude):
    /// u_bar = C * velocity_singular
    Vector<double> u_bar(const Vector<double>& x)
    {
      // Find the value of C
      double c = internal_data_pt(0)->value(0);

      // Initialise the velocity vector
      Vector<double> u = velocity_singular_function(x);

      // Find the dimension of the problem
      unsigned cached_dim = 2;

      // Multiply the components of the velocity vector by the unknown C
      for (unsigned d = 0; d < cached_dim; d++)
      {
        u[d] *= c;
      }

      // Value of u_bar at the position x
      return u;
    }

    /// Evaluate gradient of velocity singular function
    /// (including its amplitude):
    /// grad_u_bar = C * grad_velocity_singular;
    /// grad[i][j] = du_i/dx_j
    Vector<Vector<double>> grad_u_bar(const Vector<double>& x)
    {
      // Find the value of C
      double c = internal_data_pt(0)->value(0);

      // Initialise the gradient of velocity vector
      Vector<Vector<double>> grad_u = grad_velocity_singular_function(x);

      // Find the dimension of the problem
      unsigned cached_dim = 2;

      // Multiply the components of the gradient of velocity by the unknown C
      for (unsigned d = 0; d < cached_dim; d++)
      {
        for (unsigned i = 0; i < cached_dim; i++)
        {
          grad_u[d][i] *= c;
        }
      }

      // Value of grad_u_bar at the position x
      return grad_u;
    }

    /// Evaluate pressure singular function (including its amplitude):
    /// p_bar = C * pressure_singular
    double p_bar(const Vector<double>& x)
    {
      // Find the value of C
      double c = internal_data_pt(0)->value(0);

      // Value of p_bar at the position x
      return c * pressure_singular_function(x);
    }

    /// Evaluate gradient of pressure singular function
    /// (including its amplitude):
    /// grad_p_bar = C * grad_pressure_singular
    Vector<double> grad_p_bar(const Vector<double>& x)
    {
      // Find the value of C
      double c = internal_data_pt(0)->value(0);

      // Initialise the gradient of pressure
      Vector<double> grad_p = grad_pressure_singular_function(x);

      // Find the dimension of the problem
      unsigned cached_dim = 2;

      // Multiply the components of the gradient of pressure by the unknown C
      for (unsigned d = 0; d < cached_dim; d++)
      {
        grad_p[d] *= c;
      }

      // Value of grad_p_bar at the position x
      return grad_p;
    }

    /// Overload the output function
    void output(std::ostream& outfile)
    {
      outfile << std::fixed
              << std::setprecision(std::numeric_limits<double>::max_digits10);
      outfile << internal_data_pt(0)->value(0) << endl;
    }


    /// Compute residual
    void fill_in_contribution_to_residuals(Vector<double>& residuals) {}

    // Compute local residual and jacobian
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      fill_in_contribution_to_jacobian(residuals, jacobian);
    }

  }; // End of SingularNavierStokesSolutionElement class
} // namespace oomph
#endif
