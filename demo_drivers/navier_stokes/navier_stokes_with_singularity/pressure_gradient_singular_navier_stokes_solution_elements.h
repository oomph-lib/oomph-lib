// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_PRESSURE_GRADIENT_SINGULAR_NAVIER_STOKES_SOLUTION_ELEMENTS_HEADER
#define OOMPH_PRESSURE_GRADIENT_SINGULAR_NAVIER_STOKES_SOLUTION_ELEMENTS_HEADER

#include "navier_stokes.h"

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
  class PressureGradientSingularNavierStokesSolutionElement
    : public virtual SingularNavierStokesSolutionElement<
        WRAPPED_NAVIER_STOKES_ELEMENT>
  {
  public:
    /// Compute residual
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      fill_in_generic_contribution_to_residuals(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    // Compute local residual and jacobian
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      fill_in_generic_contribution_to_residuals(residuals, jacobian, 1);
    }

  private:
    /// Compute local residual, and, if flag=1, local jacobian matrix
    void fill_in_generic_contribution_to_residuals(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      // Get the local eqn number of our one-and-only
      // unknown
      int eqn_number = this->internal_local_eqn(0, 0);
      if (eqn_number >= 0)
      {
        residuals[eqn_number] =
          this->wrapped_navier_stokes_el_pt()->dpdx_fe_only(
            this->s_in_wrapped_navier_stokes_element(), this->direction_pt());

        // Do we want the Jacobian too?
        if (flag)
        {
          // Find the number of pressure dofs in the wrapped Navier-Stokes
          // element pointed by
          // the SingularNavierStokesSolutionElement class
          unsigned n_pres = this->wrapped_navier_stokes_el_pt()->npres_nst();

          // Find the dimension of the problem
          unsigned cached_dim = this->wrapped_navier_stokes_el_pt()->dim();

          // Set up memory for the pressure shape functions and their
          // derivatives
          Shape psip(n_pres), testp(n_pres);
          DShape dpsipdx(n_pres, cached_dim), dtestpdx(n_pres, cached_dim);

          // Compute the pressure shape functions and their derivatives
          // at the local coordinate S_in_wrapped_navier_stokes_element
          // (Test fcts not really needed but nobody's got around to writing
          // a fct that only picks out the basis fcts.
          this->wrapped_navier_stokes_el_pt()->dpshape_and_dptest_eulerian_nst(
            this->s_in_wrapped_navier_stokes_element(),
            psip,
            dpsipdx,
            testp,
            dtestpdx);

          // Derivs
          for (unsigned j = 0; j < n_pres; j++)
          {
            // Unknown
            int local_unknown =
              this->wrapped_navier_stokes_el_pt()->p_local_eqn(j);

            // If not pinned
            if (local_unknown >= 0)
            {
              // Add the contribution of the node to the local jacobian
              jacobian(eqn_number, local_unknown) =
                dpsipdx(j, *this->direction_pt());
            }
          }
        }
      }
    }
  };

} // namespace oomph

#endif
