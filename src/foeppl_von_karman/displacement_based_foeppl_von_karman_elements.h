// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision$
// LIC//
// LIC// $LastChangedDate$
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
// Header file for FoepplvonKarman elements
#ifndef OOMPH_FOEPPLVONKARMAN_DISPLACEMENT_ELEMENTS_HEADER
#define OOMPH_FOEPPLVONKARMAN_DISPLACEMENT_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include <sstream>

// OOMPH-LIB headers
#include "../generic/projection.h"
#include "../generic/nodes.h"
#include "../generic/oomph_utilities.h"

namespace oomph
{
  //=============================================================
  /// A class for all isoparametric elements that solve the
  /// Foeppl von Karman equations.
  /// \f[
  /// \nabla^4 w -  \eta  \frac{\partial}{\partial x_{\beta}}
  ///                     \left( \sigma^{\alpha \beta}
  ///                            \frac{\partial w}{\partial x_{\alpha}}
  ///                     \right) = p(x,y)
  /// \f]
  /// and
  /// \f[
  /// \frac{\sigma^{\alpha \beta}}{\partial x_{\beta}} = \tau_\alpha
  /// \f]
  /// This contains the generic maths. Shape functions, geometric
  /// mapping etc. must get implemented in derived class.
  //=============================================================
  class DisplacementBasedFoepplvonKarmanEquations : public virtual FiniteElement
  {
  public:
    /// \short Function pointer to pressure function fct(x,f(x)) --
    /// x is a Vector!
    typedef void (*FoepplvonKarmanPressureFctPt)(const Vector<double> &x,
                                                 double &f);

    /// \short Function pointer to traction function fct(x,f(x)) --
    /// x is a Vector!
    typedef void (*FoepplvonKarmanTractionFctPt)(const Vector<double> &x,
                                                 Vector<double> &f);

    /// \short Constructor (must initialise the Pressure_fct_pt and the
    /// Traction_fct_pt. Also set physical parameters to their default
    /// values.
    DisplacementBasedFoepplvonKarmanEquations() :
      Pressure_fct_pt(0), Traction_fct_pt(0)
    {
      // Set default
      Nu_pt = &Default_Nu_Value;

      // Set all the physical constants to the default value (zero)
      Eta_pt = &Default_Physical_Constant_Value;

      // Linear bending model?
      Linear_bending_model = false;
    }

    /// Broken copy constructor
    DisplacementBasedFoepplvonKarmanEquations(
      const DisplacementBasedFoepplvonKarmanEquations &dummy)
    {
      BrokenCopy::broken_copy("DisplacementBasedFoepplvonKarmanEquations");
    }

    /// Broken assignment operator
    void operator=(const DisplacementBasedFoepplvonKarmanEquations &)
    {
      BrokenCopy::broken_assign("DisplacementBasedFoepplvonKarmanEquations");
    }

    /// Poisson's ratio
    const double &nu() const
    {
      return *Nu_pt;
    }

    /// Pointer to Poisson's ratio
    double *&nu_pt()
    {
      return Nu_pt;
    }

    /// Eta
    const double &eta() const
    {
      return *Eta_pt;
    }

    /// Pointer to eta
    double *&eta_pt()
    {
      return Eta_pt;
    }

    /// \short Return the index at which the i-th unknown value
    /// is stored. The default value, i, is appropriate for single-physics
    /// problems. By default, these are:
    /// 0: w
    /// 1: laplacian w
    /// 2: u_x
    /// 3: u_y
    /// In derived multi-physics elements, this function should be overloaded
    /// to reflect the chosen storage scheme. Note that these equations require
    /// that the unknown is always stored at the same index at each node.
    virtual inline unsigned nodal_index_fvk(const unsigned &i = 0) const
    {
      return i;
    }

    /// Output with default number of plot points
    void output(std::ostream &outfile)
    {
      const unsigned n_plot = 5;
      output(outfile, n_plot);
    }

    /// \short Output FE representation of soln: x,y,w at
    /// n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot);

    /// C_style output with default number of plot points
    void output(FILE *file_pt)
    {
      const unsigned n_plot = 5;
      output(file_pt, n_plot);
    }

    /// \short C-style output FE representation of soln: x,y,w at
    /// n_plot^DIM plot points
    void output(FILE *file_pt, const unsigned &n_plot);

    /// Output exact soln: x,y,w_exact at n_plot^DIM plot points
    void output_fct(std::ostream &outfile,
                    const unsigned &n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

    /// \short Output exact soln: x,y,w_exact at
    /// n_plot^DIM plot points (dummy time-dependent version to
    /// keep intel compiler happy)
    virtual void output_fct(
      std::ostream &outfile,
      const unsigned &n_plot,
      const double &time,
      FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      throw OomphLibError(
        "There is no time-dependent output_fct() for Foeppl von Karman"
        "elements ",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Get error against and norm of exact solution
    void compute_error(std::ostream &outfile,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       double &error,
                       double &norm);

    /// Dummy, time dependent error checker
    void compute_error(std::ostream &outfile,
                       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                       const double &time,
                       double &error,
                       double &norm)
    {
      throw OomphLibError(
        "There is no time-dependent compute_error() for Foeppl von Karman"
        "elements",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Access function: Pointer to pressure function
    FoepplvonKarmanPressureFctPt &pressure_fct_pt()
    {
      return Pressure_fct_pt;
    }

    /// Access function: Pointer to pressure function. Const version
    FoepplvonKarmanPressureFctPt pressure_fct_pt() const
    {
      return Pressure_fct_pt;
    }

    /// Access function: Pointer to in-plane traction function
    FoepplvonKarmanTractionFctPt &traction_fct_pt()
    {
      return Traction_fct_pt;
    }

    /// Access function: Pointer to in-plane traction function. Const version
    FoepplvonKarmanTractionFctPt traction_fct_pt() const
    {
      return Traction_fct_pt;
    }

    /// \short Get pressure term at (Eulerian) position x. This function
    /// is virtual to allow overloading in multi-physics problems where
    /// the strength of the pressure function might be determined by
    /// another system of equations.
    inline virtual void get_pressure_fvk(const unsigned &ipt,
                                         const Vector<double> &x,
                                         double &pressure) const
    {
      // If no pressure function has been set, return zero
      if (Pressure_fct_pt == 0)
      {
        pressure = 0.0;
      }
      else
      {
        // Get pressure strength
        (*Pressure_fct_pt)(x, pressure);
      }
    }

    /// \short Get in-plane traction term at (Eulerian) position x.
    inline virtual void get_traction_fvk(Vector<double> &x,
                                         Vector<double> &traction) const
    {
      // If no pressure function has been set, return zero
      if (Traction_fct_pt == 0)
      {
        traction[0] = 0.0;
        traction[1] = 0.0;
      }
      else
      {
        // Get traction
        (*Traction_fct_pt)(x, traction);
      }
    }

    /// Get gradient of deflection: gradient[i] = dw/dx_i
    void get_gradient_of_deflection(const Vector<double> &s,
                                    Vector<double> &gradient) const
    {
      // Find out how many nodes there are in the element
      const unsigned n_node = nnode();

      // Get the index at which the unknown is stored
      // Indexes for unknows are
      // 0: W (vertical displacement)
      // 1: L (laplacian of W)
      // 2: Ux (In plane displacement x)
      // 3: Uy (In plane displacement y)
      unsigned w_nodal_index = nodal_index_fvk(0);

      // Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidx(n_node, 2);

      // Call the derivatives of the shape and test functions
      dshape_eulerian(s, psi, dpsidx);

      // Initialise to zero
      for (unsigned j = 0; j < 2; j++)
      {
        gradient[j] = 0.0;
      }

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over derivative directions
        for (unsigned j = 0; j < 2; j++)
        {
          gradient[j] += this->nodal_value(l, w_nodal_index) * dpsidx(l, j);
        }
      }
    }

    /// Get gradient of field: gradient[i] = d[.]/dx_i,
    // Indices for fields are
    // 0: W (vertical displacement)
    // 1: L (laplacian of W)
    // 2: Ux (In plane displacement x)
    // 3: Uy (In plane displacement y)
    void get_gradient_of_field(const Vector<double> &s,
                               Vector<double> &gradient,
                               const unsigned &index) const
    {
      // Find out how many nodes there are in the element
      const unsigned n_node = nnode();

      // Get the index at which the unknown is stored
      // Indexes for unknows are
      // 0: W (vertical displacement)
      // 1: L (laplacian of W)
      // 2: Ux (In plane displacement x)
      // 3: Uy (In plane displacement y)
      const unsigned w_nodal_index = nodal_index_fvk(index);

      // Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidx(n_node, 2);

      // Call the derivatives of the shape and test functions
      dshape_eulerian(s, psi, dpsidx);

      // Initialise to zero
      for (unsigned j = 0; j < 2; j++)
      {
        gradient[j] = 0.0;
      }

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over derivative directions
        for (unsigned j = 0; j < 2; j++)
        {
          gradient[j] += this->nodal_value(l, w_nodal_index) * dpsidx(l, j);
        }
      }
    }

    /// Get the in-plane stress (sigma) as a fct of the pre=computed
    /// displcement derivatives
    void get_sigma(DenseMatrix<double> &sigma,
                   const Vector<double> &interpolated_dwdx,
                   const Vector<double> &interpolated_duxdx,
                   const Vector<double> &interpolated_duydx)
    {
      // The strain tensor values
      double e_xx = interpolated_duxdx[0];
      double e_yy = interpolated_duydx[1];
      double e_xy = 0.5 * (interpolated_duxdx[1] + interpolated_duydx[0]);
      e_xx += 0.5 * interpolated_dwdx[0] * interpolated_dwdx[0];
      e_yy += 0.5 * interpolated_dwdx[1] * interpolated_dwdx[1];
      e_xy += 0.5 * interpolated_dwdx[0] * interpolated_dwdx[1];

      // Get Poisson's ratio
      const double _nu = nu();

      double inv_denom = 1.0 / (1.0 - _nu * _nu);

      // The stress tensor values
      // sigma_xx
      sigma(0, 0) = (e_xx + _nu * e_yy) * inv_denom;

      // sigma_yy
      sigma(1, 1) = (e_yy + _nu * e_xx) * inv_denom;

      // sigma_xy = sigma_yx
      sigma(0, 1) = sigma(1, 0) = e_xy / (1.0 + _nu);
    }

    // Get the stress for output
    void get_stress_and_strain_for_output(const Vector<double> &s,
                                          DenseMatrix<double> &sigma,
                                          DenseMatrix<double> &strain)
    {
      // Find out how many nodes there are
      const unsigned n_node = nnode();

      // Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidx(n_node, 2);

      // Local shape function
      dshape_eulerian(s, psi, dpsidx);

      // Get the derivatives of the shape and test functions
      // double J = dshape_and_dtest_eulerian_fvk(s, psi, dpsidx, test,
      // dtestdx);

      // Indices at which the unknowns are stored
      const unsigned w_nodal_index = nodal_index_fvk(0);
      const unsigned u_x_nodal_index = nodal_index_fvk(2);
      const unsigned u_y_nodal_index = nodal_index_fvk(3);

      // Allocate and initialise to zero storage for the interpolated values

      // The variables
      Vector<double> interpolated_dwdx(2, 0.0);
      Vector<double> interpolated_duxdx(2, 0.0);
      Vector<double> interpolated_duydx(2, 0.0);

      //--------------------------------------------
      // Calculate function values and derivatives:
      //--------------------------------------------
      Vector<double> nodal_value(4, 0.0);
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the nodal values
        nodal_value[0] = this->nodal_value(l, w_nodal_index);
        nodal_value[2] = this->nodal_value(l, u_x_nodal_index);
        nodal_value[3] = this->nodal_value(l, u_y_nodal_index);

        // Add contributions from current node/shape function

        // Loop over directions
        for (unsigned j = 0; j < 2; j++)
        {
          interpolated_dwdx[j] += nodal_value[0] * dpsidx(l, j);
          interpolated_duxdx[j] += nodal_value[2] * dpsidx(l, j);
          interpolated_duydx[j] += nodal_value[3] * dpsidx(l, j);

        } // Loop over directions for (j<2)

      } // Loop over nodes for (l<n_node)

      // Get in-plane stress
      get_sigma(
        sigma, interpolated_dwdx, interpolated_duxdx, interpolated_duydx);

      // The strain tensor values
      // E_xx
      strain(0, 0) = interpolated_duxdx[0];
      strain(0, 0) += 0.5 * interpolated_dwdx[0] * interpolated_dwdx[0];

      // E_yy
      strain(1, 1) = interpolated_duydx[1];
      strain(1, 1) += 0.5 * interpolated_dwdx[1] * interpolated_dwdx[1];

      // E_xy
      strain(0, 1) = 0.5 * (interpolated_duxdx[1] + interpolated_duydx[0]);
      strain(0, 1) += 0.5 * interpolated_dwdx[0] * interpolated_dwdx[1];

      // E_yx
      strain(1, 0) = strain(0, 1);
    }

    /// hierher dummy
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double> &residuals,
      DenseMatrix<double> &jacobian,
      DenseMatrix<double> &mass_matrix)
    {
      // Get Jacobian from base class (FD-ed)
      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);

      // Dummy diagonal (won't result in global unit matrix but
      // doesn't matter for zero eigenvalue/eigenvector
      unsigned ndof = mass_matrix.nrow();
      for (unsigned i = 0; i < ndof; i++)
      {
        mass_matrix(i, i) += 1.0;
      }
    }

    /// Fill in the residuals with this element's contribution
    void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      // Find out how many nodes there are
      const unsigned n_node = nnode();

      // Set up memory for the shape and test functions
      Shape psi(n_node), test(n_node);
      DShape dpsidx(n_node, 2), dtestdx(n_node, 2);

      // Indices at which the unknowns are stored
      const unsigned w_nodal_index = nodal_index_fvk(0);
      const unsigned laplacian_w_nodal_index = nodal_index_fvk(1);
      const unsigned u_x_nodal_index = nodal_index_fvk(2);
      const unsigned u_y_nodal_index = nodal_index_fvk(3);

      // Set the value of n_intpt
      const unsigned n_intpt = integral_pt()->nweight();

      // Integers to store the local equation numbers
      int local_eqn = 0;

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Call the derivatives of the shape and test functions
        double J = dshape_and_dtest_eulerian_at_knot_fvk(
          ipt, psi, dpsidx, test, dtestdx);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Allocate and initialise to zero storage for the interpolated values
        Vector<double> interpolated_x(2, 0.0);

        // The variables
        double interpolated_laplacian_w = 0;
        Vector<double> interpolated_dwdx(2, 0.0);
        Vector<double> interpolated_dlaplacian_wdx(2, 0.0);
        Vector<double> interpolated_duxdx(2, 0.0);
        Vector<double> interpolated_duydx(2, 0.0);

        //--------------------------------------------
        // Calculate function values and derivatives:
        //--------------------------------------------
        Vector<double> nodal_value(4, 0.0);
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Get the nodal values
          nodal_value[0] = this->nodal_value(l, w_nodal_index);
          nodal_value[1] = this->nodal_value(l, laplacian_w_nodal_index);
          nodal_value[2] = this->nodal_value(l, u_x_nodal_index);
          nodal_value[3] = this->nodal_value(l, u_y_nodal_index);

          // Add contributions from current node/shape function
          interpolated_laplacian_w += nodal_value[1] * psi(l);

          // Loop over directions
          for (unsigned j = 0; j < 2; j++)
          {
            interpolated_x[j] += nodal_position(l, j) * psi(l);
            interpolated_dwdx[j] += nodal_value[0] * dpsidx(l, j);
            interpolated_dlaplacian_wdx[j] += nodal_value[1] * dpsidx(l, j);
            interpolated_duxdx[j] += nodal_value[2] * dpsidx(l, j);
            interpolated_duydx[j] += nodal_value[3] * dpsidx(l, j);

          } // Loop over directions for (j<2)

        } // Loop over nodes for (l<n_node)

        // Get in-plane stress
        DenseMatrix<double> sigma(2, 2, 0.0);

        // Stress not used if we have the linear bending model
        if (!Linear_bending_model)
        {
          // Get the value of in plane stress at the integration
          // point
          get_sigma(
            sigma, interpolated_dwdx, interpolated_duxdx, interpolated_duydx);
        }

        // Get pressure function
        //-------------------
        double pressure = 0.0;
        get_pressure_fvk(ipt, interpolated_x, pressure);

        // Assemble residuals and Jacobian
        //--------------------------------

        // Loop over the test functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Get the local equation (First equation)
          local_eqn = nodal_local_eqn(l, w_nodal_index);

          // IF it's not a boundary condition
          if (local_eqn >= 0)
          {
            residuals[local_eqn] += pressure * test(l) * W;

            // Reduced order biharmonic operator
            for (unsigned k = 0; k < 2; k++)
            {
              residuals[local_eqn] +=
                interpolated_dlaplacian_wdx[k] * dtestdx(l, k) * W;
            }

            // Sigma_alpha_beta part
            if (!Linear_bending_model)
            {
              // Alpha loop
              for (unsigned alpha = 0; alpha < 2; alpha++)
              {
                // Beta loop
                for (unsigned beta = 0; beta < 2; beta++)
                {
                  residuals[local_eqn] -= eta() * sigma(alpha, beta) *
                                          interpolated_dwdx[alpha] *
                                          dtestdx(l, beta) * W;
                }
              }
            } // if(!Linear_bending_model)
          }

          // Get the local equation (Second equation)
          local_eqn = nodal_local_eqn(l, laplacian_w_nodal_index);

          // IF it's not a boundary condition
          if (local_eqn >= 0)
          {
            // The coupled Poisson equations for the biharmonic operator
            residuals[local_eqn] += interpolated_laplacian_w * test(l) * W;

            for (unsigned k = 0; k < 2; k++)
            {
              residuals[local_eqn] += interpolated_dwdx[k] * dtestdx(l, k) * W;
            }
          }

          // Get in plane traction
          Vector<double> traction(2, 0.0);
          get_traction_fvk(interpolated_x, traction);

          // Get the local equation (Third equation)
          local_eqn = nodal_local_eqn(l, u_x_nodal_index);

          // IF it's not a boundary condition
          if (local_eqn >= 0)
          {
            // tau_x
            residuals[local_eqn] += traction[0] * test(l) * W;

            // r_{\alpha = x}
            for (unsigned beta = 0; beta < 2; beta++)
            {
              residuals[local_eqn] += sigma(0, beta) * dtestdx(l, beta) * W;
            }
          }

          // Get the local equation (Fourth equation)
          local_eqn = nodal_local_eqn(l, u_y_nodal_index);

          // IF it's not a boundary condition
          if (local_eqn >= 0)
          {
            // tau_y
            residuals[local_eqn] += traction[1] * test(l) * W;

            // r_{\alpha = y}
            for (unsigned beta = 0; beta < 2; beta++)
            {
              residuals[local_eqn] += sigma(1, beta) * dtestdx(l, beta) * W;
            }
          }

        } // End loop over nodes or test functions (l<n_node)

      } // End of loop over integration points
    }

    /// \short Return FE representation of function value w_fvk(s)
    /// at local coordinate s (by default - if index > 0, returns
    /// FE representation of valued stored at index^th nodal index
    inline double interpolated_w_fvk(const Vector<double> &s,
                                     unsigned index = 0) const
    {
      // Find number of nodes
      const unsigned n_node = nnode();

      // Get the index at which the poisson unknown is stored
      const unsigned w_nodal_index = nodal_index_fvk(index);

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise value of u
      double interpolated_w = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_w += this->nodal_value(l, w_nodal_index) * psi[l];
      }

      return (interpolated_w);
    }

    /// \short Self-test: Return 0 for OK
    unsigned self_test();

    /// \short Sets a flag to signify that we are solving the linear, pure
    /// bending equations, and pin all the nodal values that will not be used in
    /// this case
    void use_linear_bending_model()
    {
      // Set the boolean flag
      Linear_bending_model = true;

      // Get the index of the first FvK nodal value
      unsigned first_fvk_nodal_index = nodal_index_fvk();

      // Get the total number of FvK nodal values (assuming they are stored
      // contiguously) at node 0 (it's the same at all nodes anyway)
      unsigned total_fvk_nodal_indices = 4;

      // Get the number of nodes in this element
      unsigned n_node = nnode();

      // Loop over the appropriate nodal indices
      for (unsigned index = first_fvk_nodal_index + 2; // because [2] is u_x
                                                       // and [3] is u_y
           index < first_fvk_nodal_index + total_fvk_nodal_indices;
           index++)
      {
        // Loop over the nodes in the element
        for (unsigned inod = 0; inod < n_node; inod++)
        {
          // Pin the nodal value at the current index
          node_pt(inod)->pin(index);
        }
      }
    }

  protected:
    /// \short Shape/test functions and derivs w.r.t. to global coords at
    /// local coord. s; return  Jacobian of mapping
    virtual double dshape_and_dtest_eulerian_fvk(const Vector<double> &s,
                                                 Shape &psi,
                                                 DShape &dpsidx,
                                                 Shape &test,
                                                 DShape &dtestdx) const = 0;

    /// \short Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return  Jacobian of mapping
    virtual double dshape_and_dtest_eulerian_at_knot_fvk(
      const unsigned &ipt,
      Shape &psi,
      DShape &dpsidx,
      Shape &test,
      DShape &dtestdx) const = 0;

    /// Pointer to global Poisson's ratio
    double *Nu_pt;

    /// Pointer to global eta
    double *Eta_pt;

    /// Pointer to pressure function:
    FoepplvonKarmanPressureFctPt Pressure_fct_pt;

    /// Pointer to traction function:
    FoepplvonKarmanTractionFctPt Traction_fct_pt;

  private:
    /// Default value for Poisson's ratio
    static double Default_Nu_Value;

    /// Default value for physical constants
    static double Default_Physical_Constant_Value;

    /// \short Flag which stores whether we are using a linear, pure
    /// bending model instead of the full non-linear Foeppl-von Karman
    bool Linear_bending_model;
  };

  //==========================================================
  /// Foeppl von Karman upgraded to become projectable
  //==========================================================
  template<class FVK_ELEMENT>
  class ProjectableDisplacementBasedFoepplvonKarmanElement :
    public virtual ProjectableElement<FVK_ELEMENT>
  {
  public:
    /// \short Specify the values associated with field fld.  The
    /// information is returned in a vector of pairs which comprise the
    /// Data object and the value within it, that correspond to field
    /// fld.
    Vector<std::pair<Data *, unsigned>> data_values_of_field(
      const unsigned &fld)
    {
#ifdef PARANOID
      if (fld > 3)
      {
        std::stringstream error_stream;
        error_stream
          << "Foeppl von Karman elements only store four fields so fld must be"
          << "0 to 3 rather than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Create the vector
      unsigned nnod = this->nnode();
      Vector<std::pair<Data *, unsigned>> data_values(nnod);

      // Loop over all nodes
      for (unsigned j = 0; j < nnod; j++)
      {
        // Add the data value associated field: The node itself
        data_values[j] = std::make_pair(this->node_pt(j), fld);
      }

      // Return the vector
      return data_values;
    }

    /// \short Number of fields to be projected: Just two
    unsigned nfields_for_projection()
    {
      return 4;
    }

    /// \short Number of history values to be stored for fld-th field.
    /// (Note: count includes current value!)
    unsigned nhistory_values_for_projection(const unsigned &fld)
    {
#ifdef PARANOID
      if (fld > 3)
      {
        std::stringstream error_stream;
        error_stream
          << "Foeppl von Karman elements only store four fields so fld must be"
          << "0 to 3 rather than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return this->node_pt(0)->ntstorage();
    }

    ///\short Number of positional history values
    /// (Note: count includes current value!)
    unsigned nhistory_values_for_coordinate_projection()
    {
      return this->node_pt(0)->position_time_stepper_pt()->ntstorage();
    }

    /// \short Return Jacobian of mapping and shape functions of field fld
    /// at local coordinate s
    double jacobian_and_shape_of_field(const unsigned &fld,
                                       const Vector<double> &s,
                                       Shape &psi)
    {
#ifdef PARANOID
      if (fld > 3)
      {
        std::stringstream error_stream;
        error_stream
          << "Foeppl von Karman elements only store four fields so fld must be"
          << "0 to 3 rather than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      unsigned n_dim = this->dim();
      unsigned n_node = this->nnode();
      Shape test(n_node);
      DShape dpsidx(n_node, n_dim), dtestdx(n_node, n_dim);
      double J =
        this->dshape_and_dtest_eulerian_fvk(s, psi, dpsidx, test, dtestdx);
      return J;
    }

    /// \short Return interpolated field fld at local coordinate s, at
    /// time level t (t=0: present; t>0: history values)
    double get_field(const unsigned &t,
                     const unsigned &fld,
                     const Vector<double> &s)
    {
#ifdef PARANOID
      if (fld > 3)
      {
        std::stringstream error_stream;
        error_stream
          << "Foeppl von Karman elements only store four fields so fld must be"
          << "0 to 3 rather than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // Find the index at which the variable is stored
      unsigned w_nodal_index = this->nodal_index_fvk(fld);

      // Local shape function
      unsigned n_node = this->nnode();
      Shape psi(n_node);

      // Find values of shape function
      this->shape(s, psi);

      // Initialise value of u
      double interpolated_w = 0.0;

      // Sum over the local nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_w += this->nodal_value(t, l, w_nodal_index) * psi[l];
      }
      return interpolated_w;
    }

    /// Return number of values in field fld: One per node
    unsigned nvalue_of_field(const unsigned &fld)
    {
#ifdef PARANOID
      if (fld > 3)
      {
        std::stringstream error_stream;
        error_stream
          << "Foeppl von Karman elements only store four fields so fld must be"
          << "0 to 3 rather than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return this->nnode();
    }

    /// Return local equation number of field fld of node j.
    int local_equation(const unsigned &fld, const unsigned &j)
    {
#ifdef PARANOID
      if (fld > 3)
      {
        std::stringstream error_stream;
        error_stream
          << "Foeppl von Karman elements only store four fields so fld must be"
          << "0 to 3 rather than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      const unsigned w_nodal_index = this->nodal_index_fvk(fld);
      return this->nodal_local_eqn(j, w_nodal_index);
    }
  };

  //=======================================================================
  /// Face geometry for element is the same as that for the underlying
  /// wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<
    ProjectableDisplacementBasedFoepplvonKarmanElement<ELEMENT>> :
    public virtual FaceGeometry<ELEMENT>
  {
  public:
    FaceGeometry() : FaceGeometry<ELEMENT>() {}
  };

  //=======================================================================
  /// Face geometry of the Face Geometry for element is the same as
  /// that for the underlying wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<
    FaceGeometry<ProjectableDisplacementBasedFoepplvonKarmanElement<ELEMENT>>> :
    public virtual FaceGeometry<FaceGeometry<ELEMENT>>
  {
  public:
    FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT>>() {}
  };

} // namespace oomph

#endif
