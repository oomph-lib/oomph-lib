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
// Header for a multi-physics elements that couple the Navier--Stokes
// and advection diffusion elements via a multi domain approach,
// giving Boussinesq convection

#ifndef OOMPH_BOUSSINESQ_ELEMENTS_HEADER
#define OOMPH_BOUSSINESQ_ELEMENTS_HEADER

// Oomph-lib headers, we require the generic, advection-diffusion
// and navier-stokes elements.
#include "generic.h"
#include "advection_diffusion.h"
#include "navier_stokes.h"


namespace oomph
{
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////


  //======================class definition==============================
  /// A class that solves the Boussinesq approximation of the Navier--Stokes
  /// and energy equations by coupling two pre-existing classes.
  /// The QAdvectionDiffusionElement with bi-quadratic interpolation for the
  /// scalar variable (temperature) and
  /// QCrouzeixRaviartElement which solves the Navier--Stokes equations
  /// using bi-quadratic interpolation for the velocities and a discontinuous
  /// bi-linear interpolation for the pressure. Note that we are free to
  /// choose the order in which we store the variables at the nodes. In this
  /// case we choose to store the variables in the order fluid velocities
  /// followed by temperature. We must, therefore, overload the function
  /// AdvectionDiffusionEquations<DIM>::u_index_adv_diff() to indicate that
  /// the temperature is stored at the DIM-th position not the 0-th. We do not
  /// need to overload the corresponding function in the
  /// NavierStokesEquations<DIM> class because the velocities are stored
  /// first.
  //=========================================================================
  template<unsigned DIM>
  class BuoyantQCrouzeixRaviartElement
    : public virtual QAdvectionDiffusionElement<DIM, 3>,
      public virtual QCrouzeixRaviartElement<DIM>
  {
  private:
    /// Pointer to a private data member, the Rayleigh number
    double* Ra_pt;

    /// The static default value of the Rayleigh number
    static double Default_Physical_Constant_Value;

  public:
    /// Constructor: call the underlying constructors and
    /// initialise the pointer to the Rayleigh number to point
    /// to the default value of 0.0.
    BuoyantQCrouzeixRaviartElement()
      : QAdvectionDiffusionElement<DIM, 3>(), QCrouzeixRaviartElement<DIM>()
    {
      Ra_pt = &Default_Physical_Constant_Value;
    }

    /// Unpin p_dof-th pressure dof
    void unfix_pressure(const unsigned& p_dof)
    {
      this->internal_data_pt(this->P_nst_internal_index)->unpin(p_dof);
    }

    /// The required number of values stored at the nodes is the sum of
    /// the required values of the two single-physics  elements. Note that this
    /// step is generic for any multi-physics element of this type.
    unsigned required_nvalue(const unsigned& n) const
    {
      return (QAdvectionDiffusionElement<DIM, 3>::required_nvalue(n) +
              QCrouzeixRaviartElement<DIM>::required_nvalue(n));
    }

    /// Access function for the Rayleigh number (const version)
    const double& ra() const
    {
      return *Ra_pt;
    }

    /// Access function for the pointer to the Rayleigh number
    double*& ra_pt()
    {
      return Ra_pt;
    }

    /// Final override for disable ALE
    void disable_ALE()
    {
      // Disable ALE in both sets of equations
      NavierStokesEquations<DIM>::disable_ALE();
      AdvectionDiffusionEquations<DIM>::disable_ALE();
    }

    /// Final override for enable ALE
    void enable_ALE()
    {
      // Enable ALE in both sets of equations
      NavierStokesEquations<DIM>::enable_ALE();
      AdvectionDiffusionEquations<DIM>::enable_ALE();
    }

    /// Number of scalars/fields output by this element. Broken
    /// virtual. Needs to be implemented for each new specific element type.
    /// Temporary dummy
    unsigned nscalar_paraview() const
    {
      throw OomphLibError(
        "This function hasn't been implemented for this element",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);

      // Dummy unsigned
      return 0;
    }

    /// Write values of the i-th scalar field at the plot points. Broken
    /// virtual. Needs to be implemented for each new specific element type.
    /// Temporary dummy
    void scalar_value_paraview(std::ofstream& file_out,
                               const unsigned& i,
                               const unsigned& nplot) const
    {
      throw OomphLibError(
        "This function hasn't been implemented for this element",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }


    /// Name of the i-th scalar field. Default implementation
    /// returns V1 for the first one, V2 for the second etc. Can (should!) be
    /// overloaded with more meaningful names.
    std::string scalar_name_paraview(const unsigned& i) const
    {
      return "V" + StringConversion::to_string(i);
    }


    ///  Overload the standard output function with the broken default
    void output(std::ostream& outfile)
    {
      FiniteElement::output(outfile);
    }

    /// Output function:
    ///  Output x, y, u, v, p, theta at Nplot^DIM plot points
    // Start of output function
    void output(std::ostream& outfile, const unsigned& nplot)
    {
      // vector of local coordinates
      Vector<double> s(DIM);

      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, nplot, s);

        // Output the position of the plot point
        for (unsigned i = 0; i < DIM; i++)
        {
          outfile << this->interpolated_x(s, i) << " ";
        }

        // Output the fluid velocities at the plot point
        for (unsigned i = 0; i < DIM; i++)
        {
          outfile << this->interpolated_u_nst(s, i) << " ";
        }

        // Output the fluid pressure at the plot point
        outfile << this->interpolated_p_nst(s) << " ";

        // Output the temperature (the advected variable) at the plot point
        outfile << this->interpolated_u_adv_diff(s) << std::endl;
      }
      outfile << std::endl;

      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);
    } // End of output function


    /// C-style output function: Broken default
    void output(FILE* file_pt)
    {
      FiniteElement::output(file_pt);
    }

    ///  C-style output function: Broken default
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FiniteElement::output(file_pt, n_plot);
    }

    /// Output function for an exact solution: Broken default
    void output_fct(std::ostream& outfile,
                    const unsigned& Nplot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      FiniteElement::output_fct(outfile, Nplot, exact_soln_pt);
    }


    /// Output function for a time-dependent exact solution:
    /// Broken default.
    void output_fct(std::ostream& outfile,
                    const unsigned& Nplot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      FiniteElement::output_fct(outfile, Nplot, time, exact_soln_pt);
    }

    /// Overload the index at which the temperature
    /// variable is stored. We choose to store it after the fluid velocities.
    inline unsigned u_index_adv_diff() const
    {
      return DIM;
    }

    /// Validate against exact solution at given time
    /// Solution is provided via function pointer.
    /// Plot at a given number of plot points and compute L2 error
    /// and L2 norm of velocity solution over element
    /// Call the broken default
    void compute_error(std::ostream& outfile,
                       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                       const double& time,
                       double& error,
                       double& norm)
    {
      FiniteElement::compute_error(outfile, exact_soln_pt, time, error, norm);
    }

    /// Validate against exact solution.
    /// Solution is provided via function pointer.
    /// Plot at a given number of plot points and compute L2 error
    /// and L2 norm of velocity solution over element
    /// Call the broken default
    void compute_error(std::ostream& outfile,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       double& error,
                       double& norm)
    {
      FiniteElement::compute_error(outfile, exact_soln_pt, error, norm);
    }

    /// Overload the wind function in the advection-diffusion equations.
    /// This provides the coupling from the Navier--Stokes equations to the
    /// advection-diffusion equations because the wind is the fluid velocity.
    void get_wind_adv_diff(const unsigned& ipt,
                           const Vector<double>& s,
                           const Vector<double>& x,
                           Vector<double>& wind) const
    {
      // The wind function is simply the velocity at the points
      this->interpolated_u_nst(s, wind);
    }


    /// Overload the body force in the Navier-Stokes equations
    /// This provides the coupling from the advection-diffusion equations
    /// to the Navier--Stokes equations, the body force is the
    /// temperature multiplied by the Rayleigh number acting in the
    /// direction opposite to gravity.
    void get_body_force_nst(const double& time,
                            const unsigned& ipt,
                            const Vector<double>& s,
                            const Vector<double>& x,
                            Vector<double>& result)
    {
      // Get the temperature
      const double interpolated_t = this->interpolated_u_adv_diff(s);

      // Get vector that indicates the direction of gravity from
      // the Navier-Stokes equations
      Vector<double> gravity(NavierStokesEquations<DIM>::g());

      // Temperature-dependent body force:
      for (unsigned i = 0; i < DIM; i++)
      {
        result[i] = -gravity[i] * interpolated_t * ra();
      }
    } // end of get_body_force

    /// Calculate the element's contribution to the residual vector.
    /// Recall that fill_in_* functions MUST NOT initialise the entries
    /// in the vector to zero. This allows us to call the
    /// fill_in_* functions of the constituent single-physics elements
    /// sequentially, without wiping out any previously computed entries.
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Fill in the residuals of the Navier-Stokes equations
      NavierStokesEquations<DIM>::fill_in_contribution_to_residuals(residuals);

      // Fill in the residuals of the advection-diffusion eqautions
      AdvectionDiffusionEquations<DIM>::fill_in_contribution_to_residuals(
        residuals);
    }


//-----------Finite-difference the entire jacobian-----------------------
//-----------------------------------------------------------------------
#ifdef USE_FD_JACOBIAN_FOR_BUOYANT_Q_ELEMENT


    /// Compute the element's residual vector and the Jacobian matrix.
    /// Jacobian is computed by finite-differencing.
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // This function computes the Jacobian by finite-differencing
      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    /// Add the element's contribution to its residuals vector,
    /// jacobian matrix and mass matrix
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // Call the standard (Broken) function
      // which will prevent these elements from being used
      // in eigenproblems until replaced.
      FiniteElement::fill_in_contribution_to_jacobian_and_mass_matrix(
        residuals, jacobian, mass_matrix);
    }

#else
//--------Finite-difference off-diagonal-blocks in jacobian--------------
//-----------------------------------------------------------------------
#ifdef USE_OFF_DIAGONAL_FD_JACOBIAN_FOR_BUOYANT_Q_ELEMENT

    /// Helper function to get the off-diagonal blocks of the Jacobian
    /// matrix by finite differences
    void fill_in_off_diagonal_jacobian_blocks_by_fd(
      Vector<double>& residuals, DenseMatrix<double>& jacobian)
    {
      // Local storage for the index in the nodes at which the
      // Navier-Stokes velocities are stored (we know that this should be 0,1,2)
      unsigned u_nodal_nst[DIM];
      for (unsigned i = 0; i < DIM; i++)
      {
        u_nodal_nst[i] = this->u_index_nst(i);
      }

      // Local storage for the  index at which the temperature is stored
      unsigned u_nodal_adv_diff = this->u_index_adv_diff();

      // Find the total number of unknowns in the elements
      unsigned n_dof = this->ndof();

      // Temporary storage for residuals
      Vector<double> newres(n_dof);

      // Storage for local unknown and local equation
      int local_unknown = 0, local_eqn = 0;

      // Set the finite difference step
      double fd_step = FiniteElement::Default_fd_jacobian_step;

      // Find the number of nodes
      unsigned n_node = this->nnode();

      // Calculate the contribution of the Navier--Stokes velocities to the
      // advection-diffusion equations

      // Loop over the nodes
      for (unsigned n = 0; n < n_node; n++)
      {
        // There are DIM values of the  velocities
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the local velocity equation number
          local_unknown = this->nodal_local_eqn(n, u_nodal_nst[i]);

          // If it's not pinned
          if (local_unknown >= 0)
          {
            // Get a pointer to the velocity value
            double* value_pt = this->node_pt(n)->value_pt(u_nodal_nst[i]);

            // Save the old value
            double old_var = *value_pt;

            // Increment the value
            *value_pt += fd_step;

            // Get the altered advection-diffusion residuals.
            // which must be done using fill_in_contribution because
            // get_residuals will always return the full residuals
            // because the appropriate fill_in function is overloaded above
            for (unsigned m = 0; m < n_dof; m++)
            {
              newres[m] = 0.0;
            }
            AdvectionDiffusionEquations<DIM>::fill_in_contribution_to_residuals(
              newres);

            // AdvectionDiffusionEquations<DIM>::get_residuals(newres);

            // Now fill in the Advection-Diffusion sub-block
            // of the jacobian
            for (unsigned m = 0; m < n_node; m++)
            {
              // Local equation for temperature
              local_eqn = this->nodal_local_eqn(m, u_nodal_adv_diff);

              // If it's not a boundary condition
              if (local_eqn >= 0)
              {
                double sum =
                  (newres[local_eqn] - residuals[local_eqn]) / fd_step;
                jacobian(local_eqn, local_unknown) = sum;
              }
            }

            // Reset the nodal data
            *value_pt = old_var;
          }
        }


        // Calculate the contribution of the temperature to the Navier--Stokes
        // equations
        {
          // Get the local equation number
          local_unknown = this->nodal_local_eqn(n, u_nodal_adv_diff);

          // If it's not pinned
          if (local_unknown >= 0)
          {
            // Get a pointer to the concentration value
            double* value_pt = this->node_pt(n)->value_pt(u_nodal_adv_diff);

            // Save the old value
            double old_var = *value_pt;

            // Increment the value (Really need access)
            *value_pt += fd_step;

            // Get the altered Navier--Stokes residuals
            // which must be done using fill_in_contribution because
            // get_residuals will always return the full residuals
            // because the appropriate fill_in function is overloaded above
            for (unsigned m = 0; m < n_dof; m++)
            {
              newres[m] = 0.0;
            }
            NavierStokesEquations<DIM>::fill_in_contribution_to_residuals(
              newres);

            // NavierStokesEquations<DIM>::get_residuals(newres);

            // Now fill in the Navier-Stokes sub-block
            for (unsigned m = 0; m < n_node; m++)
            {
              // Loop over the fluid velocities
              for (unsigned j = 0; j < DIM; j++)
              {
                // Local fluid equation
                local_eqn = this->nodal_local_eqn(m, u_nodal_nst[j]);
                if (local_eqn >= 0)
                {
                  double sum =
                    (newres[local_eqn] - residuals[local_eqn]) / fd_step;
                  jacobian(local_eqn, local_unknown) = sum;
                }
              }
            }

            // Reset the nodal data
            *value_pt = old_var;
          }
        }

      } // End of loop over nodes
    }


    /// Compute the element's residual Vector and the Jacobian matrix.
    /// Use finite-differencing only for the off-diagonal blocks.
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Calculate the Navier-Stokes contributions (diagonal block and
      // residuals)
      NavierStokesEquations<DIM>::fill_in_contribution_to_jacobian(residuals,
                                                                   jacobian);

      // Calculate the advection-diffusion contributions
      //(diagonal block and residuals)
      AdvectionDiffusionEquations<DIM>::fill_in_contribution_to_jacobian(
        residuals, jacobian);

      // We now fill in the off-diagonal (interaction) blocks by finite
      // differences.
      fill_in_off_diagonal_jacobian_blocks_by_fd(residuals, jacobian);
    } // End of jacobian calculation


    /// Add the element's contribution to its residuals vector,
    /// jacobian matrix and mass matrix
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // Get the analytic diagonal terms
      NavierStokesEquations<
        DIM>::fill_in_contribution_to_jacobian_and_mass_matrix(residuals,
                                                               jacobian,
                                                               mass_matrix);

      AdvectionDiffusionEquations<
        DIM>::fill_in_contribution_to_jacobian_and_mass_matrix(residuals,
                                                               jacobian,
                                                               mass_matrix);

      // Now fill in the off-diagonal blocks
      fill_in_off_diagonal_jacobian_blocks_by_fd(residuals, jacobian);
    }


    //--------------------Analytic jacobian---------------------------------
//-----------------------------------------------------------------------
#else

    /// Helper function to get the off-diagonal blocks of the Jacobian
    /// matrix analytically
    void fill_in_off_diagonal_jacobian_blocks_analytic(
      Vector<double>& residuals, DenseMatrix<double>& jacobian)
    {
      // We now fill in the off-diagonal (interaction) blocks analytically
      // This requires knowledge of exactly how the residuals are assembled
      // within the parent elements and involves yet another loop over
      // the integration points!

      // Local storage for the index in the nodes at which the
      // Navier-Stokes velocities are stored (we know that this should be 0,1,2)
      unsigned u_nodal_nst[DIM];
      for (unsigned i = 0; i < DIM; i++)
      {
        u_nodal_nst[i] = this->u_index_nst(i);
      }

      // Local storage for the  index at which the temperature is stored
      const unsigned u_nodal_adv_diff = this->u_index_adv_diff();

      // Find out how many nodes there are
      const unsigned n_node = this->nnode();

      // Set up memory for the shape and test functions and their derivatives
      Shape psif(n_node), testf(n_node);
      DShape dpsifdx(n_node, DIM), dtestfdx(n_node, DIM);

      // Number of integration points
      const unsigned n_intpt = this->integral_pt()->nweight();

      // Get Physical Variables from Element
      double Ra = this->ra();
      double Pe = this->pe();
      Vector<double> gravity = this->g();

      // Integers to store the local equations and unknowns
      int local_eqn = 0, local_unknown = 0;

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the integral weight
        double w = this->integral_pt()->weight(ipt);

        // Call the derivatives of the shape and test functions
        double J = this->dshape_and_dtest_eulerian_at_knot_nst(
          ipt, psif, dpsifdx, testf, dtestfdx);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Calculate local values of temperature derivatives
        // Allocate
        Vector<double> interpolated_du_adv_diff_dx(DIM, 0.0);

        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Get the nodal value
          double u_value = this->raw_nodal_value(l, u_nodal_adv_diff);
          // Loop over the derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_du_adv_diff_dx[j] += u_value * dpsifdx(l, j);
          }
        }

        // Assemble the jacobian terms

        // Loop over the test functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Assemble the contributions of the temperature to
          // the Navier--Stokes equations (which arise through the buoyancy
          // body-force term)

          // Loop over the velocity components in the Navier--Stokes equtions
          for (unsigned i = 0; i < DIM; i++)
          {
            // If it's not a boundary condition
            local_eqn = this->nodal_local_eqn(l, u_nodal_nst[i]);
            if (local_eqn >= 0)
            {
              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // We have only the temperature degree of freedom to consider
                // If it's non-zero add in the contribution
                local_unknown = this->nodal_local_eqn(l2, u_nodal_adv_diff);
                if (local_unknown >= 0)
                {
                  // Add contribution to jacobian matrix
                  jacobian(local_eqn, local_unknown) +=
                    -gravity[i] * psif(l2) * Ra * testf(l) * W;
                }
              }
            }
          }

          // Assemble the contributions of the fluid velocities to the
          // advection-diffusion equation for the temperature
          {
            local_eqn = this->nodal_local_eqn(l, u_nodal_adv_diff);
            // If it's not pinned
            if (local_eqn >= 0)
            {
              // Loop over the shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Loop over the velocity degrees of freedom
                for (unsigned i2 = 0; i2 < DIM; i2++)
                {
                  // Get the local unknown
                  local_unknown = this->nodal_local_eqn(l2, u_nodal_nst[i2]);
                  // If it's not pinned
                  if (local_unknown >= 0)
                  {
                    // Add the contribution to the jacobian matrix
                    jacobian(local_eqn, local_unknown) -=
                      Pe * psif(l2) * interpolated_du_adv_diff_dx[i2] *
                      testf(l) * W;
                  }
                }
              }
            }
          }
        }
      }
    }


    /// Compute the element's residual Vector and the Jacobian matrix.
    /// Use analytic expressions for the off-diagonal blocks
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Calculate the Navier-Stokes contributions (diagonal block and
      // residuals)
      NavierStokesEquations<DIM>::fill_in_contribution_to_jacobian(residuals,
                                                                   jacobian);

      // Calculate the advection-diffusion contributions
      //(diagonal block and residuals)
      AdvectionDiffusionEquations<DIM>::fill_in_contribution_to_jacobian(
        residuals, jacobian);

      // Fill in the off diagonal blocks analytically
      fill_in_off_diagonal_jacobian_blocks_analytic(residuals, jacobian);
    }


    /// Add the element's contribution to its residuals vector,
    /// jacobian matrix and mass matrix
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // Get the analytic diagonal terms for the jacobian and mass matrix
      NavierStokesEquations<
        DIM>::fill_in_contribution_to_jacobian_and_mass_matrix(residuals,
                                                               jacobian,
                                                               mass_matrix);

      AdvectionDiffusionEquations<
        DIM>::fill_in_contribution_to_jacobian_and_mass_matrix(residuals,
                                                               jacobian,
                                                               mass_matrix);

      // Now fill in the off-diagonal blocks in the jacobian matrix
      fill_in_off_diagonal_jacobian_blocks_analytic(residuals, jacobian);
    }

#endif
#endif
  };

  //=========================================================================
  /// Set the default physical value to be zero
  //=========================================================================
  template<>
  double BuoyantQCrouzeixRaviartElement<2>::Default_Physical_Constant_Value =
    0.0;
  template<>
  double BuoyantQCrouzeixRaviartElement<3>::Default_Physical_Constant_Value =
    0.0;


  //=======================================================================
  /// Face geometry of the 2D Buoyant Crouzeix_Raviart elements
  //=======================================================================
  template<unsigned int DIM>
  class FaceGeometry<BuoyantQCrouzeixRaviartElement<DIM>>
    : public virtual QElement<DIM - 1, 3>
  {
  public:
    FaceGeometry() : QElement<DIM - 1, 3>() {}
  };

  //=======================================================================
  /// Face geometry of the Face geometry of 2D Buoyant Crouzeix_Raviart elements
  //=======================================================================
  template<>
  class FaceGeometry<FaceGeometry<BuoyantQCrouzeixRaviartElement<2>>>
    : public virtual PointElement
  {
  public:
    FaceGeometry() : PointElement() {}
  };


  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////


  //============start_element_class============================================
  /// A RefineableElement class that solves the
  /// Boussinesq approximation of the Navier--Stokes
  /// and energy equations by coupling two pre-existing classes.
  /// The RefineableQAdvectionDiffusionElement
  /// with bi-quadratic interpolation for the
  /// scalar variable (temperature) and
  /// RefineableQCrouzeixRaviartElement which solves the Navier--Stokes
  /// equations using bi-quadratic interpolation for the velocities and a
  /// discontinuous bi-linear interpolation for the pressure. Note that we are
  /// free to choose the order in which we store the variables at the nodes. In
  /// this case we choose to store the variables in the order fluid velocities
  /// followed by temperature. We must, therefore, overload the function
  /// AdvectionDiffusionEquations<DIM>::u_index_adv_diff() to indicate that
  /// the temperature is stored at the DIM-th position not the 0-th. We do not
  /// need to overload the corresponding function in the
  /// NavierStokesEquations<DIM> class because the velocities are stored
  /// first. Finally, we choose to use the flux-recovery calculation from the
  /// fluid velocities to provide the error used in the mesh adaptation.
  //==========================================================================
  template<unsigned DIM>
  class RefineableBuoyantQCrouzeixRaviartElement
    : public virtual RefineableQAdvectionDiffusionElement<DIM, 3>,
      public virtual RefineableQCrouzeixRaviartElement<DIM>
  {
  private:
    /// Pointer to a new physical variable, the Rayleigh number
    double* Ra_pt;

    /// The static default value of the Rayleigh number
    static double Default_Physical_Constant_Value;

  public:
    /// Constructor: call the underlying constructors and
    /// initialise the pointer to the Rayleigh number to address the default
    /// value of 0.0
    RefineableBuoyantQCrouzeixRaviartElement()
      : RefineableQAdvectionDiffusionElement<DIM, 3>(),
        RefineableQCrouzeixRaviartElement<DIM>()
    {
      Ra_pt = &Default_Physical_Constant_Value;
    }

    /// The required number of values stored at the nodes is
    /// the sum of the required values of the two single-physics elements. This
    /// step is generic for any composed element of this type.
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return (RefineableQAdvectionDiffusionElement<DIM, 3>::required_nvalue(n) +
              RefineableQCrouzeixRaviartElement<DIM>::required_nvalue(n));
    }

    /// Access function for the Rayleigh number (const version)
    const double& ra() const
    {
      return *Ra_pt;
    }

    /// Access function for the pointer to the Rayleigh number
    double*& ra_pt()
    {
      return Ra_pt;
    }


    /// Final override for disable ALE
    void disable_ALE()
    {
      // Disable ALE in both sets of equations
      RefineableNavierStokesEquations<DIM>::disable_ALE();
      RefineableAdvectionDiffusionEquations<DIM>::disable_ALE();
    }

    /// Final override for enable ALE
    void enable_ALE()
    {
      // Enable ALE in both sets of equations
      RefineableNavierStokesEquations<DIM>::enable_ALE();
      RefineableAdvectionDiffusionEquations<DIM>::enable_ALE();
    }


    /// Number of scalars/fields output by this element. Broken
    /// virtual. Needs to be implemented for each new specific element type.
    /// Temporary dummy
    unsigned nscalar_paraview() const
    {
      throw OomphLibError(
        "This function hasn't been implemented for this element",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);

      // Dummy unsigned
      return 0;
    }

    /// Write values of the i-th scalar field at the plot points. Broken
    /// virtual. Needs to be implemented for each new specific element type.
    /// Temporary dummy
    void scalar_value_paraview(std::ofstream& file_out,
                               const unsigned& i,
                               const unsigned& nplot) const
    {
      throw OomphLibError(
        "This function hasn't been implemented for this element",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }


    /// Name of the i-th scalar field. Default implementation
    /// returns V1 for the first one, V2 for the second etc. Can (should!) be
    /// overloaded with more meaningful names.
    std::string scalar_name_paraview(const unsigned& i) const
    {
      return "V" + StringConversion::to_string(i);
    }

    ///  Overload the standard output function with the broken default
    void output(std::ostream& outfile)
    {
      FiniteElement::output(outfile);
    }

    /// Output function:
    ///  x,y,u   or    x,y,z,u at Nplot^DIM plot points
    void output(std::ostream& outfile, const unsigned& nplot)
    {
      // vector of local coordinates
      Vector<double> s(DIM);

      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, nplot, s);

        // Output the position of the plot point
        for (unsigned i = 0; i < DIM; i++)
        {
          outfile << this->interpolated_x(s, i) << " ";
        }

        // Output the fluid velocities at the plot point
        for (unsigned i = 0; i < DIM; i++)
        {
          outfile << this->interpolated_u_nst(s, i) << " ";
        }

        // Output the fluid pressure at the plot point
        outfile << this->interpolated_p_nst(s) << " ";

        // Output the temperature at the plot point
        outfile << this->interpolated_u_adv_diff(s) << " " << std::endl;
      }
      outfile << std::endl;

      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);
    }

    /// C-style output function:  Broken default
    void output(FILE* file_pt)
    {
      FiniteElement::output(file_pt);
    }

    ///  C-style output function: Broken default
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FiniteElement::output(file_pt, n_plot);
    }

    /// Output function for an exact solution: Broken default
    void output_fct(std::ostream& outfile,
                    const unsigned& Nplot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      FiniteElement::output_fct(outfile, Nplot, exact_soln_pt);
    }


    /// Output function for a time-dependent exact solution.
    /// Broken default
    void output_fct(std::ostream& outfile,
                    const unsigned& Nplot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      FiniteElement::output_fct(outfile, Nplot, time, exact_soln_pt);
    }

    /// Overload the index at which the temperature
    /// variable is stored. We choose to store is after the fluid velocities.
    unsigned u_index_adv_diff() const
    {
      return DIM;
    }

    /// Number of vertex nodes in the element is obtained from the
    /// geometric element.
    unsigned nvertex_node() const
    {
      return QElement<DIM, 3>::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element,
    /// Call the geometric element's function.
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QElement<DIM, 3>::vertex_node_pt(j);
    }

    /// The total number of continously interpolated values is
    /// DIM+1 (DIM fluid velocities and one temperature).
    unsigned ncont_interpolated_values() const
    {
      return DIM + 1;
    }


    /// Get the continuously interpolated values at the local coordinate
    /// s. We choose to put the fluid velocities first, followed by the
    /// temperature.
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      // Storage for the fluid velocities
      Vector<double> nst_values;

      // Get the fluid velocities from the fluid element
      RefineableQCrouzeixRaviartElement<DIM>::get_interpolated_values(
        s, nst_values);

      // Storage for the temperature
      Vector<double> advection_values;

      // Get the temperature from the advection-diffusion element
      RefineableQAdvectionDiffusionElement<DIM, 3>::get_interpolated_values(
        s, advection_values);

      // Add the fluid velocities to the values vector
      for (unsigned i = 0; i < DIM; i++)
      {
        values.push_back(nst_values[i]);
      }

      // Add the concentration to the end
      values.push_back(advection_values[0]);
    }


    /// Get all continuously interpolated values at the local
    /// coordinate s at time level t (t=0: present; t>0: previous).
    /// We choose to put the fluid velocities first, followed by the
    /// temperature
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      // Storage for the fluid velocities
      Vector<double> nst_values;

      // Get the fluid velocities from the fluid element
      RefineableQCrouzeixRaviartElement<DIM>::get_interpolated_values(
        t, s, nst_values);

      // Storage for the temperature
      Vector<double> advection_values;

      // Get the temperature from the advection-diffusion element
      RefineableQAdvectionDiffusionElement<DIM, 3>::get_interpolated_values(
        s, advection_values);

      // Add the fluid velocities to the values vector
      for (unsigned i = 0; i < DIM; i++)
      {
        values.push_back(nst_values[i]);
      }

      // Add the concentration to the end
      values.push_back(advection_values[0]);

    } // end of get_interpolated_values


    /// The additional hanging node information must be set up
    /// for both single-physics elements.
    void further_setup_hanging_nodes()
    {
      RefineableQCrouzeixRaviartElement<DIM>::further_setup_hanging_nodes();
      RefineableQAdvectionDiffusionElement<DIM,
                                           3>::further_setup_hanging_nodes();
    }


    /// Call the rebuild_from_sons functions for each of the
    /// constituent multi-physics elements.
    void rebuild_from_sons(Mesh*& mesh_pt)
    {
      RefineableQAdvectionDiffusionElement<DIM, 3>::rebuild_from_sons(mesh_pt);
      RefineableQCrouzeixRaviartElement<DIM>::rebuild_from_sons(mesh_pt);
    }


    /// Call the underlying single-physics element's further_build()
    /// functions and make sure that the pointer to the Rayleigh number
    /// is passed to the sons
    void further_build()
    {
      RefineableQCrouzeixRaviartElement<DIM>::further_build();
      RefineableQAdvectionDiffusionElement<DIM, 3>::further_build();

      // Cast the pointer to the father element to the specific
      // element type
      RefineableBuoyantQCrouzeixRaviartElement<DIM>* cast_father_element_pt =
        dynamic_cast<RefineableBuoyantQCrouzeixRaviartElement<DIM>*>(
          this->father_element_pt());

      // Set the pointer to the Rayleigh number to be the same as that in
      // the father
      this->Ra_pt = cast_father_element_pt->ra_pt();
    } // end of further build


    /// The recovery order is that of the NavierStokes elements.
    unsigned nrecovery_order()
    {
      return RefineableQCrouzeixRaviartElement<DIM>::nrecovery_order();
    }

    /// The number of Z2 flux terms is the same as that in
    /// the fluid element plus that in the advection-diffusion element
    unsigned num_Z2_flux_terms()
    {
      return (
        RefineableQCrouzeixRaviartElement<DIM>::num_Z2_flux_terms() +
        RefineableQAdvectionDiffusionElement<DIM, 3>::num_Z2_flux_terms());
    }


    /// Get the Z2 flux by concatenating the fluxes from the fluid and
    /// the advection diffusion elements.
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
      // Find the number of fluid fluxes
      unsigned n_fluid_flux =
        RefineableQCrouzeixRaviartElement<DIM>::num_Z2_flux_terms();

      // Fill in the first flux entries as the velocity entries
      RefineableQCrouzeixRaviartElement<DIM>::get_Z2_flux(s, flux);

      // Find the number of temperature fluxes
      unsigned n_temp_flux =
        RefineableQAdvectionDiffusionElement<DIM, 3>::num_Z2_flux_terms();
      Vector<double> temp_flux(n_temp_flux);

      // Get the temperature flux
      RefineableQAdvectionDiffusionElement<DIM, 3>::get_Z2_flux(s, temp_flux);

      // Add the temperature flux to the end of the flux vector
      for (unsigned i = 0; i < n_temp_flux; i++)
      {
        flux[n_fluid_flux + i] = temp_flux[i];
      }

    } // end of get_Z2_flux

    /// The number of compound fluxes is two (one for the fluid and
    /// one for the temperature)
    unsigned ncompound_fluxes()
    {
      return 2;
    }

    /// Fill in which flux components are associated with the fluid
    /// measure and which are associated with the temperature measure
    void get_Z2_compound_flux_indices(Vector<unsigned>& flux_index)
    {
      // Find the number of fluid fluxes
      unsigned n_fluid_flux =
        RefineableQCrouzeixRaviartElement<DIM>::num_Z2_flux_terms();
      // Find the number of temperature fluxes
      unsigned n_temp_flux =
        RefineableQAdvectionDiffusionElement<DIM, 3>::num_Z2_flux_terms();

      // The fluid fluxes are first
      // The values of the flux_index vector are zero on entry, so we
      // could omit this line
      for (unsigned i = 0; i < n_fluid_flux; i++)
      {
        flux_index[i] = 0;
      }

      // Set the temperature fluxes (the last set of fluxes
      for (unsigned i = 0; i < n_temp_flux; i++)
      {
        flux_index[n_fluid_flux + i] = 1;
      }

    } // end of get_Z2_compound_flux_indices


    /// Validate against exact solution at given time
    /// Solution is provided via function pointer.
    /// Plot at a given number of plot points and compute L2 error
    /// and L2 norm of velocity solution over element
    /// Overload to broken default
    void compute_error(std::ostream& outfile,
                       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                       const double& time,
                       double& error,
                       double& norm)
    {
      FiniteElement::compute_error(outfile, exact_soln_pt, time, error, norm);
    }

    /// Validate against exact solution.
    /// Solution is provided via function pointer.
    /// Plot at a given number of plot points and compute L2 error
    /// and L2 norm of velocity solution over element
    /// Overload to broken default.
    void compute_error(std::ostream& outfile,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       double& error,
                       double& norm)
    {
      FiniteElement::compute_error(outfile, exact_soln_pt, error, norm);
    }

    /// Overload the wind function in the advection-diffusion equations.
    /// This provides the coupling from the Navier--Stokes equations to the
    /// advection-diffusion equations because the wind is the fluid velocity.
    void get_wind_adv_diff(const unsigned& ipt,
                           const Vector<double>& s,
                           const Vector<double>& x,
                           Vector<double>& wind) const
    {
      // The wind function is simply the velocity at the points
      this->interpolated_u_nst(s, wind);
    }


    /// Overload the body force in the navier-stokes equations
    /// This provides the coupling from the advection-diffusion equations
    /// to the Navier--Stokes equations, the body force is the
    /// temperature multiplied by the Rayleigh number acting in the
    /// direction opposite to gravity.
    void get_body_force_nst(const double& time,
                            const unsigned& ipt,
                            const Vector<double>& s,
                            const Vector<double>& x,
                            Vector<double>& result)
    {
      // Get the temperature
      const double interpolated_t = this->interpolated_u_adv_diff(s);

      // Get vector that indicates the direction of gravity from
      // the Navier-Stokes equations
      Vector<double> gravity(NavierStokesEquations<DIM>::g());

      // Temperature-dependent body force:
      for (unsigned i = 0; i < DIM; i++)
      {
        result[i] = -gravity[i] * interpolated_t * ra();
      }
    }

    /// Fill in the constituent elements' contribution to the residual vector.
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the residuals of the Navier-Stokes equations
      RefineableNavierStokesEquations<DIM>::fill_in_contribution_to_residuals(
        residuals);

      // Call the residuals of the advection-diffusion equations
      RefineableAdvectionDiffusionEquations<
        DIM>::fill_in_contribution_to_residuals(residuals);
    }


    /// Compute the element's residual Vector and the jacobian matrix
    /// using full finite differences, the default implementation
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
#ifdef USE_FD_JACOBIAN_FOR_REFINEABLE_BUOYANT_Q_ELEMENT
      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
#else
      // Calculate the Navier-Stokes contributions (diagonal block and
      // residuals)
      RefineableNavierStokesEquations<DIM>::fill_in_contribution_to_jacobian(
        residuals, jacobian);

      // Calculate the advection-diffusion contributions
      //(diagonal block and residuals)
      RefineableAdvectionDiffusionEquations<
        DIM>::fill_in_contribution_to_jacobian(residuals, jacobian);

      // We now fill in the off-diagonal (interaction) blocks analytically
      this->fill_in_off_diagonal_jacobian_blocks_analytic(residuals, jacobian);
#endif
    } // End of jacobian calculation

    /// Add the element's contribution to its residuals vector,
    /// jacobian matrix and mass matrix
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // Call the (broken) version in the base class
      FiniteElement::fill_in_contribution_to_jacobian_and_mass_matrix(
        residuals, jacobian, mass_matrix);
    }

    /// Compute the contribution of the off-diagonal blocks
    /// analytically.
    void fill_in_off_diagonal_jacobian_blocks_analytic(
      Vector<double>& residuals, DenseMatrix<double>& jacobian)
    {
      // Perform another loop over the integration loops using the information
      // from the original elements' residual assembly loops to determine
      // the conributions to the jacobian

      // Local storage for pointers to hang_info objects
      HangInfo *hang_info_pt = 0, *hang_info2_pt = 0;

      // Local storage for the index in the nodes at which the
      // Navier-Stokes velocities are stored (we know that this should be 0,1,2)
      unsigned u_nodal_nst[DIM];
      for (unsigned i = 0; i < DIM; i++)
      {
        u_nodal_nst[i] = this->u_index_nst(i);
      }

      // Local storage for the  index at which the temperature is stored
      const unsigned u_nodal_adv_diff = this->u_index_adv_diff();

      // Find out how many nodes there are
      const unsigned n_node = this->nnode();

      // Set up memory for the shape and test functions and their derivatives
      Shape psif(n_node), testf(n_node);
      DShape dpsifdx(n_node, DIM), dtestfdx(n_node, DIM);

      // Number of integration points
      const unsigned n_intpt = this->integral_pt()->nweight();

      // Get Physical Variables from Element
      double Ra = this->ra();
      double Pe = this->pe();
      Vector<double> gravity = this->g();

      // Integers to store the local equations and unknowns
      int local_eqn = 0, local_unknown = 0;

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the integral weight
        double w = this->integral_pt()->weight(ipt);

        // Call the derivatives of the shape and test functions
        double J = this->dshape_and_dtest_eulerian_at_knot_nst(
          ipt, psif, dpsifdx, testf, dtestfdx);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Calculate local values of temperature derivatives
        // Allocate
        Vector<double> interpolated_du_adv_diff_dx(DIM, 0.0);

        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Get the nodal value
          double u_value = this->nodal_value(l, u_nodal_adv_diff);
          // Loop over the derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_du_adv_diff_dx[j] += u_value * dpsifdx(l, j);
          }
        }

        // Assemble the Jacobian terms
        //---------------------------

        // Loop over the test functions/eqns
        for (unsigned l = 0; l < n_node; l++)
        {
          // Local variables to store the number of master nodes and
          // the weight associated with the shape function if the node is
          // hanging
          unsigned n_master = 1;
          double hang_weight = 1.0;

          // Local bool (is the node hanging)
          bool is_node_hanging = this->node_pt(l)->is_hanging();

          // If the node is hanging, get the number of master nodes
          if (is_node_hanging)
          {
            hang_info_pt = this->node_pt(l)->hanging_pt();
            n_master = hang_info_pt->nmaster();
          }
          // Otherwise there is just one master node, the node itself
          else
          {
            n_master = 1;
          }

          // Loop over the master nodes
          for (unsigned m = 0; m < n_master; m++)
          {
            // If the node is hanging get weight from master node
            if (is_node_hanging)
            {
              // Get the hang weight from the master node
              hang_weight = hang_info_pt->master_weight(m);
            }
            else
            {
              // Node contributes with full weight
              hang_weight = 1.0;
            }


            // Assemble derivatives of Navier Stokes momentum w.r.t. temperature
            //-----------------------------------------------------------------

            // Loop over velocity components for equations
            for (unsigned i = 0; i < DIM; i++)
            {
              // Get the equation number
              if (is_node_hanging)
              {
                // Get the equation number from the master node
                local_eqn = this->local_hang_eqn(
                  hang_info_pt->master_node_pt(m), u_nodal_nst[i]);
              }
              else
              {
                // Local equation number
                local_eqn = this->nodal_local_eqn(l, u_nodal_nst[i]);
              }

              if (local_eqn >= 0)
              {
                // Local variables to store the number of master nodes
                // and the weights associated with each hanging node
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;

                // Loop over the nodes for the unknowns
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Local bool (is the node hanging)
                  bool is_node2_hanging = this->node_pt(l2)->is_hanging();

                  // If the node is hanging, get the number of master nodes
                  if (is_node2_hanging)
                  {
                    hang_info2_pt = this->node_pt(l2)->hanging_pt();
                    n_master2 = hang_info2_pt->nmaster();
                  }
                  // Otherwise there is one master node, the node itself
                  else
                  {
                    n_master2 = 1;
                  }

                  // Loop over the master nodes
                  for (unsigned m2 = 0; m2 < n_master2; m2++)
                  {
                    if (is_node2_hanging)
                    {
                      // Read out the local unknown from the master node
                      local_unknown = this->local_hang_eqn(
                        hang_info2_pt->master_node_pt(m2), u_nodal_adv_diff);
                      // Read out the hanging weight from the master node
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    else
                    {
                      // The local unknown number comes from the node
                      local_unknown =
                        this->nodal_local_eqn(l2, u_nodal_adv_diff);
                      // The hang weight is one
                      hang_weight2 = 1.0;
                    }

                    if (local_unknown >= 0)
                    {
                      // Add contribution to jacobian matrix
                      jacobian(local_eqn, local_unknown) +=
                        -gravity[i] * psif(l2) * Ra * testf(l) * W *
                        hang_weight * hang_weight2;
                    }
                  }
                }
              }
            }


            // Assemble derivative of adv diff eqn w.r.t. fluid veloc
            //------------------------------------------------------
            {
              // Get the equation number
              if (is_node_hanging)
              {
                // Get the equation number from the master node
                local_eqn = this->local_hang_eqn(
                  hang_info_pt->master_node_pt(m), u_nodal_adv_diff);
              }
              else
              {
                // Local equation number
                local_eqn = this->nodal_local_eqn(l, u_nodal_adv_diff);
              }

              // If it's not pinned
              if (local_eqn >= 0)
              {
                // Local variables to store the number of master nodes
                // and the weights associated with each hanging node
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;

                // Loop over the nodes for the unknowns
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Local bool (is the node hanging)
                  bool is_node2_hanging = this->node_pt(l2)->is_hanging();

                  // If the node is hanging, get the number of master nodes
                  if (is_node2_hanging)
                  {
                    hang_info2_pt = this->node_pt(l2)->hanging_pt();
                    n_master2 = hang_info2_pt->nmaster();
                  }
                  // Otherwise there is one master node, the node itself
                  else
                  {
                    n_master2 = 1;
                  }

                  // Loop over the master nodes
                  for (unsigned m2 = 0; m2 < n_master2; m2++)
                  {
                    // If the node is hanging
                    if (is_node2_hanging)
                    {
                      // Read out the hanging weight from the master node
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    // If the node is not hanging
                    else
                    {
                      // The hang weight is one
                      hang_weight2 = 1.0;
                    }

                    // Loop over the velocity degrees of freedom
                    for (unsigned i2 = 0; i2 < DIM; i2++)
                    {
                      // If the node is hanging
                      if (is_node2_hanging)
                      {
                        // Read out the local unknown from the master node
                        local_unknown = this->local_hang_eqn(
                          hang_info2_pt->master_node_pt(m2), u_nodal_nst[i2]);
                      }
                      else
                      {
                        // The local unknown number comes from the node
                        local_unknown =
                          this->nodal_local_eqn(l2, u_nodal_nst[i2]);
                      }

                      // If it's not pinned
                      if (local_unknown >= 0)
                      {
                        // Add the contribution to the jacobian matrix
                        jacobian(local_eqn, local_unknown) -=
                          Pe * psif(l2) * interpolated_du_adv_diff_dx[i2] *
                          testf(l) * W * hang_weight * hang_weight2;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    } // End of function
  };


  //===================================================================
  /// Set the default value of the Rayleigh number to be zero
  //===================================================================
  template<>
  double RefineableBuoyantQCrouzeixRaviartElement<
    2>::Default_Physical_Constant_Value = 0.0;

  template<>
  double RefineableBuoyantQCrouzeixRaviartElement<
    3>::Default_Physical_Constant_Value = 0.0;


} // namespace oomph

#endif
