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
// Header file for Advection Diffusion elements
#ifndef OOMPH_GEN_ADV_DIFF_ELEMENTS_HEADER
#define OOMPH_GEN_ADV_DIFF_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// OOMPH-LIB headers
#include "../generic/nodes.h"
#include "../generic/Qelements.h"
#include "../generic/oomph_utilities.h"

namespace oomph
{
  //=============================================================
  /// \short A class for all elements that solve the Advection
  /// Diffusion equations in conservative form using isoparametric elements.
  /// \f[
  /// \frac{\partial}{\partial x_{i}}\left(
  /// Pe w_{i}(x_{k}) u - D_{ij}(x_{k})\frac{\partial u}{\partial x_{j}}\right)
  /// = f(x_{j})
  /// \f]
  /// This contains the generic maths. Shape functions, geometric
  /// mapping etc. must get implemented in derived class.
  //=============================================================
  template<unsigned DIM>
  class GeneralisedAdvectionDiffusionEquations : public virtual FiniteElement
  {
  public:
    /// \short Function pointer to source function fct(x,f(x)) --
    /// x is a Vector!
    typedef void (*GeneralisedAdvectionDiffusionSourceFctPt)(
      const Vector<double> &x, double &f);

    /// \short Function pointer to wind function fct(x,w(x)) --
    /// x is a Vector!
    typedef void (*GeneralisedAdvectionDiffusionWindFctPt)(
      const Vector<double> &x, Vector<double> &wind);

    /// \short Funciton pointer to a diffusivity function
    typedef void (*GeneralisedAdvectionDiffusionDiffFctPt)(
      const Vector<double> &x, DenseMatrix<double> &D);

    /// \short Constructor: Initialise the Source_fct_pt and Wind_fct_pt
    /// to null and set (pointer to) Peclet number to default
    GeneralisedAdvectionDiffusionEquations() :
      Source_fct_pt(0),
      Wind_fct_pt(0),
      Conserved_wind_fct_pt(0),
      Diff_fct_pt(0),
      ALE_is_disabled(false)
    {
      // Set Peclet number to default
      Pe_pt = &Default_peclet_number;
      // Set Peclet Strouhal number to default
      PeSt_pt = &Default_peclet_number;
    }

    /// Broken copy constructor
    GeneralisedAdvectionDiffusionEquations(
      const GeneralisedAdvectionDiffusionEquations &dummy)
    {
      BrokenCopy::broken_copy("GeneralisedAdvectionDiffusionEquations");
    }

    /// Broken assignment operator
    void operator=(const GeneralisedAdvectionDiffusionEquations &)
    {
      BrokenCopy::broken_assign("GeneralisedAdvectionDiffusionEquations");
    }

    /// \short Return the index at which the unknown value
    /// is stored. The default value, 0, is appropriate for single-physics
    /// problems, when there is only one variable, the value that satisfies
    /// the advection-diffusion equation.
    /// In derived multi-physics elements, this function should be overloaded
    /// to reflect the chosen storage scheme. Note that these equations require
    /// that the unknown is always stored at the same index at each node.
    virtual inline unsigned u_index_cons_adv_diff() const
    {
      return 0;
    }

    /// \short du/dt at local node n.
    /// Uses suitably interpolated value for hanging nodes.
    double du_dt_cons_adv_diff(const unsigned &n) const
    {
      // Get the data's timestepper
      TimeStepper *time_stepper_pt = this->node_pt(n)->time_stepper_pt();

      // Initialise dudt
      double dudt = 0.0;
      // Loop over the timesteps, if there is a non Steady timestepper
      if (!time_stepper_pt->is_steady())
      {
        // Find the index at which the variable is stored
        const unsigned u_nodal_index = u_index_cons_adv_diff();

        // Number of timsteps (past & present)
        const unsigned n_time = time_stepper_pt->ntstorage();

        for (unsigned t = 0; t < n_time; t++)
        {
          dudt +=
            time_stepper_pt->weight(1, t) * nodal_value(t, n, u_nodal_index);
        }
      }
      return dudt;
    }

    /// \short Disable ALE, i.e. assert the mesh is not moving -- you do this
    /// at your own risk!
    void disable_ALE()
    {
      ALE_is_disabled = true;
    }

    /// \short (Re-)enable ALE, i.e. take possible mesh motion into account
    /// when evaluating the time-derivative. Note: By default, ALE is
    /// enabled, at the expense of possibly creating unnecessary work
    /// in problems where the mesh is, in fact, stationary.
    void enable_ALE()
    {
      ALE_is_disabled = false;
    }

    /// Output with default number of plot points
    void output(std::ostream &outfile)
    {
      unsigned nplot = 5;
      output(outfile, nplot);
    }

    /// \short Output FE representation of soln: x,y,u or x,y,z,u at
    /// nplot^DIM plot points
    void output(std::ostream &outfile, const unsigned &nplot);

    /// C_style output with default number of plot points
    void output(FILE *file_pt)
    {
      unsigned n_plot = 5;
      output(file_pt, n_plot);
    }

    /// \short C-style output FE representation of soln: x,y,u or x,y,z,u at
    /// n_plot^DIM plot points
    void output(FILE *file_pt, const unsigned &n_plot);

    /// Output exact soln: x,y,u_exact or x,y,z,u_exact at nplot^DIM plot points
    void output_fct(std::ostream &outfile,
                    const unsigned &nplot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

    /// \short Output exact soln: x,y,u_exact or x,y,z,u_exact at
    /// nplot^DIM plot points (dummy time-dependent version to
    /// keep intel compiler happy)
    virtual void output_fct(
      std::ostream &outfile,
      const unsigned &nplot,
      const double &time,
      FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      throw OomphLibError("There is no time-dependent output_fct() for "
                          "Advection Diffusion elements",
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
        "No time-dependent compute_error() for Advection Diffusion elements",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// \short Integrate the concentration over the element
    double integrate_u();

    /// Access function: Pointer to source function
    GeneralisedAdvectionDiffusionSourceFctPt &source_fct_pt()
    {
      return Source_fct_pt;
    }

    /// Access function: Pointer to source function. Const version
    GeneralisedAdvectionDiffusionSourceFctPt source_fct_pt() const
    {
      return Source_fct_pt;
    }

    /// Access function: Pointer to wind function
    GeneralisedAdvectionDiffusionWindFctPt &wind_fct_pt()
    {
      return Wind_fct_pt;
    }

    /// Access function: Pointer to wind function. Const version
    GeneralisedAdvectionDiffusionWindFctPt wind_fct_pt() const
    {
      return Wind_fct_pt;
    }

    /// Access function: Pointer to additional (conservative) wind function
    GeneralisedAdvectionDiffusionWindFctPt &conserved_wind_fct_pt()
    {
      return Conserved_wind_fct_pt;
    }

    /// \short Access function: Pointer to additoinal (conservative)
    /// wind function.
    /// Const version
    GeneralisedAdvectionDiffusionWindFctPt conserved_wind_fct_pt() const
    {
      return Conserved_wind_fct_pt;
    }

    /// Access function: Pointer to diffusion  function
    GeneralisedAdvectionDiffusionDiffFctPt &diff_fct_pt()
    {
      return Diff_fct_pt;
    }

    /// Access function: Pointer to diffusion function. Const version
    GeneralisedAdvectionDiffusionDiffFctPt diff_fct_pt() const
    {
      return Diff_fct_pt;
    }

    /// Peclet number
    const double &pe() const
    {
      return *Pe_pt;
    }

    /// Pointer to Peclet number
    double *&pe_pt()
    {
      return Pe_pt;
    }

    /// Peclet number multiplied by Strouhal number
    const double &pe_st() const
    {
      return *PeSt_pt;
    }

    /// Pointer to Peclet number multipled by Strouha number
    double *&pe_st_pt()
    {
      return PeSt_pt;
    }

    /// \short Get source term at (Eulerian) position x. This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the strength of the source function might be determined by
    /// another system of equations
    inline virtual void get_source_cons_adv_diff(const unsigned &ipt,
                                                 const Vector<double> &x,
                                                 double &source) const
    {
      // If no source function has been set, return zero
      if (Source_fct_pt == 0)
      {
        source = 0.0;
      }
      else
      {
        // Get source strength
        (*Source_fct_pt)(x, source);
      }
    }

    /// \short Get wind at (Eulerian) position x and/or local coordinate s.
    /// This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the wind function might be determined by
    /// another system of equations
    inline virtual void get_wind_cons_adv_diff(const unsigned &ipt,
                                               const Vector<double> &s,
                                               const Vector<double> &x,
                                               Vector<double> &wind) const
    {
      // If no wind function has been set, return zero
      if (Wind_fct_pt == 0)
      {
        for (unsigned i = 0; i < DIM; i++)
        {
          wind[i] = 0.0;
        }
      }
      else
      {
        // Get wind
        (*Wind_fct_pt)(x, wind);
      }
    }

    /// \short Get additional (conservative)
    /// wind at (Eulerian) position x and/or local coordinate s.
    /// This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the wind function might be determined by
    /// another system of equations
    inline virtual void get_conserved_wind_cons_adv_diff(
      const unsigned &ipt,
      const Vector<double> &s,
      const Vector<double> &x,
      Vector<double> &wind) const
    {
      // If no wind function has been set, return zero
      if (Conserved_wind_fct_pt == 0)
      {
        for (unsigned i = 0; i < DIM; i++)
        {
          wind[i] = 0.0;
        }
      }
      else
      {
        // Get wind
        (*Conserved_wind_fct_pt)(x, wind);
      }
    }

    /// \short Get diffusivity tensor at (Eulerian) position
    /// x and/or local coordinate s.
    /// This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the wind function might be determined by
    /// another system of equations
    inline virtual void get_diff_cons_adv_diff(const unsigned &ipt,
                                               const Vector<double> &s,
                                               const Vector<double> &x,
                                               DenseMatrix<double> &D) const
    {
      // If no wind function has been set, return identity
      if (Diff_fct_pt == 0)
      {
        for (unsigned i = 0; i < DIM; i++)
        {
          for (unsigned j = 0; j < DIM; j++)
          {
            if (i == j)
            {
              D(i, j) = 1.0;
            }
            else
            {
              D(i, j) = 0.0;
            }
          }
        }
      }
      else
      {
        // Get diffusivity tensor
        (*Diff_fct_pt)(x, D);
      }
    }

    /// Get flux: \f$\mbox{flux}[i] = \mbox{d}u / \mbox{d}x_i \f$
    void get_flux(const Vector<double> &s, Vector<double> &flux) const
    {
      // Find out how many nodes there are in the element
      unsigned n_node = nnode();

      // Get the nodal index at which the unknown is stored
      unsigned u_nodal_index = u_index_cons_adv_diff();

      // Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidx(n_node, DIM);

      // Call the derivatives of the shape and test functions
      dshape_eulerian(s, psi, dpsidx);

      // Initialise to zero
      for (unsigned j = 0; j < DIM; j++)
      {
        flux[j] = 0.0;
      }

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over derivative directions
        for (unsigned j = 0; j < DIM; j++)
        {
          flux[j] += nodal_value(l, u_nodal_index) * dpsidx(l, j);
        }
      }
    }

    /// Get flux: \f$\mbox{flux}[i] = \mbox{d}u / \mbox{d}x_i \f$
    void get_total_flux(const Vector<double> &s,
                        Vector<double> &total_flux) const
    {
      // Find out how many nodes there are in the element
      const unsigned n_node = nnode();

      // Get the nodal index at which the unknown is stored
      const unsigned u_nodal_index = u_index_cons_adv_diff();

      // Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidx(n_node, DIM);

      // Call the derivatives of the shape and test functions
      dshape_eulerian(s, psi, dpsidx);

      // Storage for the Eulerian position
      Vector<double> interpolated_x(DIM, 0.0);
      // Storage for the concentration
      double interpolated_u = 0.0;
      // Storage for the derivatives of the concentration
      Vector<double> interpolated_dudx(DIM, 0.0);

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the value at the node
        const double u_value = this->nodal_value(l, u_nodal_index);
        interpolated_u += u_value * psi(l);
        // Loop over directions
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_x[j] += this->nodal_position(l, j) * psi(l);
          interpolated_dudx[j] += u_value * dpsidx(l, j);
        }
      }

      // Dummy integration point
      unsigned ipt = 0;

      // Get the conserved wind (non-divergence free)
      Vector<double> conserved_wind(DIM);
      get_conserved_wind_cons_adv_diff(ipt, s, interpolated_x, conserved_wind);

      // Get diffusivity tensor
      DenseMatrix<double> D(DIM, DIM);
      get_diff_cons_adv_diff(ipt, s, interpolated_x, D);

      // Calculate the total flux made up of the diffusive flux
      // and the conserved wind
      for (unsigned i = 0; i < DIM; i++)
      {
        total_flux[i] = 0.0;
        for (unsigned j = 0; j < DIM; j++)
        {
          total_flux[i] += D(i, j) * interpolated_dudx[j];
        }
        total_flux[i] -= conserved_wind[i] * interpolated_u;
      }
    }

    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      // Call the generic residuals function with flag set to 0 and using
      // a dummy matrix
      fill_in_generic_residual_contribution_cons_adv_diff(
        residuals,
        GeneralisedElement::Dummy_matrix,
        GeneralisedElement::Dummy_matrix,
        0);
    }

    /// \short Add the element's contribution to its residual vector and
    /// the element Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
    {
      // Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_cons_adv_diff(
        residuals, jacobian, GeneralisedElement::Dummy_matrix, 1);
    }

    /// Add the element's contribution to its residuals vector,
    /// jacobian matrix and mass matrix
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double> &residuals,
      DenseMatrix<double> &jacobian,
      DenseMatrix<double> &mass_matrix)
    {
      // Call the generic routine with the flag set to 2
      fill_in_generic_residual_contribution_cons_adv_diff(
        residuals, jacobian, mass_matrix, 2);
    }

    /// Return FE representation of function value u(s) at local coordinate s
    inline double interpolated_u_cons_adv_diff(const Vector<double> &s) const
    {
      // Find number of nodes
      unsigned n_node = nnode();

      // Get the nodal index at which the unknown is stored
      unsigned u_nodal_index = u_index_cons_adv_diff();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise value of u
      double interpolated_u = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += nodal_value(l, u_nodal_index) * psi[l];
      }

      return (interpolated_u);
    }

    /// \short Self-test: Return 0 for OK
    unsigned self_test();

  protected:
    /// \short Shape/test functions and derivs w.r.t. to global coords at
    /// local coord. s; return  Jacobian of mapping
    virtual double dshape_and_dtest_eulerian_cons_adv_diff(
      const Vector<double> &s,
      Shape &psi,
      DShape &dpsidx,
      Shape &test,
      DShape &dtestdx) const = 0;

    /// \short Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return  Jacobian of mapping
    virtual double dshape_and_dtest_eulerian_at_knot_cons_adv_diff(
      const unsigned &ipt,
      Shape &psi,
      DShape &dpsidx,
      Shape &test,
      DShape &dtestdx) const = 0;

    /// \short Add the element's contribution to its residual vector only
    /// (if flag=and/or element  Jacobian matrix
    virtual void fill_in_generic_residual_contribution_cons_adv_diff(
      Vector<double> &residuals,
      DenseMatrix<double> &jacobian,
      DenseMatrix<double> &mass_matrix,
      unsigned flag);

    /// Pointer to global Peclet number
    double *Pe_pt;

    /// Pointer to global Peclet number multiplied by Strouhal number
    double *PeSt_pt;

    /// Pointer to source function:
    GeneralisedAdvectionDiffusionSourceFctPt Source_fct_pt;

    /// Pointer to wind function:
    GeneralisedAdvectionDiffusionWindFctPt Wind_fct_pt;

    /// Pointer to additional (conservative) wind function:
    GeneralisedAdvectionDiffusionWindFctPt Conserved_wind_fct_pt;

    /// Pointer to diffusivity funciton
    GeneralisedAdvectionDiffusionDiffFctPt Diff_fct_pt;

    /// \short Boolean flag to indicate if ALE formulation is disabled when
    /// time-derivatives are computed. Only set to false if you're sure
    /// that the mesh is stationary.
    bool ALE_is_disabled;

  private:
    /// Static default value for the Peclet number
    static double Default_peclet_number;
  };

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  //======================================================================
  /// \short QGeneralisedAdvectionDiffusionElement elements are
  /// linear/quadrilateral/brick-shaped Advection Diffusion elements with
  /// isoparametric interpolation for the function.
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class QGeneralisedAdvectionDiffusionElement :
    public virtual QElement<DIM, NNODE_1D>,
    public virtual GeneralisedAdvectionDiffusionEquations<DIM>
  {
  private:
    /// \short Static array of ints to hold number of variables at
    /// nodes: Initial_Nvalue[n]
    static const unsigned Initial_Nvalue;

  public:
    ///\short  Constructor: Call constructors for QElement and
    /// Advection Diffusion equations
    QGeneralisedAdvectionDiffusionElement() :
      QElement<DIM, NNODE_1D>(), GeneralisedAdvectionDiffusionEquations<DIM>()
    {
    }

    /// Broken copy constructor
    QGeneralisedAdvectionDiffusionElement(
      const QGeneralisedAdvectionDiffusionElement<DIM, NNODE_1D> &dummy)
    {
      BrokenCopy::broken_copy("QGeneralisedAdvectionDiffusionElement");
    }

    /// Broken assignment operator
    void operator=(const QGeneralisedAdvectionDiffusionElement<DIM, NNODE_1D> &)
    {
      BrokenCopy::broken_assign("QGeneralisedAdvectionDiffusionElement");
    }

    /// \short  Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned &n) const
    {
      return Initial_Nvalue;
    }

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream &outfile)
    {
      GeneralisedAdvectionDiffusionEquations<DIM>::output(outfile);
    }

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot)
    {
      GeneralisedAdvectionDiffusionEquations<DIM>::output(outfile, n_plot);
    }

    /// \short C-style output function:
    ///  x,y,u   or    x,y,z,u
    void output(FILE *file_pt)
    {
      GeneralisedAdvectionDiffusionEquations<DIM>::output(file_pt);
    }

    ///  \short C-style output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(FILE *file_pt, const unsigned &n_plot)
    {
      GeneralisedAdvectionDiffusionEquations<DIM>::output(file_pt, n_plot);
    }

    /// \short Output function for an exact solution:
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
    void output_fct(std::ostream &outfile,
                    const unsigned &n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      GeneralisedAdvectionDiffusionEquations<DIM>::output_fct(
        outfile, n_plot, exact_soln_pt);
    }

    /// \short Output function for a time-dependent exact solution.
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
    /// (Calls the steady version)
    void output_fct(std::ostream &outfile,
                    const unsigned &n_plot,
                    const double &time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      GeneralisedAdvectionDiffusionEquations<DIM>::output_fct(
        outfile, n_plot, time, exact_soln_pt);
    }

  protected:
    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_cons_adv_diff(
      const Vector<double> &s,
      Shape &psi,
      DShape &dpsidx,
      Shape &test,
      DShape &dtestdx) const;

    /// \short Shape, test functions & derivs. w.r.t. to global coords. at
    /// integration point ipt. Return Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_cons_adv_diff(
      const unsigned &ipt,
      Shape &psi,
      DShape &dpsidx,
      Shape &test,
      DShape &dtestdx) const;
  };

  // Inline functions:

  //======================================================================
  /// \short Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QGeneralisedAdvectionDiffusionElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_cons_adv_diff(const Vector<double> &s,
                                            Shape &psi,
                                            DShape &dpsidx,
                                            Shape &test,
                                            DShape &dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian(s, psi, dpsidx);

    // Loop over the test functions and derivatives and set them equal to the
    // shape functions
    for (unsigned i = 0; i < NNODE_1D; i++)
    {
      test[i] = psi[i];
      for (unsigned j = 0; j < DIM; j++)
      {
        dtestdx(i, j) = dpsidx(i, j);
      }
    }

    // Return the jacobian
    return J;
  }

  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QGeneralisedAdvectionDiffusionElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_cons_adv_diff(const unsigned &ipt,
                                                    Shape &psi,
                                                    DShape &dpsidx,
                                                    Shape &test,
                                                    DShape &dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

    // Set the test functions equal to the shape functions (pointer copy)
    test = psi;
    dtestdx = dpsidx;

    // Return the jacobian
    return J;
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //=======================================================================
  /// \short Face geometry for the QGeneralisedAdvectionDiffusionElement
  /// elements: The spatial dimension of the face elements is one lower than
  /// that of the bulk element but they have the same number of points along
  /// their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QGeneralisedAdvectionDiffusionElement<DIM, NNODE_1D>> :
    public virtual QElement<DIM - 1, NNODE_1D>
  {
  public:
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : QElement<DIM - 1, NNODE_1D>() {}
  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //=======================================================================
  /// Face geometry for the 1D QGeneralisedAdvectionDiffusion elements: Point
  /// elements
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<QGeneralisedAdvectionDiffusionElement<1, NNODE_1D>> :
    public virtual PointElement
  {
  public:
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : PointElement() {}
  };

} // namespace oomph

#endif
