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
// Header file for axisymmetric FoepplvonKarman elements
#ifndef OOMPH_AXISYM_DISPL_BASED_FOEPPLVONKARMAN_ELEMENTS_HEADER
#define OOMPH_AXISYM_DISPL_BASED_FOEPPLVONKARMAN_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// OOMPH-LIB headers
#include "generic/nodes.h"
#include "generic/Qelements.h"
#include "generic/oomph_utilities.h"

namespace oomph
{
  //=============================================================
  /// A class for all isoparametric elements that solve the
  /// axisYm Foeppl von Karman equations in a displacement based formulation.
  ///
  /// This contains the generic maths. Shape functions, geometric
  /// mapping etc. must get implemented in derived class.
  //=============================================================
  class AxisymFoepplvonKarmanEquations : public virtual FiniteElement
  {
  public:
    /// \short Function pointer to pressure function fct(r,f(r)) --
    /// r is a Vector!
    typedef void (*AxisymFoepplvonKarmanPressureFctPt)(const double& r,
                                                       double& f);

    /// \short Constructor (must initialise the Pressure_fct_pt). Also
    /// set physical parameters to their default values.
    AxisymFoepplvonKarmanEquations() : Pressure_fct_pt(0), Nu_pt(0) {}

    /// Broken copy constructor
    AxisymFoepplvonKarmanEquations(const AxisymFoepplvonKarmanEquations& dummy)
    {
      BrokenCopy::broken_copy("AxisymFoepplvonKarmanEquations");
    }

    /// Broken assignment operator
    void operator=(const AxisymFoepplvonKarmanEquations&)
    {
      BrokenCopy::broken_assign("AxisymFoepplvonKarmanEquations");
    }

    /// Poisson's ratio
    const double& nu() const
    {
#ifdef PARANOID
      if (Nu_pt == 0)
      {
        std::stringstream error_stream;
        error_stream << "Nu has not yet been set!" << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return *Nu_pt;
    }

    /// Pointer to Poisson's ratio
    double*& nu_pt()
    {
      return Nu_pt;
    }

    /// FvK parameter
    const double& eta() const
    {
      return *Eta_pt;
    }

    /// Pointer to FvK parameter
    double*& eta_pt()
    {
      return Eta_pt;
    }

    /// \short Return the index at which the i-th unknown value
    /// is stored. The default value, i, is appropriate for single-physics
    /// problems. By default, these are:
    /// 0: transverse displacement w
    /// 1: laplacian w
    /// 2: radial displacement u
    /// In derived multi-physics elements, this function should be overloaded
    /// to reflect the chosen storage scheme. Note that these equations require
    /// that the unknown is always stored at the same index at each node.
    virtual inline unsigned nodal_index_fvk(const unsigned& i = 0) const
    {
      return i;
    }

    /// Output with default number of plot points
    void output(std::ostream& outfile)
    {
      const unsigned n_plot = 5;
      output(outfile, n_plot);
    }

    /// \short Output FE representation of soln: r,w,sigma_r_r,sigma_phi_phi
    /// at n_plot plot points
    void output(std::ostream& outfile, const unsigned& n_plot);

    /// C_style output with default number of plot points
    void output(FILE* file_pt)
    {
      const unsigned n_plot = 5;
      output(file_pt, n_plot);
    }

    /// \short C-style output FE representation of soln: r,w at
    /// n_plot plot points
    void output(FILE* file_pt, const unsigned& n_plot);

    /// Output exact soln: r,w_exact at n_plot plot points
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

    /// \short Output exact soln: r,w_exact at
    /// n_plot plot points (dummy time-dependent version to
    /// keep intel compiler happy)
    virtual void output_fct(
      std::ostream& outfile,
      const unsigned& n_plot,
      const double& time,
      FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      throw OomphLibError(
        "There is no time-dependent output_fct() for Foeppl von Karman"
        "elements ",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Get error against and norm of exact solution
    void compute_error(std::ostream& outfile,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       double& error,
                       double& norm);

    /// Dummy, time dependent error checker
    void compute_error(std::ostream& outfile,
                       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                       const double& time,
                       double& error,
                       double& norm)
    {
      throw OomphLibError(
        "There is no time-dependent compute_error() for Foeppl von Karman"
        "elements",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Access function: Pointer to pressure function
    AxisymFoepplvonKarmanPressureFctPt& pressure_fct_pt()
    {
      return Pressure_fct_pt;
    }

    /// Access function: Pointer to pressure function. Const version
    AxisymFoepplvonKarmanPressureFctPt pressure_fct_pt() const
    {
      return Pressure_fct_pt;
    }

    /// \short Get pressure term at (Eulerian) position r. This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the strength of the pressure function might be determined by
    /// another system of equations.
    inline virtual void get_pressure_fvk(const unsigned& ipt,
                                         const double& r,
                                         double& pressure) const
    {
      // If no pressure function has been set, return zero
      if (Pressure_fct_pt == 0)
      {
        pressure = 0.0;
      }
      else
      {
        // Get pressure strength
        (*Pressure_fct_pt)(r, pressure);
      }
    }

    /// Get gradient of deflection: gradient[i] = dw/dr_i */
    void get_gradient_of_deflection(const Vector<double>& s,
                                    Vector<double>& gradient) const
    {
      // Find out how many nodes there are in the element
      const unsigned n_node = nnode();

      // Get the index at which the unknown is stored
      const unsigned w_nodal_index = nodal_index_fvk(0);

      // Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidr(n_node, 1);

      // Call the derivatives of the shape and test functions
      dshape_eulerian(s, psi, dpsidr);

      // Initialise to zero
      gradient[0] = 0.0;

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        gradient[0] += this->nodal_value(l, w_nodal_index) * dpsidr(l, 0);
      }
    }

    /// Fill in the residuals with this element's contribution
    void fill_in_contribution_to_residuals(Vector<double>& residuals);

    // hierher Jacobian not yet implemented
    // void fill_in_contribution_to_jacobian(Vector<double> &residuals,
    //                                      DenseMatrix<double> &jacobian);

    /// \short Return FE representation of transverse displacement
    inline double interpolated_w_fvk(const Vector<double>& s) const
    {
      // Find number of nodes
      const unsigned n_node = nnode();

      // Get the index at which the transverse displacement unknown is stored
      const unsigned w_nodal_index = nodal_index_fvk(0);

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

    /// \short Return FE representation of radial displacement
    inline double interpolated_u_fvk(const Vector<double>& s) const
    {
      // Find number of nodes
      const unsigned n_node = nnode();

      // Get the index at which the radial displacement unknown is stored
      const unsigned u_nodal_index = nodal_index_fvk(2);

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise value of u
      double interpolated_u = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += this->nodal_value(l, u_nodal_index) * psi[l];
      }

      return (interpolated_u);
    }

    /// \short Compute in-plane stresses. Return boolean to indicate success
    /// (false if attempt to evaluate stresses at zero radius)
    bool interpolated_stress(const Vector<double>& s,
                             double& sigma_r_r,
                             double& sigma_phi_phi) const;

    /// \short Self-test: Return 0 for OK
    unsigned self_test();

    // switch back on and test!

    /// \short Sets a flag to signify that we are solving the linear,
    /// pure bending equations, and pin all the nodal values that will
    /// not be used in this case
    void use_linear_bending_model()
    {
      // Set the boolean flag
      Linear_bending_model = true;

      // Get the index of the first FvK nodal value
      unsigned first_fvk_nodal_index = nodal_index_fvk();

      // Get the total number of FvK nodal values (assuming they are stored
      // contiguously) at node 0 (it's the same at all nodes anyway)
      unsigned total_fvk_nodal_indices = 3;

      // Get the number of nodes in this element
      unsigned n_node = nnode();

      // Loop over the appropriate nodal indices
      for (unsigned index = first_fvk_nodal_index + 2;
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
    virtual double dshape_and_dtest_eulerian_axisym_fvk(
      const Vector<double>& s,
      Shape& psi,
      DShape& dpsidr,
      Shape& test,
      DShape& dtestdr) const = 0;

    /// \short Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return  Jacobian of mapping
    virtual double dshape_and_dtest_eulerian_at_knot_axisym_fvk(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidr,
      Shape& test,
      DShape& dtestdr) const = 0;

    /// Pointer to FvK parameter
    double* Eta_pt;

    /// Pointer to pressure function:
    AxisymFoepplvonKarmanPressureFctPt Pressure_fct_pt;

    /// Pointer to Poisson's ratio
    double* Nu_pt;

    /// \short Flag which stores whether we are using a linear,
    /// pure bending model instead of the full non-linear Foeppl-von Karman
    bool Linear_bending_model;

  private:
    /// Default value for physical constants
    static double Default_Physical_Constant_Value;
  };

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  //======================================================================
  /// Axisym FoepplvonKarmanElement elements are 1D
  /// Foeppl von Karman elements with isoparametric interpolation for the
  /// function.
  //======================================================================
  template<unsigned NNODE_1D>
  class AxisymFoepplvonKarmanElement :
    public virtual QElement<1, NNODE_1D>,
    public virtual AxisymFoepplvonKarmanEquations
  {
  private:
    /// \short Static int that holds the number of variables at
    /// nodes: always the same
    static const unsigned Initial_Nvalue;

  public:
    ///\short  Constructor: Call constructors for QElement and
    /// AxisymFoepplvonKarmanEquations
    AxisymFoepplvonKarmanElement() :
      QElement<1, NNODE_1D>(), AxisymFoepplvonKarmanEquations()
    {
    }

    /// Broken copy constructor
    AxisymFoepplvonKarmanElement(
      const AxisymFoepplvonKarmanElement<NNODE_1D>& dummy)
    {
      BrokenCopy::broken_copy("AxisymFoepplvonKarmanElement");
    }

    /// Broken assignment operator
    void operator=(const AxisymFoepplvonKarmanElement<NNODE_1D>&)
    {
      BrokenCopy::broken_assign("AxisymFoepplvonKarmanElement");
    }

    /// \short  Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue;
    }

    /// \short Output function:
    ///  r,w,sigma_r_r,sigma_phi_phi
    void output(std::ostream& outfile)
    {
      AxisymFoepplvonKarmanEquations::output(outfile);
    }

    ///  \short Output function:
    ///   r,w,sigma_r_r,sigma_phi_phi at n_plot plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      AxisymFoepplvonKarmanEquations::output(outfile, n_plot);
    }

    /// \short C-style output function:
    ///  r,w
    void output(FILE* file_pt)
    {
      AxisymFoepplvonKarmanEquations::output(file_pt);
    }

    ///  \short C-style output function:
    ///   r,w at n_plot plot points
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      AxisymFoepplvonKarmanEquations::output(file_pt, n_plot);
    }

    /// \short Output function for an exact solution:
    ///  r,w_exact at n_plot plot points
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      AxisymFoepplvonKarmanEquations::output_fct(
        outfile, n_plot, exact_soln_pt);
    }

    /// \short Output function for a time-dependent exact solution.
    ///  r,w_exact at n_plot plot points
    /// (Calls the steady version)
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      AxisymFoepplvonKarmanEquations::output_fct(
        outfile, n_plot, time, exact_soln_pt);
    }

  protected:
    /// \short Shape, test functions & derivs. w.r.t. to global coords.
    /// Return Jacobian.
    inline double dshape_and_dtest_eulerian_axisym_fvk(const Vector<double>& s,
                                                       Shape& psi,
                                                       DShape& dpsidr,
                                                       Shape& test,
                                                       DShape& dtestdr) const;

    /// \short Shape, test functions & derivs. w.r.t. to global coords. at
    /// integration point ipt. Return Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_axisym_fvk(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidr,
      Shape& test,
      DShape& dtestdr) const;
  };

  // Inline functions:

  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned NNODE_1D>
  double AxisymFoepplvonKarmanElement<
    NNODE_1D>::dshape_and_dtest_eulerian_axisym_fvk(const Vector<double>& s,
                                                    Shape& psi,
                                                    DShape& dpsidr,
                                                    Shape& test,
                                                    DShape& dtestdr) const

  {
    // Call the geometrical shape functions and derivatives
    const double J = this->dshape_eulerian(s, psi, dpsidr);

    // Set the test functions equal to the shape functions
    test = psi;
    dtestdr = dpsidr;

    // Return the jacobian
    return J;
  }

  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned NNODE_1D>
  double AxisymFoepplvonKarmanElement<NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_axisym_fvk(const unsigned& ipt,
                                                 Shape& psi,
                                                 DShape& dpsidr,
                                                 Shape& test,
                                                 DShape& dtestdr) const
  {
    // Call the geometrical shape functions and derivatives
    const double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidr);

    // Set the pointers of the test functions
    test = psi;
    dtestdr = dpsidr;

    // Return the jacobian
    return J;
  }

} // namespace oomph

#endif
