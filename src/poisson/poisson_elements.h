// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
// Header file for Poisson elements
#ifndef OOMPH_POISSON_ELEMENTS_HEADER
#define OOMPH_POISSON_ELEMENTS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include <sstream>

// OOMPH-LIB headers
#include "generic/projection.h"
#include "generic/nodes.h"
#include "generic/Qelements.h"
#include "generic/oomph_utilities.h"


namespace oomph
{
  //=============================================================
  /// A class for all isoparametric elements that solve the
  /// Poisson equations.
  /// \f[
  /// \frac{\partial^2 u}{\partial x_i^2} = f(x_j)
  /// \f]
  /// This contains the generic maths. Shape functions, geometric
  /// mapping etc. must get implemented in derived class.
  //=============================================================
  template<unsigned DIM>
  class PoissonEquations : public virtual FiniteElement
  {
  public:
    /// \short Function pointer to source function fct(x,f(x)) --
    /// x is a Vector!
    typedef void (*PoissonSourceFctPt)(const Vector<double>& x, double& f);


    /// \short Function pointer to gradient of source function  fct(x,g(x)) --
    /// x is a Vector!
    typedef void (*PoissonSourceFctGradientPt)(const Vector<double>& x,
                                               Vector<double>& gradient);


    /// Constructor (must initialise the Source_fct_pt to null)
    PoissonEquations() : Source_fct_pt(0), Source_fct_gradient_pt(0) {}

    /// Broken copy constructor
    PoissonEquations(const PoissonEquations& dummy)
    {
      BrokenCopy::broken_copy("PoissonEquations");
    }

    /// Broken assignment operator
    void operator=(const PoissonEquations&)
    {
      BrokenCopy::broken_assign("PoissonEquations");
    }

    /// \short Return the index at which the unknown value
    /// is stored. The default value, 0, is appropriate for single-physics
    /// problems, when there is only one variable, the value that satisfies
    /// the poisson equation.
    /// In derived multi-physics elements, this function should be overloaded
    /// to reflect the chosen storage scheme. Note that these equations require
    /// that the unknown is always stored at the same index at each node.
    virtual inline unsigned u_index_poisson() const
    {
      return 0;
    }


    /// \short Output solution in data vector at local cordinates s:
    /// x,y [,z], u
    void point_output_data(const Vector<double>& s, Vector<double>& data)
    {
      // Dimension
      unsigned dim = s.size();

      // Resize data for values
      data.resize(dim + 1);

      // Write values in the vector
      for (unsigned i = 0; i < dim; i++)
      {
        data[i] = interpolated_x(s, i);
      }
      data[dim] = this->interpolated_u_poisson(s);
    }


    /// \short Number of scalars/fields output by this element. Reimplements
    /// broken virtual function in base class.
    unsigned nscalar_paraview() const
    {
      return 1;
    }

    /// \short Write values of the i-th scalar field at the plot points. Needs
    /// to be implemented for each new specific element type.
    void scalar_value_paraview(std::ofstream& file_out,
                               const unsigned& i,
                               const unsigned& nplot) const
    {
#ifdef PARANOID
      if (i != 0)
      {
        std::stringstream error_stream;
        error_stream
          << "Poisson elements only store a single field so i must be 0 rather"
          << " than " << i << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      unsigned local_loop = this->nplot_points_paraview(nplot);
      for (unsigned j = 0; j < local_loop; j++)
      {
        // Get the local coordinate of the required plot point
        Vector<double> s(DIM);
        this->get_s_plot(j, nplot, s);

        file_out << this->interpolated_u_poisson(s) << std::endl;
      }
    }

    /// \short Name of the i-th scalar field. Default implementation
    /// returns V1 for the first one, V2 for the second etc. Can (should!) be
    /// overloaded with more meaningful names in specific elements.
    std::string scalar_name_paraview(const unsigned& i) const
    {
#ifdef PARANOID
      if (i != 0)
      {
        std::stringstream error_stream;
        error_stream
          << "Poisson elements only store a single field so i must be 0 rather"
          << " than " << i << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      return "Poisson solution";
    }


    /// \short Write values of the i-th scalar field at the plot points. Needs
    /// to be implemented for each new specific element type.
    void scalar_value_fct_paraview(
      std::ofstream& file_out,
      const unsigned& i,
      const unsigned& nplot,
      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt) const
    {
#ifdef PARANOID
      if (i != 0)
      {
        std::stringstream error_stream;
        error_stream << "Poisson equation has only one field. Can't call "
                     << "this function for value " << i << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Vector of local coordinates
      Vector<double> s(DIM);

      // Vector for coordinates
      Vector<double> x(DIM);

      // Exact solution Vector
      Vector<double> exact_soln(1);

      // Loop over plot points
      unsigned num_plot_points = nplot_points_paraview(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        get_s_plot(iplot, nplot, s);

        // Get x position as Vector
        interpolated_x(s, x);

        // Get exact solution at this point
        (*exact_soln_pt)(x, exact_soln);

        // Write it
        file_out << exact_soln[0] << std::endl;
      }
    }

    /// Output with default number of plot points
    void output(std::ostream& outfile)
    {
      const unsigned n_plot = 5;
      output(outfile, n_plot);
    }

    /// \short Output FE representation of soln: x,y,u or x,y,z,u at
    /// n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& n_plot);

    /// C_style output with default number of plot points
    void output(FILE* file_pt)
    {
      const unsigned n_plot = 5;
      output(file_pt, n_plot);
    }

    /// \short C-style output FE representation of soln: x,y,u or x,y,z,u at
    /// n_plot^DIM plot points
    void output(FILE* file_pt, const unsigned& n_plot);

    /// Output exact soln: x,y,u_exact or x,y,z,u_exact at n_plot^DIM plot
    /// points
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

    /// \short Output exact soln: x,y,u_exact or x,y,z,u_exact at
    /// n_plot^DIM plot points (dummy time-dependent version to
    /// keep intel compiler happy)
    virtual void output_fct(
      std::ostream& outfile,
      const unsigned& n_plot,
      const double& time,
      FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      throw OomphLibError(
        "There is no time-dependent output_fct() for Poisson elements ",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }


    /// \short Compute norm of solution: square of the L2 norm
    void compute_norm(double& norm);

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
        "There is no time-dependent compute_error() for Poisson elements",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Access function: Pointer to source function
    PoissonSourceFctPt& source_fct_pt()
    {
      return Source_fct_pt;
    }

    /// Access function: Pointer to source function. Const version
    PoissonSourceFctPt source_fct_pt() const
    {
      return Source_fct_pt;
    }

    /// Access function: Pointer to gradient of source function
    PoissonSourceFctGradientPt& source_fct_gradient_pt()
    {
      return Source_fct_gradient_pt;
    }

    /// Access function: Pointer to gradient source function. Const version
    PoissonSourceFctGradientPt source_fct_gradient_pt() const
    {
      return Source_fct_gradient_pt;
    }


    /// Get source term at (Eulerian) position x. This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the strength of the source function might be determined by
    /// another system of equations.
    inline virtual void get_source_poisson(const unsigned& ipt,
                                           const Vector<double>& x,
                                           double& source) const
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


    /// Get gradient of source term at (Eulerian) position x. This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the strength of the source function might be determined by
    /// another system of equations. Computed via function pointer
    /// (if set) or by finite differencing (default)
    inline virtual void get_source_gradient_poisson(
      const unsigned& ipt,
      const Vector<double>& x,
      Vector<double>& gradient) const
    {
      // If no gradient function has been set, FD it
      if (Source_fct_gradient_pt == 0)
      {
        // Reference value
        double source = 0.0;
        get_source_poisson(ipt, x, source);

        // FD it
        double eps_fd = GeneralisedElement::Default_fd_jacobian_step;
        double source_pls = 0.0;
        Vector<double> x_pls(x);
        for (unsigned i = 0; i < DIM; i++)
        {
          x_pls[i] += eps_fd;
          get_source_poisson(ipt, x_pls, source_pls);
          gradient[i] = (source_pls - source) / eps_fd;
          x_pls[i] = x[i];
        }
      }
      else
      {
        // Get gradient
        (*Source_fct_gradient_pt)(x, gradient);
      }
    }


    /// Get flux: flux[i] = du/dx_i
    void get_flux(const Vector<double>& s, Vector<double>& flux) const
    {
      // Find out how many nodes there are in the element
      const unsigned n_node = nnode();

      // Get the index at which the unknown is stored
      const unsigned u_nodal_index = u_index_poisson();

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
          flux[j] += this->nodal_value(l, u_nodal_index) * dpsidx(l, j);
        }
      }
    }


    /// Get derivative of flux w.r.t. to nodal values:
    /// dflux_dnodal_u[i][j] = d ( du/dx_i ) / dU_j
    void get_dflux_dnodal_u(const Vector<double>& s,
                            Vector<Vector<double>>& dflux_dnodal_u) const
    {
      // Find out how many nodes there are in the element
      const unsigned n_node = nnode();

      // Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidx(n_node, DIM);

      // Call the derivatives of the shape and test functions
      dshape_eulerian(s, psi, dpsidx);

      // And here it is...
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = 0; j < n_node; j++)
        {
          dflux_dnodal_u[i][j] = dpsidx(j, i);
        }
      }
    }


    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      fill_in_generic_residual_contribution_poisson(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }


    /// Add the element's contribution to its residual vector and
    /// element Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_poisson(residuals, jacobian, 1);
    }


    /// \short Return FE representation of function value u_poisson(s)
    /// at local coordinate s
    virtual inline double interpolated_u_poisson(const Vector<double>& s) const
    {
      // Find number of nodes
      const unsigned n_node = nnode();

      // Get the index at which the poisson unknown is stored
      const unsigned u_nodal_index = u_index_poisson();

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


    /// \short Compute derivatives of elemental residual vector with respect
    /// to nodal coordinates. Overwrites default implementation in
    /// FiniteElement base class.
    /// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
    virtual void get_dresidual_dnodal_coordinates(
      RankThreeTensor<double>& dresidual_dnodal_coordinates);

    /// \short Self-test: Return 0 for OK
    unsigned self_test();


  protected:
    /// \short Shape/test functions and derivs w.r.t. to global coords at
    /// local coord. s; return  Jacobian of mapping
    virtual double dshape_and_dtest_eulerian_poisson(const Vector<double>& s,
                                                     Shape& psi,
                                                     DShape& dpsidx,
                                                     Shape& test,
                                                     DShape& dtestdx) const = 0;


    /// \short Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return  Jacobian of mapping
    virtual double dshape_and_dtest_eulerian_at_knot_poisson(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const = 0;

    /// \short Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return Jacobian of mapping (J). Also compute
    /// derivatives of dpsidx, dtestdx and J w.r.t. nodal coordinates.
    virtual double dshape_and_dtest_eulerian_at_knot_poisson(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      RankFourTensor<double>& d_dpsidx_dX,
      Shape& test,
      DShape& dtestdx,
      RankFourTensor<double>& d_dtestdx_dX,
      DenseMatrix<double>& djacobian_dX) const = 0;

    /// \short Compute element residual Vector only (if flag=and/or element
    /// Jacobian matrix
    virtual void fill_in_generic_residual_contribution_poisson(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag);

    /// Pointer to source function:
    PoissonSourceFctPt Source_fct_pt;

    /// Pointer to gradient of source function
    PoissonSourceFctGradientPt Source_fct_gradient_pt;
  };


  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// QPoissonElement elements are linear/quadrilateral/brick-shaped
  /// Poisson elements with isoparametric interpolation for the function.
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class QPoissonElement : public virtual QElement<DIM, NNODE_1D>,
                          public virtual PoissonEquations<DIM>
  {
  private:
    /// \short Static int that holds the number of variables at
    /// nodes: always the same
    static const unsigned Initial_Nvalue;

  public:
    ///\short  Constructor: Call constructors for QElement and
    /// Poisson equations
    QPoissonElement() : QElement<DIM, NNODE_1D>(), PoissonEquations<DIM>() {}

    /// Broken copy constructor
    QPoissonElement(const QPoissonElement<DIM, NNODE_1D>& dummy)
    {
      BrokenCopy::broken_copy("QPoissonElement");
    }

    /// Broken assignment operator
    void operator=(const QPoissonElement<DIM, NNODE_1D>&)
    {
      BrokenCopy::broken_assign("QPoissonElement");
    }

    /// \short  Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue;
    }

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream& outfile)
    {
      PoissonEquations<DIM>::output(outfile);
    }


    ///  \short Output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      PoissonEquations<DIM>::output(outfile, n_plot);
    }


    /// \short C-style output function:
    ///  x,y,u   or    x,y,z,u
    void output(FILE* file_pt)
    {
      PoissonEquations<DIM>::output(file_pt);
    }


    ///  \short C-style output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      PoissonEquations<DIM>::output(file_pt, n_plot);
    }


    /// \short Output function for an exact solution:
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      PoissonEquations<DIM>::output_fct(outfile, n_plot, exact_soln_pt);
    }


    /// \short Output function for a time-dependent exact solution.
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
    /// (Calls the steady version)
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      PoissonEquations<DIM>::output_fct(outfile, n_plot, time, exact_soln_pt);
    }


  protected:
    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_poisson(const Vector<double>& s,
                                                    Shape& psi,
                                                    DShape& dpsidx,
                                                    Shape& test,
                                                    DShape& dtestdx) const;


    /// \short Shape, test functions & derivs. w.r.t. to global coords. at
    /// integration point ipt. Return Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_poisson(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const;

    /// \short Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return Jacobian of mapping (J). Also compute
    /// derivatives of dpsidx, dtestdx and J w.r.t. nodal coordinates.
    inline double dshape_and_dtest_eulerian_at_knot_poisson(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      RankFourTensor<double>& d_dpsidx_dX,
      Shape& test,
      DShape& dtestdx,
      RankFourTensor<double>& d_dtestdx_dX,
      DenseMatrix<double>& djacobian_dX) const;
  };


  // Inline functions:


  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QPoissonElement<DIM, NNODE_1D>::dshape_and_dtest_eulerian_poisson(
    const Vector<double>& s,
    Shape& psi,
    DShape& dpsidx,
    Shape& test,
    DShape& dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    const double J = this->dshape_eulerian(s, psi, dpsidx);

    // Set the test functions equal to the shape functions
    test = psi;
    dtestdx = dpsidx;

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
  double QPoissonElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_poisson(const unsigned& ipt,
                                              Shape& psi,
                                              DShape& dpsidx,
                                              Shape& test,
                                              DShape& dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    const double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

    // Set the pointers of the test functions
    test = psi;
    dtestdx = dpsidx;

    // Return the jacobian
    return J;
  }


  //======================================================================
  /// Define the shape functions (psi) and test functions (test) and
  /// their derivatives w.r.t. global coordinates (dpsidx and dtestdx)
  /// and return Jacobian of mapping (J). Additionally compute the
  /// derivatives of dpsidx, dtestdx and J w.r.t. nodal coordinates.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QPoissonElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_poisson(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      RankFourTensor<double>& d_dpsidx_dX,
      Shape& test,
      DShape& dtestdx,
      RankFourTensor<double>& d_dtestdx_dX,
      DenseMatrix<double>& djacobian_dX) const
  {
    // Call the geometrical shape functions and derivatives
    const double J = this->dshape_eulerian_at_knot(
      ipt, psi, dpsidx, djacobian_dX, d_dpsidx_dX);

    // Set the pointers of the test functions
    test = psi;
    dtestdx = dpsidx;
    d_dtestdx_dX = d_dpsidx_dX;

    // Return the jacobian
    return J;
  }


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the QPoissonElement elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QPoissonElement<DIM, NNODE_1D>>
    : public virtual QElement<DIM - 1, NNODE_1D>
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
  /// Face geometry for the 1D QPoissonElement elements: Point elements
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<QPoissonElement<1, NNODE_1D>> : public virtual PointElement
  {
  public:
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : PointElement() {}
  };


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //==========================================================
  /// Poisson upgraded to become projectable
  //==========================================================
  template<class POISSON_ELEMENT>
  class ProjectablePoissonElement
    : public virtual ProjectableElement<POISSON_ELEMENT>
  {
  public:
    /// \short Specify the values associated with field fld.
    /// The information is returned in a vector of pairs which comprise
    /// the Data object and the value within it, that correspond to field fld.
    Vector<std::pair<Data*, unsigned>> data_values_of_field(const unsigned& fld)
    {
#ifdef PARANOID
      if (fld != 0)
      {
        std::stringstream error_stream;
        error_stream << "Poisson elements only store a single field so fld "
                        "must be 0 rather"
                     << " than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Create the vector
      unsigned nnod = this->nnode();
      Vector<std::pair<Data*, unsigned>> data_values(nnod);

      // Loop over all nodes
      for (unsigned j = 0; j < nnod; j++)
      {
        // Add the data value associated field: The node itself
        data_values[j] = std::make_pair(this->node_pt(j), fld);
      }

      // Return the vector
      return data_values;
    }

    /// \short Number of fields to be projected: Just one
    unsigned nfields_for_projection()
    {
      return 1;
    }

    /// \short Number of history values to be stored for fld-th field
    /// (includes current value!)
    unsigned nhistory_values_for_projection(const unsigned& fld)
    {
#ifdef PARANOID
      if (fld != 0)
      {
        std::stringstream error_stream;
        error_stream << "Poisson elements only store a single field so fld "
                        "must be 0 rather"
                     << " than " << fld << std::endl;
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
    double jacobian_and_shape_of_field(const unsigned& fld,
                                       const Vector<double>& s,
                                       Shape& psi)
    {
#ifdef PARANOID
      if (fld != 0)
      {
        std::stringstream error_stream;
        error_stream << "Poisson elements only store a single field so fld "
                        "must be 0 rather"
                     << " than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      unsigned n_dim = this->dim();
      unsigned n_node = this->nnode();
      Shape test(n_node);
      DShape dpsidx(n_node, n_dim), dtestdx(n_node, n_dim);
      double J =
        this->dshape_and_dtest_eulerian_poisson(s, psi, dpsidx, test, dtestdx);
      return J;
    }


    /// \short Return interpolated field fld at local coordinate s, at time
    /// level t (t=0: present; t>0: history values)
    double get_field(const unsigned& t,
                     const unsigned& fld,
                     const Vector<double>& s)
    {
#ifdef PARANOID
      if (fld != 0)
      {
        std::stringstream error_stream;
        error_stream << "Poisson elements only store a single field so fld "
                        "must be 0 rather"
                     << " than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // Find the index at which the variable is stored
      unsigned u_nodal_index = this->u_index_poisson();

      // Local shape function
      unsigned n_node = this->nnode();
      Shape psi(n_node);

      // Find values of shape function
      this->shape(s, psi);

      // Initialise value of u
      double interpolated_u = 0.0;

      // Sum over the local nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += this->nodal_value(l, u_nodal_index) * psi[l];
      }
      return interpolated_u;
    }


    /// Return number of values in field fld: One per node
    unsigned nvalue_of_field(const unsigned& fld)
    {
#ifdef PARANOID
      if (fld != 0)
      {
        std::stringstream error_stream;
        error_stream << "Poisson elements only store a single field so fld "
                        "must be 0 rather"
                     << " than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return this->nnode();
    }


    /// Return local equation number of value j in field fld.
    int local_equation(const unsigned& fld, const unsigned& j)
    {
#ifdef PARANOID
      if (fld != 0)
      {
        std::stringstream error_stream;
        error_stream << "Poisson elements only store a single field so fld "
                        "must be 0 rather"
                     << " than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      const unsigned u_nodal_index = this->u_index_poisson();
      return this->nodal_local_eqn(j, u_nodal_index);
    }
  };


  //=======================================================================
  /// Face geometry for element is the same as that for the underlying
  /// wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<ProjectablePoissonElement<ELEMENT>>
    : public virtual FaceGeometry<ELEMENT>
  {
  public:
    FaceGeometry() : FaceGeometry<ELEMENT>() {}
  };


  //=======================================================================
  /// Face geometry of the Face Geometry for element is the same as
  /// that for the underlying wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<FaceGeometry<ProjectablePoissonElement<ELEMENT>>>
    : public virtual FaceGeometry<FaceGeometry<ELEMENT>>
  {
  public:
    FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT>>() {}
  };


} // namespace oomph

#endif
