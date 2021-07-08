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
// Header file for scalar advection elements

#ifndef OOMPH_SCALAR_ADVECTION_ELEMENTS_HEADER
#define OOMPH_SCALAR_ADVECTION_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include "flux_transport_elements.h"
#include "generic/Qelements.h"
#include "generic/Qspectral_elements.h"
#include "generic/dg_elements.h"

namespace oomph
{
  //==============================================================
  /// Base class for advection equations
  //=============================================================
  template<unsigned DIM>
  class ScalarAdvectionEquations : public FluxTransportEquations<DIM>
  {
    /// \short Typedef for a wind function as a possible function of position
    typedef void (*ScalarAdvectionWindFctPt)(const Vector<double>& x,
                                             Vector<double>& wind);

    /// \short Function pointer to the wind function
    ScalarAdvectionWindFctPt Wind_fct_pt;

  protected:
    /// \short A single flux is interpolated
    inline unsigned nflux() const
    {
      return 1;
    }

    /// \short Return the flux as a function of the unknown
    void flux(const Vector<double>& u, DenseMatrix<double>& f);

    /// \short Return the flux derivatives as a function of the unknowns
    void dflux_du(const Vector<double>& u, RankThreeTensor<double>& df_du);

    ///\short Return the wind at a given position
    inline virtual void get_wind_scalar_adv(const unsigned& ipt,
                                            const Vector<double>& s,
                                            const Vector<double>& x,
                                            Vector<double>& wind) const
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

  public:
    /// Constructor
    ScalarAdvectionEquations() : FluxTransportEquations<DIM>(), Wind_fct_pt(0)
    {
    }

    /// Access function: Pointer to wind function
    ScalarAdvectionWindFctPt& wind_fct_pt()
    {
      return Wind_fct_pt;
    }

    /// Access function: Pointer to wind function. Const version
    ScalarAdvectionWindFctPt wind_fct_pt() const
    {
      return Wind_fct_pt;
    }

    ///\short The number of unknowns at each node is the number of values
    unsigned required_nvalue(const unsigned& n) const
    {
      return 1;
    }

    /// Compute the error and norm of solution integrated over the element
    /// Does not plot the error in the outfile
    void compute_error(
      std::ostream& outfile,
      FiniteElement::UnsteadyExactSolutionFctPt initial_condition_pt,
      const double& t,
      Vector<double>& error,
      Vector<double>& norm)
    {
      // Find the number of fluxes
      const unsigned n_flux = this->nflux();
      // Find the number of nodes
      const unsigned n_node = this->nnode();
      // Storage for the shape function and derivatives of shape function
      Shape psi(n_node);
      DShape dpsidx(n_node, DIM);

      // Find the number of integration points
      unsigned n_intpt = this->integral_pt()->nweight();

      error.resize(n_flux);
      norm.resize(n_flux);
      for (unsigned i = 0; i < n_flux; i++)
      {
        error[i] = 0.0;
        norm[i] = 0.0;
      }

      Vector<double> s(DIM);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the shape functions at the knot
        double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

        // Get the integral weight
        double W = this->integral_pt()->weight(ipt) * J;

        // Storage for the local functions
        Vector<double> interpolated_x(DIM, 0.0);
        Vector<double> interpolated_u(n_flux, 0.0);

        // Loop over the shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Locally cache the shape function
          const double psi_ = psi(l);
          for (unsigned i = 0; i < DIM; i++)
          {
            interpolated_x[i] += this->nodal_position(l, i) * psi_;
          }

          for (unsigned i = 0; i < n_flux; i++)
          {
            // Calculate the velocity and tangent vector
            interpolated_u[i] += this->nodal_value(l, i) * psi_;
          }
        }

        // Get the global wind
        Vector<double> wind(DIM);
        this->get_wind_scalar_adv(ipt, s, interpolated_x, wind);

        // Rescale the position
        for (unsigned i = 0; i < DIM; i++)
        {
          interpolated_x[i] -= wind[i] * t;
        }

        // Now get the initial condition at this value of x
        Vector<double> exact_u(n_flux, 0.0);
        (*initial_condition_pt)(0.0, interpolated_x, exact_u);

        // Loop over the unknowns
        for (unsigned i = 0; i < n_flux; i++)
        {
          error[i] += pow((interpolated_u[i] - exact_u[i]), 2.0) * W;
          norm[i] += interpolated_u[i] * interpolated_u[i] * W;
        }
      }
    }
  };

  template<unsigned DIM, unsigned NNODE_1D>
  class QSpectralScalarAdvectionElement :
    public virtual QSpectralElement<DIM, NNODE_1D>,
    public virtual ScalarAdvectionEquations<DIM>
  {
  public:
    ///\short  Constructor: Call constructors for QElement and
    /// Advection Diffusion equations
    QSpectralScalarAdvectionElement() :
      QSpectralElement<DIM, NNODE_1D>(), ScalarAdvectionEquations<DIM>()
    {
    }

    /// Broken copy constructor
    QSpectralScalarAdvectionElement(
      const QSpectralScalarAdvectionElement<DIM, NNODE_1D>& dummy)
    {
      BrokenCopy::broken_copy("QSpectralScalarAdvectionElement");
    }

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(
     const QSpectralScalarAdvectionElement<DIM,NNODE_1D>&)
     {
      BrokenCopy::broken_assign("QSpectralScalarAdvectionElement");
      }*/

    /// \short  Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return 1;
    }

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream& outfile)
    {
      ScalarAdvectionEquations<DIM>::output(outfile);
    }

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      ScalarAdvectionEquations<DIM>::output(outfile, n_plot);
    }

    /*/// \short C-style output function:
   ///  x,y,u   or    x,y,z,u
   void output(FILE* file_pt)
    {
     ScalarAdvectionEquations<NFLUX,DIM>::output(file_pt);
    }

   ///  \short C-style output function:
   ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
   void output(FILE* file_pt, const unsigned &n_plot)
    {
     ScalarAdvectionEquations<NFLUX,DIM>::output(file_pt,n_plot);
    }

   /// \short Output function for an exact solution:
   ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
   void output_fct(std::ostream &outfile, const unsigned &n_plot,
                   FiniteElement::SteadyExactSolutionFctPt
                   exact_soln_pt)
    {
     ScalarAdvectionEquations<NFLUX,DIM>::
      output_fct(outfile,n_plot,exact_soln_pt);}


   /// \short Output function for a time-dependent exact solution.
   ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
   /// (Calls the steady version)
   void output_fct(std::ostream &outfile, const unsigned &n_plot,
                   const double& time,
                   FiniteElement::UnsteadyExactSolutionFctPt
                   exact_soln_pt)
    {
     ScalarAdvectionEquations<NFLUX,DIM>::
      output_fct(outfile,n_plot,time,exact_soln_pt);
      }*/

  protected:
    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_flux_transport(
      const Vector<double>& s,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const;

    /// \short Shape, test functions & derivs. w.r.t. to global coords. at
    /// integration point ipt. Return Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_flux_transport(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const;
  };

  // Inline functions:

  //======================================================================
  /// \short Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QSpectralScalarAdvectionElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_flux_transport(const Vector<double>& s,
                                             Shape& psi,
                                             DShape& dpsidx,
                                             Shape& test,
                                             DShape& dtestdx) const
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
  double QSpectralScalarAdvectionElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_flux_transport(const unsigned& ipt,
                                                     Shape& psi,
                                                     DShape& dpsidx,
                                                     Shape& test,
                                                     DShape& dtestdx) const
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
  /// \short Face geometry for the QScalarAdvectionElement elements:
  /// The spatial dimension of the face elements is one lower than that
  /// of the bulk element but they have the same number of points along
  /// their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QSpectralScalarAdvectionElement<DIM, NNODE_1D>> :
    public virtual QSpectralElement<DIM - 1, NNODE_1D>
  {
  public:
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : QSpectralElement<DIM - 1, NNODE_1D>() {}
  };

  //======================================================================
  /// FaceElement for Discontinuous Galerkin Problems
  //======================================================================
  template<class ELEMENT>
  class DGScalarAdvectionFaceElement :
    public virtual FaceGeometry<ELEMENT>,
    public virtual DGFaceElement
  {
  public:
    /// Constructor
    DGScalarAdvectionFaceElement(FiniteElement* const& element_pt,
                                 const int& face_index) :
      FaceGeometry<ELEMENT>(), DGFaceElement()
    {
      // Attach geometric information to the element
      // N.B. This function also assigns nbulk_value from required_nvalue
      // of the bulk element.
      element_pt->build_face_element(face_index, this);
    }

    // There is a single required n_flux
    unsigned required_nflux()
    {
      return 1;
    }

    /// Calculate the normal flux, so this is the dot product of the
    /// numerical flux with n_in
    void numerical_flux(const Vector<double>& n_out,
                        const Vector<double>& u_int,
                        const Vector<double>& u_ext,
                        Vector<double>& flux)
    {
      const unsigned dim = this->nodal_dimension();
      Vector<double> Wind(dim);
      Vector<double> s, x;
      // Dummy integration point
      unsigned ipt = 0;
      dynamic_cast<ELEMENT*>(this->bulk_element_pt())
        ->get_wind_scalar_adv(ipt, s, x, Wind);

      // Now we can work this out for standard upwind advection
      double dot = 0.0;
      for (unsigned i = 0; i < dim; i++)
      {
        dot += Wind[i] * n_out[i];
      }

      const unsigned n_value = this->required_nflux();
      if (dot >= 0.0)
      {
        for (unsigned n = 0; n < n_value; n++)
        {
          flux[n] = dot * u_int[n];
        }
      }
      else
      {
        for (unsigned n = 0; n < n_value; n++)
        {
          flux[n] = dot * u_ext[n];
        }
      }

      // flux[0] = 0.5*(u_int[0]+u_ext[0])*n_out[0];
    }
  };

  //=================================================================
  /// General DGScalarAdvectionClass. Establish the template parameters
  //===================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class DGSpectralScalarAdvectionElement
  {
  };

  //==================================================================
  // Specialization for 1D DG Advection element
  //==================================================================
  template<unsigned NNODE_1D>
  class DGSpectralScalarAdvectionElement<1, NNODE_1D> :
    public QSpectralScalarAdvectionElement<1, NNODE_1D>,
    public DGElement
  {
    friend class DGScalarAdvectionFaceElement<
      DGSpectralScalarAdvectionElement<1, NNODE_1D>>;

    Vector<double> Inverse_mass_diagonal;

  public:
    // There is a single required n_flux
    unsigned required_nflux()
    {
      return 1;
    }

    // Calculate averages
    void calculate_element_averages(double*& average_value)
    {
      FluxTransportEquations<1>::calculate_element_averages(average_value);
    }

    // Constructor
    DGSpectralScalarAdvectionElement() :
      QSpectralScalarAdvectionElement<1, NNODE_1D>(), DGElement()
    {
    }

    ~DGSpectralScalarAdvectionElement() {}

    void build_all_faces()
    {
      // Make the two faces
      Face_element_pt.resize(2);
      // Make the face on the left
      Face_element_pt[0] = new DGScalarAdvectionFaceElement<
        DGSpectralScalarAdvectionElement<1, NNODE_1D>>(this, -1);
      // Make the face on the right
      Face_element_pt[1] = new DGScalarAdvectionFaceElement<
        DGSpectralScalarAdvectionElement<1, NNODE_1D>>(this, +1);
    }

    ///\short Compute the residuals for the Navier--Stokes equations;
    /// flag=1(or 0): do (or don't) compute the Jacobian as well.
    void fill_in_generic_residual_contribution_flux_transport(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
    {
      QSpectralScalarAdvectionElement<1, NNODE_1D>::
        fill_in_generic_residual_contribution_flux_transport(
          residuals, jacobian, mass_matrix, flag);

      this->add_flux_contributions_to_residuals(residuals, jacobian, flag);
    }

    //============================================================================
    /// Function that returns the current value of the residuals
    /// multiplied by the inverse mass matrix (virtual so that it can be
    /// overloaded specific elements in which time and memory saving tricks can
    /// be applied)
    //============================================================================
    void get_inverse_mass_matrix_times_residuals(Vector<double>& minv_res)
    {
      // If there are external data this is not going to work
      if (nexternal_data() > 0)
      {
        std::ostringstream error_stream;
        error_stream
          << "Cannot use a discontinuous formulation for the mass matrix when\n"
          << "there are external data.\n "
          << "Do not call Problem::enable_discontinuous_formulation()\n";

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }

      // Now let's assemble stuff
      const unsigned n_dof = this->ndof();

      // Resize and initialise the vector that will holds the residuals
      minv_res.resize(n_dof);
      for (unsigned n = 0; n < n_dof; n++)
      {
        minv_res[n] = 0.0;
      }

      // If we are recycling the mass matrix
      if (Mass_matrix_reuse_is_enabled && Mass_matrix_has_been_computed)
      {
        // Get the residuals
        this->fill_in_contribution_to_residuals(minv_res);
      }
      // Otherwise
      else
      {
        // Temporary mass matrix
        DenseDoubleMatrix M(n_dof, n_dof, 0.0);

        // Get the local mass matrix and residuals
        this->fill_in_contribution_to_mass_matrix(minv_res, M);

        // Store the diagonal entries
        Inverse_mass_diagonal.clear();
        for (unsigned n = 0; n < n_dof; n++)
        {
          Inverse_mass_diagonal.push_back(1.0 / M(n, n));
        }

        // The mass matrix has been computed
        Mass_matrix_has_been_computed = true;
      }

      for (unsigned n = 0; n < n_dof; n++)
      {
        minv_res[n] *= Inverse_mass_diagonal[n];
      }
    }
  };

  //=======================================================================
  /// Face geometry of the 1D  DG elements
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<DGSpectralScalarAdvectionElement<1, NNODE_1D>> :
    public virtual PointElement
  {
  public:
    FaceGeometry() : PointElement() {}
  };

  //==================================================================
  /// Specialisation for 2D DG Elements
  //==================================================================
  template<unsigned NNODE_1D>
  class DGSpectralScalarAdvectionElement<2, NNODE_1D> :
    public QSpectralScalarAdvectionElement<2, NNODE_1D>,
    public DGElement
  {
    friend class DGScalarAdvectionFaceElement<
      DGSpectralScalarAdvectionElement<2, NNODE_1D>>;

  public:
    // Calculate averages
    void calculate_element_averages(double*& average_value)
    {
      FluxTransportEquations<2>::calculate_element_averages(average_value);
    }

    // There is a single required n_flux
    unsigned required_nflux()
    {
      return 1;
    }

    // Constructor
    DGSpectralScalarAdvectionElement() :
      QSpectralScalarAdvectionElement<2, NNODE_1D>(), DGElement()
    {
    }

    ~DGSpectralScalarAdvectionElement() {}

    void build_all_faces()
    {
      Face_element_pt.resize(4);
      Face_element_pt[0] = new DGScalarAdvectionFaceElement<
        DGSpectralScalarAdvectionElement<2, NNODE_1D>>(this, 2);
      Face_element_pt[1] = new DGScalarAdvectionFaceElement<
        DGSpectralScalarAdvectionElement<2, NNODE_1D>>(this, 1);
      Face_element_pt[2] = new DGScalarAdvectionFaceElement<
        DGSpectralScalarAdvectionElement<2, NNODE_1D>>(this, -2);
      Face_element_pt[3] = new DGScalarAdvectionFaceElement<
        DGSpectralScalarAdvectionElement<2, NNODE_1D>>(this, -1);
    }

    ///\short Compute the residuals for the Navier--Stokes equations;
    /// flag=1(or 0): do (or don't) compute the Jacobian as well.
    void fill_in_generic_residual_contribution_flux_transport(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
    {
      QSpectralScalarAdvectionElement<2, NNODE_1D>::
        fill_in_generic_residual_contribution_flux_transport(
          residuals, jacobian, mass_matrix, flag);

      this->add_flux_contributions_to_residuals(residuals, jacobian, flag);
    }
  };

  //=======================================================================
  /// Face geometry of the DG elements
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<DGSpectralScalarAdvectionElement<2, NNODE_1D>> :
    public virtual QSpectralElement<1, NNODE_1D>
  {
  public:
    FaceGeometry() : QSpectralElement<1, NNODE_1D>() {}
  };

  //=============================================================
  /// Non-spectral version of the classes
  //============================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class QScalarAdvectionElement :
    public virtual QElement<DIM, NNODE_1D>,
    public virtual ScalarAdvectionEquations<DIM>
  {
  public:
    ///\short  Constructor: Call constructors for QElement and
    /// Advection Diffusion equations
    QScalarAdvectionElement() :
      QElement<DIM, NNODE_1D>(), ScalarAdvectionEquations<DIM>()
    {
    }

    /// Broken copy constructor
    QScalarAdvectionElement(const QScalarAdvectionElement<DIM, NNODE_1D>& dummy)
    {
      BrokenCopy::broken_copy("QScalarAdvectionElement");
    }

    /// Broken assignment operator
    /*void operator=(
     const QScalarAdvectionElement<DIM,NNODE_1D>&)
     {
      BrokenCopy::broken_assign("QScalarAdvectionElement");
      }*/

    /// \short  Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return 1;
    }

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream& outfile)
    {
      ScalarAdvectionEquations<DIM>::output(outfile);
    }

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      ScalarAdvectionEquations<DIM>::output(outfile, n_plot);
    }

    /*/// \short C-style output function:
   ///  x,y,u   or    x,y,z,u
   void output(FILE* file_pt)
    {
     ScalarAdvectionEquations<NFLUX,DIM>::output(file_pt);
    }

   ///  \short C-style output function:
   ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
   void output(FILE* file_pt, const unsigned &n_plot)
    {
     ScalarAdvectionEquations<NFLUX,DIM>::output(file_pt,n_plot);
    }

   /// \short Output function for an exact solution:
   ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
   void output_fct(std::ostream &outfile, const unsigned &n_plot,
                   FiniteElement::SteadyExactSolutionFctPt
                   exact_soln_pt)
    {
     ScalarAdvectionEquations<NFLUX,DIM>::
      output_fct(outfile,n_plot,exact_soln_pt);}


   /// \short Output function for a time-dependent exact solution.
   ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
   /// (Calls the steady version)
   void output_fct(std::ostream &outfile, const unsigned &n_plot,
                   const double& time,
                   FiniteElement::UnsteadyExactSolutionFctPt
                   exact_soln_pt)
    {
     ScalarAdvectionEquations<NFLUX,DIM>::
      output_fct(outfile,n_plot,time,exact_soln_pt);
      }*/

  protected:
    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_flux_transport(
      const Vector<double>& s,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const;

    /// \short Shape, test functions & derivs. w.r.t. to global coords. at
    /// integration point ipt. Return Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_flux_transport(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const;
  };

  // Inline functions:

  //======================================================================
  /// \short Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QScalarAdvectionElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_flux_transport(const Vector<double>& s,
                                             Shape& psi,
                                             DShape& dpsidx,
                                             Shape& test,
                                             DShape& dtestdx) const
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
  double QScalarAdvectionElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_flux_transport(const unsigned& ipt,
                                                     Shape& psi,
                                                     DShape& dpsidx,
                                                     Shape& test,
                                                     DShape& dtestdx) const
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
  /// \short Face geometry for the QScalarAdvectionElement elements:
  /// The spatial dimension of the face elements is one lower than that
  /// of the bulk element but they have the same number of points along
  /// their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QScalarAdvectionElement<DIM, NNODE_1D>> :
    public virtual QElement<DIM - 1, NNODE_1D>
  {
  public:
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : QElement<DIM - 1, NNODE_1D>() {}
  };

  //=================================================================
  /// General DGScalarAdvectionClass. Establish the template parameters
  //===================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class DGScalarAdvectionElement
  {
  };

  //==================================================================
  // Specialization for 1D DG Advection element
  //==================================================================
  template<unsigned NNODE_1D>
  class DGScalarAdvectionElement<1, NNODE_1D> :
    public QScalarAdvectionElement<1, NNODE_1D>,
    public DGElement
  {
    friend class DGScalarAdvectionFaceElement<
      DGScalarAdvectionElement<1, NNODE_1D>>;

  public:
    // There is a single required n_flux
    unsigned required_nflux()
    {
      return 1;
    }

    // Calculate averages
    void calculate_element_averages(double*& average_value)
    {
      FluxTransportEquations<1>::calculate_element_averages(average_value);
    }

    // Constructor
    DGScalarAdvectionElement() :
      QScalarAdvectionElement<1, NNODE_1D>(), DGElement()
    {
    }

    ~DGScalarAdvectionElement() {}

    void build_all_faces()
    {
      // Make the two faces
      Face_element_pt.resize(2);
      // Make the face on the left
      Face_element_pt[0] =
        new DGScalarAdvectionFaceElement<DGScalarAdvectionElement<1, NNODE_1D>>(
          this, -1);
      // Make the face on the right
      Face_element_pt[1] =
        new DGScalarAdvectionFaceElement<DGScalarAdvectionElement<1, NNODE_1D>>(
          this, +1);
    }

    ///\short Compute the residuals for the Navier--Stokes equations;
    /// flag=1(or 0): do (or don't) compute the Jacobian as well.
    void fill_in_generic_residual_contribution_flux_transport(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
    {
      QScalarAdvectionElement<1, NNODE_1D>::
        fill_in_generic_residual_contribution_flux_transport(
          residuals, jacobian, mass_matrix, flag);

      this->add_flux_contributions_to_residuals(residuals, jacobian, flag);
    }
  };

  //=======================================================================
  /// Face geometry of the 1D  DG elements
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<DGScalarAdvectionElement<1, NNODE_1D>> :
    public virtual PointElement
  {
  public:
    FaceGeometry() : PointElement() {}
  };

  //==================================================================
  /// Specialisation for 2D DG Elements
  //==================================================================
  template<unsigned NNODE_1D>
  class DGScalarAdvectionElement<2, NNODE_1D> :
    public QScalarAdvectionElement<2, NNODE_1D>,
    public DGElement
  {
    friend class DGScalarAdvectionFaceElement<
      DGScalarAdvectionElement<2, NNODE_1D>>;

  public:
    // There is a single required n_flux
    unsigned required_nflux()
    {
      return 1;
    }

    // Calculate averages
    void calculate_element_averages(double*& average_value)
    {
      FluxTransportEquations<2>::calculate_element_averages(average_value);
    }

    // Constructor
    DGScalarAdvectionElement() :
      QScalarAdvectionElement<2, NNODE_1D>(), DGElement()
    {
    }

    ~DGScalarAdvectionElement() {}

    void build_all_faces()
    {
      Face_element_pt.resize(4);
      Face_element_pt[0] =
        new DGScalarAdvectionFaceElement<DGScalarAdvectionElement<2, NNODE_1D>>(
          this, 2);
      Face_element_pt[1] =
        new DGScalarAdvectionFaceElement<DGScalarAdvectionElement<2, NNODE_1D>>(
          this, 1);
      Face_element_pt[2] =
        new DGScalarAdvectionFaceElement<DGScalarAdvectionElement<2, NNODE_1D>>(
          this, -2);
      Face_element_pt[3] =
        new DGScalarAdvectionFaceElement<DGScalarAdvectionElement<2, NNODE_1D>>(
          this, -1);
    }

    ///\short Compute the residuals for the Navier--Stokes equations;
    /// flag=1(or 0): do (or don't) compute the Jacobian as well.
    void fill_in_generic_residual_contribution_flux_transport(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
    {
      QScalarAdvectionElement<2, NNODE_1D>::
        fill_in_generic_residual_contribution_flux_transport(
          residuals, jacobian, mass_matrix, flag);

      this->add_flux_contributions_to_residuals(residuals, jacobian, flag);
    }
  };

  //=======================================================================
  /// Face geometry of the DG elements
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<DGScalarAdvectionElement<2, NNODE_1D>> :
    public virtual QElement<1, NNODE_1D>
  {
  public:
    FaceGeometry() : QElement<1, NNODE_1D>() {}
  };

} // namespace oomph

#endif
